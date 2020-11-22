#include "graphics.h"

Matrix4f cameraLookAt(Vec3f cameraPosition, // origin in old coordinates
		      Vec3f cameraDirection, // z unit vector in old coordinates
		      Vec3f cameraUp) { // y unit vector in old coordinates
  /*
    first translate to the new cameraPosition
    then rotate to the new frame
    y' = cameraUp
    z' = cameraDirection
    x ' = y' x z'
   */
  Vec3f& y = cameraUp;
  Vec3f& z = cameraDirection;
  Vec3f x = cross(y,z);
  x.normalize();
  y.normalize();
  z.normalize();

  Matrix4f rot = Matrix4f::identity();
  for(std::size_t i = 0 ; i < 3 ; i++){
    rot(0, i) = x[i];
    rot(1, i) = y[i];
    rot(2, i) = z[i];      
  }
  
  Matrix4f tr = Matrix4f::identity();
  for(std::size_t i = 0 ; i < 3 ; i++) tr(i,3) = - cameraPosition[i];
  return rot * tr;
}

Matrix4f projection() {
  /*
    not using this anymore, is it correct (or does it give negative coordinates?)
   */
  Matrix4f m = Matrix4f::identity();
  //m(2,2) = 0;
  m(3,3) = 0;
  m(3,2) = -1;
  return m;
}

Matrix4f viewport(int width, int height) {
  /*
    not using anymore but correct for 90degree viewing angle
   */
  Matrix4f m = Matrix4f::identity();
  m(0,0) *= 0.5 * width;
  m(1,1) *= 0.5 * width;
  m(0,3) = 0.5 * width - 1;
  m(1,3) = 0.5 * height - 1;
  return m;
}

/*
  map a point to barycentric coordinates for a given triangle
  
  returns {u*denom, v*denom, denom} (to keep everything in terms of integers)
  the uvd coordinates of the point are (1-u-v, u, v)???
 */
int getBarycentricDenom(const PixelTriangle& t){
  return (t[1].x - t[0].x) * (t[2].y  - t[0].y )
    - (t[2].x-t[0].x) * (t[1].y  - t[0].y );
}
Vec3i getBarycentricUVD(const Point& p, const PixelTriangle& t){
  // convert to barycentric coordinates

  // (yC - yA) * (x - xA) - (y - yA) * (xC - xA)
  float u = (t[2].y-t[0].y) * (p.x - t[0].x)
    - (p.y - t[0].y) * (t[2].x - t[0].x);

  // (y - yA) * (xB - xA) - (yB - yA) * (x - xA)
  float v =  (p.y - t[0].y) * (t[1].x - t[0].x)
    - (t[1].y-t[0].y) * (p.x - t[0].x);
  // when is denom small?
  // (xB - XA) * (yC - yA) - (xC - xA) * (yB - yA)
  float denom = getBarycentricDenom(t);

  if (denom < 0){
    u *= -1;
    v *= -1;
    denom *= -1;
  }
  
  return {u,v,denom};
  // true u and v will be u/denom, v/denom
  // u = v = 0 corrsponds to vertex A
  // u = 1, v = 0 is vertex B
  // u = 0, v = 1 is vertex C
}

// returns (dU/dx, dV/dx, 0) to update uvd per pixel (NOT DIVIDED BY DENOM!)
Vec3i getBarycentricDeltaX(const PixelTriangle& t){
  return {t[2].y - t[0].y, t[0].y - t[1].y, 0};
}

// returns (dU/dy, dV/dy, 0) (NOT DIVIDED BY DENOM!)
Vec3i getBarycentricDeltaY(const PixelTriangle& t){
  return {t[0].x - t[2].x, t[1].x - t[0].x, 0};
}

/*
  A point is inside a triangle if:
  0 <= u
  0 <= v
  u + v <= 1

  params uvd = {u*denom, v*denom, denom}
  
 */
bool isPointInTriangle(Vec3i uvd){
  //if(uvd[2] < 0)
  // return (uvd[0] & uvd[1] & (uvd[2]-uvd[0]-uvd[1])) < 0;
    
  // the sign bit is 1 if negative
  return ((int) uvd.u | (int) uvd.v |(int) (uvd.d-uvd.u-uvd.v)) >= 0;
}

void Graphics::drawTriangle(const Triangle& t){
  // find bounding box
  Vec2i bboxmin = Point{width - 1, height - 1};
  Vec2i bboxmax = Point{0,0};
  for (int i = 0 ; i < 3 ; i ++){
    bboxmin.x = std::max(0, std::min(bboxmin[0], (int) t[i].x));
    bboxmin.y = std::max(0, std::min(bboxmin[1], (int) t[i].y));
    bboxmax.x = std::min(width-1, std::max(bboxmax[0], (int) t[i].x));
    bboxmax.y = std::min(height-1, std::max(bboxmax[1], (int) t[i].y));
  }

  auto uvd = getBarycentricUVD(Point{bboxmin.x, bboxmin.y}, t);
  auto dx = getBarycentricDeltaX(t);
  auto dy = getBarycentricDeltaY(t);

  // if the denom was negative, UVD got flipped, so need to reverse du/dx, etc
  // TODO refactor
  if(getBarycentricDenom(t) < 0){
    dx = - dx;
    dy = - dy;
  }  
  for (int x = bboxmin.x ; x <= bboxmax.x ; x++){
    auto uvd_start_column = uvd;
    for(int y = bboxmin.y ; y <= bboxmax.y ; y++){
      if (isPointInTriangle(uvd)){
	float z = t[0].z * (uvd.d-uvd.u-uvd.v) + t[1].z * uvd.u + t[2].z * uvd.v;
	z /= uvd[2]; // denom
	
	if (z < zbuffer[x*height + y]){
	  SDL_RenderDrawPoint(renderer, x, height - 1 - y);
	  zbuffer[x*height + y] = z;
	}
      }
      uvd += dy;
    }
    uvd = uvd_start_column + dx;
  }
}

void drawTriangleFilled(Triangle& t, SDL_Renderer* renderer){
  // from bottom point to top point
  std::sort(&t[0], &t[0] + 3);

  // draw horizontal lines from left to right

  // so ideally instead of the slope/float code below I would do something
  // similar to the drawLine rasterization
 
  float p02slope = ((float) (t[2].y - t[0].y)) / (t[2].x - t[0].x);
  float p01slope = ((float) (t[1].y - t[0].y)) / (t[1].x - t[0].x);
  float p12slope = ((float) (t[2].y - t[1].y)) / (t[2].x - t[1].x);
  int x_incr = t[1].x < t[0].x + (t[1].y-t[0].y)/p02slope ? -1 : 1;
  
  for(int y = t[0].y ; y <= t[2].y ; y++){
    
    int x_start = t[0].x + (y - t[0].y)/p02slope;
    int x_end;
    if (y <= t[1].y)
      x_end = t[0].x + (y - t[0].y)/p01slope;
    else
      x_end = t[1].x + (y - t[1].y)/p12slope;
    
    while(x_start != x_end){
      // height - 1 - y
      SDL_RenderDrawPoint(renderer, x_start, y);
      x_start += x_incr;
    }
  }
}


void drawLine(Point p0, Point p1, SDL_Renderer* renderer){
  bool steep = false;
  
  // find steeper side and transpose drawing
  // increment along less steep side to avoid gaps
  if (std::abs(p1[1]-p0[1]) > std::abs(p1[0]-p0[0])){
    steep = true;
    std::swap(p0[0], p0[1]);
    std::swap(p1[0], p1[1]);
  }

  // draw from left to right
  if (p0[0] > p1[0])
    std::swap(p0, p1);

  int dx = p1[0] - p0[0];
  int dy = p1[1] - p0[1];

  int derror2 = std::abs(dy) * 2;
  int error2 = 0;
  int yincr = p1[1] > p0[1] ? 1 : -1; //pos/neg slope

  // height - 1 - y
  for(Point p(p0) ; p[0] <= p1[0] ; p[0]++){
    if(steep) // modern compilers can optimize this out of the loop
      SDL_RenderDrawPoint(renderer, p[1], p[0]);
    else
      SDL_RenderDrawPoint(renderer, p[0], p[1]);
    
    error2 += derror2;
    if(error2 > dx){
      p[1] += yincr;
      error2 -= dx*2;
    }
  }
}

void drawTriangleBoundary(const Triangle& t, SDL_Renderer* renderer){
  drawLine({t[0].x, t[0].y}, {t[1].x, t[1].y}, renderer);
  drawLine({t[1].x, t[1].y}, {t[2].x, t[2].y}, renderer);
  drawLine({t[2].x, t[2].y}, {t[0].x, t[0].y}, renderer);
}
