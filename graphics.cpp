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
  /*
    [x y z]

    // todo can calculate first column and then use it in the determinant calculation
   */
  //float det = x[0] * (y[1]*z[2] - y[2]*z[1]) + y[0]*(z[1]*x[2] - z[2]*x[1]) +z[0]*(x[1]*y[2]-x[2]*y[1]);
  for(std::size_t i = 0 ; i < 3 ; i++){
    /*rot(0, row) = (y[(row+1)%3]*z[(row+2)%3] - y[(row+2)%3]*z[(row+1)%3])/det;
    rot(1, row) = (z[(row+1)%3]*x[(row+2)%3] - z[(row+2)%3]*x[(row+1)%3])/det;
    rot(2, row) = (x[(row+1)%3]*y[(row+2)%3] - x[(row+2)%3]*y[(row+1)%3])/det;*/
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
bool isTriangleVisible(const Triangle& t){
  return t[0][2] > Graphics::NEAR_CLIP_DIST
    || t[1][2] > Graphics::NEAR_CLIP_DIST
    || t[2][2] > Graphics::NEAR_CLIP_DIST;
    }*/

/*
  map a point to barycentric coordinates for a given triangle
  
  returns {u*denom, v*denom, denom} (to keep everything in terms of integers)
  the uvd coordinates of the point are (1-u-v, u, v)???
 */
Vec3i getBarycentricUVD(const Point& p, const Triangle& ps){
  // convert to barycentric coordinates

  // (yC - yA) * (x - xA) - (y - yA) * (xC - xA)
  int u = (ps[2][1]-ps[0][1]) * (p[0] - ps[0][0])
    - (p[1] - ps[0][1]) * (ps[2][0] - ps[0][0]);

  // (y - yA) * (xB - xA) - (yB - yA) * (x - xA)
  int v =  (p[1] - ps[0][1]) * (ps[1][0] - ps[0][0])
    - (ps[1][1]-ps[0][1]) * (p[0] - ps[0][0]);
  // when is denom small?
  // (xB - XA) * (yC - yA) - (xC - xA) * (yB - yA)
  int denom = (ps[1][0] - ps[0][0]) * (ps[2][1] - ps[0][1])
    - (ps[2][0]-ps[0][0]) * (ps[1][1] - ps[0][1]);

  /*if (denom < 0){
    u *= -1;
    v *= -1;
    denom *= -1;
    }*/
  
  return {u,v,denom}; // true u and v will be u/denom, v/denom
  // coordinates are (1-u-v, u, v)/denom
}

Vec3i getBarycentricDeltaX(const Triangle& ps){
  return {ps[2][1], -ps[1][1], 0};
}
Vec3i getBarycentricDeltaY(const Triangle& ps){
  return {-ps[2][0], ps[1][0], 0};
}

/*
  A point is inside a triangle if:
  0 <= u
  0 <= v
  u + v <= 1

  params uvd = {u*denom, v*denom, denom}
  
 */
bool isPointInTriangle(Vec3i uvd){
  if(uvd[2] < 0)
   return (uvd[0] & uvd[1] & (uvd[2]-uvd[0]-uvd[1])) <= 0;
    
  // the sign bit is 1 if negative
  return (uvd[0] | uvd[1] | (uvd[2]-uvd[0]-uvd[1])) >= 0;
}

void Graphics::drawTriangle(const Triangle& ps){
  // find bounding box
  Point bboxmin = Point{width - 1, height - 1};
  Point bboxmax = Point{0,0};
  for (int i = 0 ; i < 3 ; i ++){
    bboxmin[0] = std::max(0, std::min(bboxmin[0], (int) ps[i][0]));
    bboxmin[1] = std::max(0, std::min(bboxmin[1], (int) ps[i][1]));
    bboxmax[0] = std::min(width-1, std::max(bboxmax[0], (int) ps[i][0]));
    bboxmax[1] = std::min(height-1, std::max(bboxmax[1], (int) ps[i][1]));
  }

  auto uvd = getBarycentricUVD(Point{bboxmin[0], bboxmin[1]}, ps);
  auto dx = getBarycentricDeltaX(ps);
  auto dy = getBarycentricDeltaY(ps);
  
  for (int x = bboxmin[0] ; x <= bboxmax[0] ; x++){
    auto uvd_start_column = uvd;
    for(int y = bboxmin[1] ; y <= bboxmax[1] ; y++){
      //auto uvd = getBarycentricUVD(Point{x,y}, ps);
      if (isPointInTriangle(uvd)){
	/*
	  z is now guaranteed to be in front of viewing plane so don't need to do any extra checks
	*/
	
	float z = ps[0][2] * (uvd[2]-uvd[0]-uvd[1]) + ps[1][2] * uvd[0] + ps[2][2] * uvd[1];
	z /= uvd[2]; // denom
	
	if (z < zbuffer[x*height + y]){
	  SDL_RenderDrawPoint(renderer, x, height - 1 - y);
	  zbuffer[x*height + y] = z;
	}
      }
      uvd += dy;
    }
    uvd = uvd_start_column;
    uvd += dx;
  }
}

void drawTriangleFilled(Triangle& ps, SDL_Renderer* renderer){
  // from bottom point to top point
  std::sort(&ps[0], &ps[0] + 3);

  // draw horizontal lines from left to right

  // so ideally instead of the slope/float code below I would do something
  // similar to the drawLine rasterization
 
  float p02slope = ((float) (ps[2][1] - ps[0][1])) / (ps[2][0] - ps[0][0]);
  float p01slope = ((float) (ps[1][1] - ps[0][1])) / (ps[1][0] - ps[0][0]);
  float p12slope = ((float) (ps[2][1] - ps[1][1])) / (ps[2][0] - ps[1][0]);
  int x_incr = ps[1][0] < ps[0][0] + (ps[1][1]-ps[0][1])/p02slope ? -1 : 1;
  
  for(int y = ps[0][1] ; y <= ps[2][1] ; y++){
    
    int x_start = ps[0][0] + (y - ps[0][1])/p02slope;
    int x_end;
    if (y <= ps[1][1])
      x_end = ps[0][0] + (y - ps[0][1])/p01slope;
    else
      x_end = ps[1][0] + (y - ps[1][1])/p12slope;
    
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

void drawTriangleBoundary(const Triangle& ps, SDL_Renderer* renderer){
  drawLine({ps[0][0], ps[0][1]}, {ps[1][0], ps[1][1]}, renderer);
  drawLine({ps[1][0], ps[1][1]}, {ps[2][0], ps[2][1]}, renderer);
  drawLine({ps[2][0], ps[2][1]}, {ps[0][0], ps[0][1]}, renderer);
}
