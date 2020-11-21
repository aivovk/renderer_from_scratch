#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <SDL2/SDL.h>

#include <vector>
#include <functional>

#include "geometry.h"
#include "scene.h"

/* Converts a point's x,y coordinates to a triangles barycentric coordinate
 * system (based on x,y only)
 *
 * Params:
 * p - the point
 * ps - a Triangle ( 3 points )
 *
 * Returns:
 * (U,V,D) where the real coordinates of the point are: (1-U/D-V/D) * A + U/D * B + V/D * C
 * where A,B,C are the x,y coordinates of the triangle's points
 */
Vec3i getBarycentricUVD(const Point& p, const Vec4f (&ps)[3]);

/*
 * See getBarycentricUVD for definition of uvd
 */
inline bool isPointInTriangle(Vec3i uvd);

//void drawTriangle(const Triangle&, SDL_Renderer* renderer);
void drawTriangleBoundary(const Triangle& ps, SDL_Renderer* renderer);
void drawTriangleFilled(Triangle& ps, SDL_Renderer* renderer);

// should be named isAtLeastOnePointInFrontOfTheCamera
//bool isTriangleVisible(const Triangle& t);

/*
 * Find the matrix that transforms world coordinates into a frame where the
 * camera is at (0,0,0) looking along +z with +y being up
 */
Matrix4f cameraLookAt(Vec3f cameraPosition,
		      Vec3f cameraDirection,
		      Vec3f cameraUp);

/*
 * Vanishing point projection along the z axis
 */
Matrix4f projection();

/*
 * Convert to screen coordinates
 */
Matrix4f viewport(int width, int height);

class Graphics{
private:
  enum class Visibility{INVISIBLE, VISIBLE, PARTIAL};

 
  SDL_Renderer* renderer;
  const Matrix4f& camera_transform;
  
  int width;
  int height;
  std::vector<float> zbuffer;

  std::vector<std::reference_wrapper<const Model>> models;
  
public:
  static constexpr float NEAR_CLIP_DIST = 0.1f;
  static constexpr float HORIZONTAL_VIEWING_ANGLE_RAD = 0.5 * M_PI;
  static constexpr float VIEWING_FACTOR = 1.0f / tan(0.5f * HORIZONTAL_VIEWING_ANGLE_RAD);
  
  void clearZBuffer(){
    zbuffer.assign(width*height, std::numeric_limits<float>::max());
  }

  void convertToScreen(Vec4f& v) const {
    // when x/z = tan(view_angle / 2)
    v[0] = width/2 - 1 + width/2 * v[0]/v[2] * VIEWING_FACTOR;
    v[1] = height/2 - 1 + width/2 * v[1]/v[2] * VIEWING_FACTOR;
  }
  
  void convertToScreen(Triangle& t) const {
    for(size_t i = 0 ; i < 3 ; i++)
      convertToScreen(t[i]);
  }
  
  Graphics()=delete;
  Graphics(SDL_Renderer* renderer, Matrix4f& camera_transform):
    renderer(renderer),
    camera_transform(camera_transform) {
    SDL_GetRendererOutputSize(renderer, &width, &height);
    clearZBuffer();
    camera_transform = cameraLookAt({0,0,0}, {0,0,1}, {0,1,0});
    //viewport_and_projection = viewport(width, height) * projection();
  }
  void drawTriangle(const Triangle& ps);

  void addModel(const Model& m){
    models.push_back(std::cref(m));
  }

  void addScene(const Scene& s){
    for(const auto& m : s.getModels())
      addModel(m);
  }

  void render(){
    /*
      do the following steps:

      - transform from World to Camera coordinates
      - where do i want to store these? graphics should keep a list of the models it wants to draw an array for the triangles and their visibility

      - clip the triangles : mark the triangles that intersect the z axis as invisible (or the view volume)
      - recalculate new vertices and triangles matching those

      - do i even need the arrays?

      - just do it all in one step:
      transform the triangle
      check if it's visible 
      do necessary calculations to recompute vertices
      draw it
      

     */
    for (size_t i_model = 0 ; i_model < models.size() ; i_model++){
      const auto& triangles = models[i_model].get().getTriangles();
      const auto& normals = models[i_model].get().getNormals();
      for(size_t i_triangle = 0 ; i_triangle < triangles.size() ; i_triangle++){
	const auto& t = triangles[i_triangle];

	float lightIntensity = 255.0f;
	SDL_SetRenderDrawColor(renderer,
			       lightIntensity,
			       lightIntensity,
			       lightIntensity,
			       255);
	
	Triangle screen_triangle;
	int num_visible_vertices = 0;
	for(size_t i_vertex = 0 ; i_vertex < 3 ; i_vertex++){
	  screen_triangle[i_vertex] = camera_transform * t[i_vertex];
	  if (screen_triangle[i_vertex][2] > NEAR_CLIP_DIST){
	    num_visible_vertices++;
	  }
	}	
	if (num_visible_vertices > 0){
	  // recalculate the normal (seems faster than transforming?)
	  auto n = cross(Vec3f{screen_triangle[1] - screen_triangle[0]},
			 Vec3f{screen_triangle[2] - screen_triangle[1]}).normalize();

	  
	  // camera and normal point in the same direction if face is invisible
	  /* 
	     this is, I think too much, all visible faces that should be drawn
	     are drawn, as well as some that shouldn't be (i think based on how
	     far it is from the viewing box edges)

	     do i have to actually transform the normals to the projected
	     coordinates?
	   */
	  if(n * Vec3f{0, 0, 1} > acos(0.5 * HORIZONTAL_VIEWING_ANGLE_RAD))
	    continue;		 
	  
	  if(num_visible_vertices == 1){
	    /*
	      0 1 2
	      - - +
	      + - -
	      - + -
	      3 cases , we find the first positive index, then calculate the intersections between the line segments between it and the other two points and and the near_clip_dist plane
	     */
	    //std::cout<<"1 visible\n";
	    
	    int i_pos = 0;
	    while (screen_triangle[i_pos][2] <= NEAR_CLIP_DIST) i_pos++;
	    // first vertex
	    // x
	    screen_triangle[(i_pos+1)%3][0] = screen_triangle[i_pos][0]
	      + (screen_triangle[(i_pos+1)%3][0]-screen_triangle[i_pos][0])
	      * (NEAR_CLIP_DIST-screen_triangle[i_pos][2])
	      / (screen_triangle[(i_pos+1)%3][2]-screen_triangle[i_pos][2]);
	    // y
	    screen_triangle[(i_pos+1)%3][1] = screen_triangle[i_pos][1]
	      + (screen_triangle[(i_pos+1)%3][1]-screen_triangle[i_pos][1])
	      * (NEAR_CLIP_DIST-screen_triangle[i_pos][2])
	      / (screen_triangle[(i_pos+1)%3][2]-screen_triangle[i_pos][2]);
	    // z
	    screen_triangle[(i_pos+1)%3][2] = NEAR_CLIP_DIST;
	    
	    // second vertex
	    // x
	    screen_triangle[(i_pos+2)%3][0] = screen_triangle[i_pos][0]
	      + (screen_triangle[(i_pos+2)%3][0]-screen_triangle[i_pos][0])
	      * (NEAR_CLIP_DIST-screen_triangle[i_pos][2])
	      / (screen_triangle[(i_pos+2)%3][2]-screen_triangle[i_pos][2]);
	    // y
	    screen_triangle[(i_pos+2)%3][1] = screen_triangle[i_pos][1]
	      + (screen_triangle[(i_pos+2)%3][1]-screen_triangle[i_pos][1])
	      * (NEAR_CLIP_DIST-screen_triangle[i_pos][2])
	      / (screen_triangle[(i_pos+2)%3][2]-screen_triangle[i_pos][2]);
	    screen_triangle[(i_pos+2)%3][2] = NEAR_CLIP_DIST;

	    SDL_SetRenderDrawColor(renderer, 155, 155, 155, 255);
	    convertToScreen(screen_triangle);
	    drawTriangle(screen_triangle);
	   
	  } else if (num_visible_vertices == 2) {
	    /*
	      0 1 2
	      + + -
	      - + +
	      + - +
	     */
	    //std::cout<<"2 visible\n";
	    // find invisible vertex
	    int i_pos = 0;
	    while (screen_triangle[i_pos][2] > NEAR_CLIP_DIST) i_pos++;
	    
	    Triangle extra_triangle;
	    extra_triangle[0] = screen_triangle[(i_pos+1)%3];

	    // first new vertex
	    //x
	    extra_triangle[1][0] = extra_triangle[0][0]
	      + (screen_triangle[i_pos][0]-extra_triangle[0][0])
	      * (NEAR_CLIP_DIST-extra_triangle[0][2])
	      / (screen_triangle[i_pos][2]-extra_triangle[0][2]);
	    // y
	    extra_triangle[1][1] = extra_triangle[0][1]
	      + (screen_triangle[i_pos][1]-extra_triangle[0][1])
	      * (NEAR_CLIP_DIST-extra_triangle[0][2])
	      / (screen_triangle[i_pos][2]-extra_triangle[0][2]);
	    // z
	    extra_triangle[1][2] = NEAR_CLIP_DIST;

	    //second new vertex
	    //x
	    extra_triangle[2][0] = screen_triangle[(i_pos+2)%3][0]
	      + (screen_triangle[i_pos][0]-screen_triangle[(i_pos+2)%3][0])
	      * (NEAR_CLIP_DIST-screen_triangle[(i_pos+2)%3][2])
	      / (screen_triangle[i_pos][2]-screen_triangle[(i_pos+2)%3][2]);
	    // y
	    extra_triangle[2][1] = screen_triangle[(i_pos+2)%3][1]
	      + (screen_triangle[i_pos][1]-screen_triangle[(i_pos+2)%3][1])
	      * (NEAR_CLIP_DIST-screen_triangle[(i_pos+2)%3][2])
	      / (screen_triangle[i_pos][2]-screen_triangle[(i_pos+2)%3][2]);
	    // z
	    extra_triangle[2][2] = NEAR_CLIP_DIST;

	    // correct other triangle
	    screen_triangle[i_pos] = extra_triangle[2];

	    SDL_SetRenderDrawColor(renderer, 0, 255, 255, 255);
	    convertToScreen(screen_triangle);
	    drawTriangle(screen_triangle);
	    SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255);
	    convertToScreen(extra_triangle);
	    drawTriangle(extra_triangle);
	    	    
	  } else {
	    //std::cout<<"3 visible\n";
	    convertToScreen(screen_triangle);
	    drawTriangle(screen_triangle);
	  }
	}
      }
    }
    

  }
};



#endif
