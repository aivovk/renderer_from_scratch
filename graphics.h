#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <SDL2/SDL.h>

#include <vector>
#include <functional>
#include <cassert>

#include "geometry.h"
#include "scene.h"

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

/* Converts a point's x,y pixel coordinates to a triangle's barycentric
 * coordinate system (based on x,y only)
 *
 * Params:
 * p - the point
 * t - a Triangle ( 3 points )
 *
 * Returns:
 * (U,V,D) where the real coordinates of the point are:
 * (1-U/D-V/D) * A + U/D * B + V/D * C
 * where A,B,C are the x,y coordinates of the triangle's points
 */
Vec3i getBarycentrivUVD(const Point& p, const PixelTriangle& t);
// just the denominator (to test sign)
int getBarycentricDenom(const PixelTriangle& t);
// per pixel changes in U and V
Vec3i getBarycentrivDeltaX(const PixelTriangle& t);
Vec3i getBarycentrivDeltaY(const PixelTriangle& t);

/*
 * See getBarycentricUVD for definition of uvd
 */
inline bool isPointInTriangle(Vec3i uvd);

Matrix4f cameraLookAt(Vec3f cameraPosition,
		      Vec3f cameraDirection,
		      Vec3f cameraUp);


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
  
  void clearZBuffer();
  void convertToScreen(Vec4f& v) const;
  void convertToScreen(Triangle& t) const;
  
  Graphics()=delete;
  Graphics(SDL_Renderer* renderer, Matrix4f& camera_transform):
    renderer(renderer),
    camera_transform(camera_transform) {
    SDL_GetRendererOutputSize(renderer, &width, &height);
    clearZBuffer();
  }

  void addModel(const Model& m);
  void addScene(const Scene& s);

  // render models
  void render();

  void drawLine(Point p0, Point p1);
  void drawTriangle(const Triangle&, const Color vertexColors[3]);
  void drawTriangleBoundary(const Triangle&);  
};



#endif
