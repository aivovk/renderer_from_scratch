#ifndef MODEL_H
#define MODEL_H
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include "graphics.h"

class Model{
private:
  std::vector<Triangle> triangles;

  // the normal comes out of the plane where vertices are in a clockwise order
  std::vector<Vec3f> normals;
  // TODO
  // bool twoSided = false;

  // TODO: eventually will want to store the rotation/scale axes
  
  std::vector<Color> triangle_colors;
public:
  Model() = delete;
  Model(std::string filename){
    std::ifstream f (filename);
    std::string line;
    if (f.is_open())
      {
	size_t num_lines = 0;
	Triangle t;
	while ( getline (f,line) )
	  {
	    if(line[0] == '#')
	      continue;
	    
	    std::stringstream ss(line);
	    float x, y, z;
	    ss >> x >> y >> z;
	    t[num_lines%3] = {x, y, z, 1.0f};
	    num_lines++;
	    if(num_lines % 3 == 0){
	      triangles.push_back(t);
	      normals.push_back(cross(Vec3f{t[1] - t[0]}, Vec3f{t[2] - t[1]}).normalize());
	    }
	  }
	f.close();
      }
  }
  const std::vector<Triangle>& getTriangles() const { return triangles; }
  const std::vector<Color>& getColors() const { return triangle_colors; }
  const std::vector<Vec3f>& getNormals() const { return normals; }

  void scale(const Vec3f& s){
    for (auto& t : triangles)
      for(int i = 0 ; i < 3 ; i++)
	t[i] = mult(t[i], s);
  }
  void rotate(const Vec3f& r){
    
  }
  void translate(const Vec3f& t){
    for(auto& tri : triangles)
      for(int i = 0 ; i < 3 ; i++)
	tri[i] += t;
  }
};

#endif
