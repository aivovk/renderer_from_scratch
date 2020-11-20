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
  std::vector<Vec3f> normals;
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
	    std::stringstream ss(line);
	    float x, y, z;
	    ss >> x >> y >> z;
	    t[num_lines%3] = {x, y, z, 1.0f};
	    num_lines++;
	    if(num_lines % 3 == 0)
	      triangles.push_back(t);
	  }
	f.close();
      }
  }
  const std::vector<Triangle>& getTriangles() const { return triangles; }
  const std::vector<Color>& getColors() const { return triangle_colors; }
};

#endif
