/*
 * A model is a 3D shape defined by its vertices, triangles, vertex normals, and
 * possibly colours, textures, material properties...
 * 
 * Currently loaded from a .xyz file or .obj file
 * 
 * Can be rotated, translated, and scaled as a unit.
 *
 */
#ifndef MODEL_H
#define MODEL_H
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <cmath>
#include <initializer_list>
#include "geometry.h"
#include "quaternion.h"

// indices to 3 vertices in the vertices array
struct TriangleIndex{
  std::size_t data[3];
  TriangleIndex(){}
  TriangleIndex(std::initializer_list<std::size_t> data_list) {
    int i = 0;
    for (auto iter = data_list.begin() ; i < 3 ; i++, iter++)
      data[i] = *iter;
  }

  inline std::size_t& operator[](std::size_t i){
    return data[i];
  }
  inline const std::size_t& operator[](std::size_t i) const{
    return data[i];
  }
};

class Model{
private:
  std::vector<Vec4f> vertices;
  std::vector<Vec3f> vertexNormals;

  std::vector<TriangleIndex> triangles;
  std::vector<TriangleIndex> triangleNormals;
  
  // std::vector<Color> triangle_colors;
  // textures
  
  void loadXYZ(std::string filename);
  void loadOBJ(std::string filename);
  
  // bool twoSided = false;

  // TODO: store translate amount/origin to rescale/rotate again
public:
  
  Model() = delete;
  Model(std::string filename){
    std::cout<<"Loading model: "<< filename;
    if(filename.compare(filename.size() - 3, 3, "xyz") == 0){
      std::cout<<" (XYZ file)\n";
      loadXYZ(filename);
    }
    else if(filename.compare(filename.size() - 3, 3, "obj") == 0){
      std::cout<<" (OBJ file)\n";
      loadOBJ(filename); 
    }else
      std::cout<<", NO MODEL LOADED (unknown file extension)\n";
  }
  std::size_t size() const {return triangles.size();}
  Triangle getTriangle(std::size_t i) const {
    Triangle t;
    t[0] = vertices[triangles[i][0]];
    t[1] = vertices[triangles[i][1]];
    t[2] = vertices[triangles[i][2]];
    return t;
  }
  //const std::vector<Color>& getColors() const { return triangle_colors; }
  Triangle getTriangleNormals(std::size_t i) const {
    Triangle tn;
    tn[0] = vertexNormals[triangleNormals[i][0]];
    tn[1] = vertexNormals[triangleNormals[i][1]];
    tn[2] = vertexNormals[triangleNormals[i][2]];
    return tn;
  }

  void scale(const Vec3f& s){
    for (auto& v : vertices)
      v = mult(v, s);
  }
  void rotate(const Vec3f& r){
    for(auto& v : vertices){
	v.rotate(Vec3f{1,0,0}, r.x);
	v.rotate(Vec3f{0,1,0}, r.y);
	v.rotate(Vec3f{0,0,1}, r.z);
    }
    for(auto& vn : vertexNormals){
	vn.rotate(Vec3f{1,0,0}, r.x);
	vn.rotate(Vec3f{0,1,0}, r.y);
	vn.rotate(Vec3f{0,0,1}, r.z);
    }

  }
  void translate(const Vec3f& t){
    for(auto& v : vertices)
      v += t;
  }
};

#endif
