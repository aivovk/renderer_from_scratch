/*
 * A model is a collection of triangles, loaded from a .xyz file
 * 
 * Can be rotated, translated, and scaled.
 *
 * TODO:
 * - colours (or textures)
 * - normals
 * - material properties
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
  /*
    PLAN: 
    store the vertices once
    store the vertex normals once
    
    triangles are indices into vertices (and vertex normals)
    this is slowish, since we need to look up the vertices and copy the triangle
    space for n triangles:
    // models are solid and vertices are not repeated
    n * sizeof(vertex) // for vertices
    n * sizeof(normal) // for normals
    n * 3 * sizeof(index) // ideally we would store the vertices and corresponding
                          // normals, but might need two arrays

    Total size, assuming sizeof(vertex)== 4, sizeof(normal)==3, sizeof(index) = 1
    10n
			  
    alternative:
    store by triangle, if we order the normals the same way triangle, then we will need
    n * 3 * sizeof(vertex)
    n * 3 * sizeof(normal)
    Total size:
    21n

    Another consideration:
    each step we will be doing
    cameraTransorm * vertex
    and 
    cameraTransorm * normal

    where do we want to do this?
    ideally the model class should just store the indices
    the graphics class stores a space efficient representation of the vertices and normals

    so will need to change the constructors:
    Scene
    g.loadScene() // this function returns it's vertices to graphics
    loadModel() // returns it's triangles, which get fed to graphics

    scene just stores it's models
    and models just store the indices of their vertices and normals (ideally these are ordered the same)
    graphics stores:
    vertices
    vertices in camera coordinates
    normals
    normals in camera coordinates

    graphics pipeline:
    cameraMatrix * vertices
    cameraMatrix^{-1}^{T} * normals (don't need to do this? lights are in world coordinates? also, don't need to do this unless doing uneven scaling)
    
    for each triangle:
    cull it by location (or split it into several)
    cull it by normal direction (calculated from vertices)

    vertex shader
    transorm triangle to screen coordinates (or do we want to mark the vertices and do it there?)
    rasterizer
    fragment shader + texture lookup

    a problem though : say i want to deform a model, i can only look up the
    triangles, meaning i will transform the vertices approximately 3 times each
    unless i mark them somehow

    what is the alternative? just store the triangles in the model, and then do
    the transformations per triangle this uses more space and requires 3
    operations per vertex (doesn't require finding the model), will figure out
    how to optimize this later
    
   */
  std::vector<Vec4f> vertices;
  std::vector<Vec3f> vertexNormals;

  std::vector<TriangleIndex> triangles;
  std::vector<TriangleIndex> triangleNormals; //vertex normals
  
  // std::vector<Color> triangle_colors;
  // textures
  
  void loadXYZ(std::string filename);
  void loadOBJ(std::string filename);

  // TODO
  // bool twoSided = false;
  // TODO: eventually will want to store the rotation/scale axes or the translate amount
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
