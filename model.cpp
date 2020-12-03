#include "model.h"

void Model::loadXYZ(std::string filename){
  std::ifstream f (filename);
  std::string line;
  if (f.is_open())
    {
      size_t num_lines = 0;
      while ( getline (f,line) )
	{
	  if(line[0] == '#')
	    continue;
	    
	  std::stringstream ss(line);
	  float x, y, z;
	  ss >> x >> y >> z;
	  vertices.push_back({x, y, z, 1.0f});
	  num_lines++;
	  if(num_lines % 3 == 0){
	    triangles.push_back(TriangleIndex{num_lines-3, num_lines-2, num_lines-1});
	    triangleNormals.push_back(TriangleIndex{num_lines-3, num_lines-2, num_lines-1});
	    auto t = getTriangle(size() - 1);
	    auto n = cross(t[1] - t[0],
			   t[2] - t[1]).normalize();
	    for(int i = 0 ; i < 3; i++)
	      vertexNormals.push_back(n);
	  }
	}
      f.close();
    }

}

void Model::loadOBJ(std::string filename){
  std::ifstream f (filename);
  std::string line;
  char junk;
  
  if (f.is_open()){
    while(getline(f, line)){
      std::stringstream ss(line);      
      // vertex
      if(line.compare(0, 2, "v ") == 0){
	ss >> junk;
	float x, y, z;
	ss >> x >> y >> z;
	vertices.push_back({x,y,z,1.0f});
      }
      // vertex normal
      else if (line.compare(0,2, "vn") == 0){
	ss >> junk >> junk;
	float x, y, z;
	ss >> x >> y >> z;
	vertexNormals.push_back({x,y,z});
      }      
      // vertex texture coordinate
      else if (line.compare(0,2, "vt") == 0){
	ss >> junk >> junk;
	
      }      
      // faces (triangles)
      else if (line.compare(0,2, "f ") == 0){
	ss >> junk;
	TriangleIndex t;
	TriangleIndex tn;
	std::size_t iVertex, iTexture, iNormal;
	for(int i = 0 ; i < 3 ; i++){
	  ss >> iVertex >> junk >> iTexture >> junk >> iNormal;
	  t[i] = iVertex - 1;
	  tn[i] = iNormal - 1;
	}
	triangles.push_back(t);
	triangleNormals.push_back(tn);
      }
    }
    f.close();
  }
}
