#ifndef SCENE_H
#define SCENE_H
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include "model.h"

class Scene{
private:
  std::vector<Model> models;
  enum LineDescription{FILENAME = 0, SCALE = 1, ROTATE = 2, TRANSLATE = 3, NUM_LINES = 4};
public:
  Scene() = delete;
  
  Scene(std::string filename){
    std::ifstream f (filename);
    /*
      file format is the following:
      # comment lines must have a # as the first character
      # the rest of the lines must be in groups of 4 and the following order
      # filename
      model.xyz
      # scale 
      sx sy sz
      # rotation about each axis in radians
      rx ry rz
      # translation
      tx ty tz
     */

    
    std::string line;
    if (f.is_open())
      {
	size_t num_lines = 0;
	while ( getline (f,line) )
	  {
	    if(line[0] == '#')
	      continue;

	    if (num_lines % NUM_LINES == FILENAME){
	      models.push_back(Model(line));
	    }else{ 
	      std::stringstream ss(line);
	      float x, y, z;
	      ss >> x >> y >> z;
	      if( num_lines % NUM_LINES == SCALE)
		models.back().scale({x, y, z});
	      else if ( num_lines % NUM_LINES == ROTATE)
		models.back().rotate({x, y, z});
	      else if ( num_lines % NUM_LINES == TRANSLATE)
		models.back().translate({x, y, z});
	    }
	    num_lines++;
	  }
	f.close();
      }
  }
  const std::vector<Model>& getModels() const { return models; }

};

#endif
