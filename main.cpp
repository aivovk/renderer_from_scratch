#include <iostream>

#include "graphics.h"

void printRenderDriverInfo(){
    SDL_RendererInfo info;
    SDL_GetRenderDriverInfo(0, &info);
    std::cout<<info.name<<std::endl;
    std::cout<<(info.flags & SDL_RENDERER_SOFTWARE)<<std::endl;
    std::cout<<(info.flags & SDL_RENDERER_ACCELERATED)<<std::endl;
    std::cout<<(info.flags & SDL_RENDERER_PRESENTVSYNC)<<std::endl;
    std::cout<<(info.flags & SDL_RENDERER_TARGETTEXTURE)<<std::endl;
    std::cout<<info.max_texture_width<<std::endl;
    std::cout<<info.max_texture_height<<std::endl;
}

int main(int argc, char** argv){
  SDL_Event event;

  // graphics will get the actual width and height from the renderer
  const int width_request = 800;
  const int height_request = 600;
  
  const int maxFPS = 24;
  const float maxFrameTicks = 1000.0/maxFPS; // milliseconds
  int currentTick, delayTicks;
  int lastTick = SDL_GetTicks();
  int frame = 0;

  SDL_Init(SDL_INIT_VIDEO);
  SDL_Window* window = SDL_CreateWindow("Window",
					0, // window location
					0,
					width_request,
					height_request,
 					SDL_WINDOW_SHOWN);// | SDL_WINDOW_FULLSCREEN_DESKTOP);
  
  SDL_Renderer* renderer = SDL_CreateRenderer(window,
					      -1, // first renderer matching flags
					      SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
  SDL_SetRelativeMouseMode(SDL_TRUE);
  //printRenderDriverInfo();

  
  bool quit = false;

  std::string filename;
  if(argc > 1)
    filename = argv[1];
  else
    filename = "s1.scene";
  Scene s(filename);
  
  Vec3f camera = {0,1,-5};
  Vec3f up = {0,1,0};
  Vec3f cameraDir = {0,0,1};
  Matrix4f camera_transform;

  Graphics g(renderer, camera_transform);
  g.addScene(s);
 
  float cameraLookSpeed = 0.005f; // radians / (pixel * frame duration)???
  float cameraMoveSpeed = 1.0f; // program distance units / frame duration?
  while(!quit){
    while(SDL_PollEvent(&event)){
      if(event.type == SDL_QUIT)
	quit = true;
      if(event.type==SDL_MOUSEMOTION){
	float xrel = event.motion.xrel * cameraLookSpeed;
	float yrel = event.motion.yrel * cameraLookSpeed;

	// left/right camera movement
	cameraDir.rotate(Vec3f{0,1,0}, xrel);
	up.rotate(Vec3f{0,1,0}, xrel);

	// up/down camera movement
	Vec3f cameraPerp = cross(up, cameraDir);
	cameraDir.rotate(cameraPerp, yrel);
	up.rotate(cameraPerp, yrel);

	cameraDir.normalize();
	up.normalize();
      }
      if(event.type == SDL_KEYDOWN){
	switch (event.key.keysym.sym)
	  {
	  case SDLK_w:
	    camera += cameraDir * cameraMoveSpeed;
	    break;
	  case SDLK_a:
	    camera -= cross(up,cameraDir) * cameraMoveSpeed;
	    break;
	  case SDLK_s:
	    camera -= cameraDir * cameraMoveSpeed;
	    break;
	  case SDLK_d:
	    camera += cross(up,cameraDir) * cameraMoveSpeed;
	    break;
	  }
      }
    }
    // clear
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    g.clearZBuffer();

    //draw
    camera_transform = cameraLookAt(camera, cameraDir, up);
    g.render();
    SDL_RenderPresent(renderer);

    // timing
    frame++;
    currentTick = SDL_GetTicks();
    delayTicks = std::max(0.0f, maxFrameTicks - (currentTick - lastTick));
    SDL_Delay(delayTicks);
    lastTick = SDL_GetTicks();
  }

  
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
