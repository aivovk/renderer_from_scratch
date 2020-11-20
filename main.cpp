#include <iostream>

#include "graphics.h"
#include "model.h"

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

int main(){
  SDL_Event event;
  const int width = 800;
  const int height = 600;
  const int FPS = 24;
  float frameTicks = 1000.0/FPS; // milliseconds
  int lastTick, currentTick, delayTicks;
  int frame = 0;

  SDL_Init(SDL_INIT_VIDEO);
  SDL_Window* window = SDL_CreateWindow("Window",
					280,
					75,
					width,
					height,
 					SDL_WINDOW_SHOWN);
  SDL_Renderer* renderer = SDL_CreateRenderer(window,
					      -1,
					      SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC) ;
  //printRenderDriverInfo();

  bool quit = false;

  Model m("model.xyz");
  
  Vec3f camera = {0, 1, -5};
  Vec3f up = {0,1,0};
  Vec3f cameraDir = {0,0,1};

  Matrix4f camera_transform;

  //Matrix4f viewport_and_projection = viewport(width, height) * projection();
  Graphics g(renderer, camera_transform);
  g.addModel(m);
 
  float mouseSpeed = 0.005;
  while(!quit){
    //std::cout<<frame<<'\n';
    lastTick = SDL_GetTicks();
    while(SDL_PollEvent(&event)){
      if(event.type == SDL_QUIT)
	quit = true;
      // TODO: EVENT CODE GOES HERE //
      if(event.type==SDL_MOUSEMOTION){
	float xrel = event.motion.xrel * mouseSpeed;
	/* rotation about y-axis
	 * if camera moves right ... 
	 *
	 */
	cameraDir[2] = cos(xrel)*cameraDir[2]-sin(xrel)*cameraDir[0];
	cameraDir[0] = sin(xrel)*cameraDir[2]+cos(xrel)*cameraDir[0];
      }
      if(event.type == SDL_KEYDOWN){
	switch (event.key.keysym.sym)
	  {
	  case SDLK_w:
	    camera += cameraDir;
	    break;
	  case SDLK_a:
	    camera -= cross(up,cameraDir);
	    break;
	  case SDLK_s:
	    camera -= cameraDir;
	    break;
	  case SDLK_d:
	    camera += cross(up,cameraDir);
	    break;
	  }
      }
      //std::cout<<"Position: "<<camera<<std::endl;
      //std::cout<<"Direction: "<<cameraDir<<std::endl;
    }
    // clear
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    g.clearZBuffer();

    /*** DRAW CODE GOES HERE ***/
    //Matrix4f world_to_screen = viewport_and_projection * cameraLookAt(camera, cameraDir, up);
    camera_transform = cameraLookAt(camera, cameraDir, up);
    
    
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    g.render();

    // present
    SDL_RenderPresent(renderer);
    frame++;

    currentTick = SDL_GetTicks();
    delayTicks = frameTicks - (currentTick - lastTick);
    if (delayTicks < 0) delayTicks = 0;
    SDL_Delay(delayTicks);

      
  }

  
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
