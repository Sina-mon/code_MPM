#include "GraphicsEngine.h"

GraphicsEngine::GraphicsEngine()
{
}
// ----------------------------------------------------------------------------
GraphicsEngine::~GraphicsEngine()
{
	if(gl_Camera != NULL)
		delete gl_Camera;

	if(gl_Light != NULL)
		delete gl_Light;

	if(gl_Diffuse_Texture != NULL)
		delete gl_Diffuse_Texture;

	if(gl_Shadow_Texture != NULL)
		delete gl_Shadow_Texture;

	for(int index = 0; index < (int)enum_Canvas::COUNT; index++)
	{
		if(v_Canvas_Texture[index] != NULL)
			delete v_Canvas_Texture[index];

		if(v_Canvas_Mesh[index] != NULL)
			delete v_Canvas_Mesh[index];
	}

	if(gl_Particle_Mesh != NULL)
		delete gl_Particle_Mesh;

	SDL_DestroyWindow(p_Window);
	SDL_Quit();
}
// ----------------------------------------------------------------------------
