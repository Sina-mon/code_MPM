#ifndef GRAPHICSENGINE_H
#define GRAPHICSENGINE_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <SDL2/SDL.h>

#include "Definitions.h"

#include "Errors.h"
#include "GLSLProgram.h"
#include "Vertex.h"
#include "Mesh.h"
#include "Texture.h"
#include "Transformation.h"
#include "Camera.h"
#include "Light.h"

#include "PhysicsEngine.h"
#include "MaterialPoint_BC.h"
#include "GridPoint.h"

enum class GameState {PLAY, EXIT};

class GraphicsEngine
{
	public:
		GraphicsEngine();
		virtual ~GraphicsEngine();

		void initializeSystems(void);
		void setPhysicsEngineReference(PhysicsEngine *mpmPhysicsEngine);
		void drawGame(void);
		void runVisualization(PhysicsEngine *pPhysicsEngine, bool bExit);
		void saveScreenshot(int x, int y, int w, int h, const char * filename);
	protected:
		//snapshot save
		float f_TimeSnapshot_Interval = 1.0; // set inside the setPhysicsEngineReference method
		float f_TimeSnapshot_LastSave = -1.0e12; // before creation
		int i_TimeSnapshotCycle = 0;

		SDL_Window *p_Window = NULL;

		float i_ScreenWidth = 0.4*4*1000;
		float i_ScreenHeight = 0.4*1*1000;

		glm::vec3 f3_World_Center = glm::vec3(0.0, 0.0, 0.0);
		glm::vec3 f3_World_Dimensions = glm::vec3(0.0, 0.0, 0.0);

		glm::vec3 f3_Camera_Position_Original = glm::vec3(0.0,0.0, 0.0);
		glm::vec3 f3_Camera_Target_Original = glm::vec3(0.0,0.0,0.0);

		GameState e_GameState = GameState::PLAY;

		// shaders
		GLSLProgram gl_BasicProgram;
		GLSLProgram gl_ShadowProgram;
		GLSLProgram gl_FinalProgram;

		Texture *gl_Diffuse_Texture = NULL;
		Texture	*gl_Shadow_Texture  = NULL;

		Camera *gl_Camera = NULL;
		Light *gl_Light = NULL;
		Mesh *gl_Particle_Mesh = NULL;


		enum class enum_Canvas : int {
			MAIN = 0,
			SOLID,
			J2_PLASTICSTRAIN,
			J2_STRESS,
			COUNT,
			MASSGRADIENT_GP,
			MASSGRADIENT_MP,
			ENERGY_STRAIN,
			ENERGY_PLASTIC,
		};

		Texture	*v_Canvas_Texture[(int)enum_Canvas::COUNT];
		Mesh	*v_Canvas_Mesh[(int)enum_Canvas::COUNT];

		float f_Time = 0.0;

		PhysicsEngine *mpm_PhysicsEngine;

		void initializeShaders(void);
		void gameLoop(bool bExit);
		void processInput(void);
	private:
};

#endif // GRAPHICSENGINE_H
