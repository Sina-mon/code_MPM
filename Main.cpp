#include <iostream>
#include <vector>

#include <omp.h>

#define STB_IMAGE_IMPLEMENTATION

#include "PhysicsEngine.h"
#include "GraphicsEngine.h"
#include "Definitions.h"

#include "ConstitutiveRelation.h"

int main (int argc, char ** argv)
{
	// physics engine initialization ------------------------------------------
	PhysicsEngine thePhysicsEngine;

//	thePhysicsEngine.initializeWorld_Classic_Ring_Fan();
//	thePhysicsEngine.initializeWorld_Classic_Cellular_Langrand();
	thePhysicsEngine.initializeWorld_Classic_Cellular_Langrand_Hexagonal();
//	thePhysicsEngine.initializeWorld_Classic_Foam();

	// graphics engine initialization -----------------------------------------
	GraphicsEngine theGraphicsEngine;

	theGraphicsEngine.runVisualization(&thePhysicsEngine);

	return(0);
}

