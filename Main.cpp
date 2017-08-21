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

//	thePhysicsEngine.initializeWorld_CPDI_HalfRing_Xiang_FullLength();
	thePhysicsEngine.initializeWorld_CPDI_HalfRing_Xiang_PlainStress();
//	thePhysicsEngine.initializeWorld_CPDI_HalfRing_Xiang_PlainStrain();
//	thePhysicsEngine.initializeWorld_Bar_CPDI();
//	thePhysicsEngine.initializeWorld_Bar();
//	thePhysicsEngine.initializeWorld_Ring();
//	thePhysicsEngine.initializeWorld_AuxeticSwisscheeseCell();
//	thePhysicsEngine.initializeWorld_AuxeticPolygonCell();

	// graphics engine initialization -----------------------------------------
	GraphicsEngine theGraphicsEngine;

	theGraphicsEngine.runVisualization(&thePhysicsEngine);

	return(0);
}

