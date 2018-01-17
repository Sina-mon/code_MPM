#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::runVisualization(PhysicsEngine *pPhysicsEngine, bool bExit = false)
{
	this->initializeSystems();
	this->setPhysicsEngineReference(pPhysicsEngine);

	this->gameLoop(bExit);
}
// ----------------------------------------------------------------------------
