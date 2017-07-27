#include "MaterialPoint_Factory_Classic_CC.h"

// ----------------------------------------------------------------------------
MaterialPoint_BC *MaterialPoint_Factory_Classic_CC::createMaterialPoint(glm::dvec3 d3Center)
{
	MaterialPoint_Classic_CC *thisMP = new MaterialPoint_Classic_CC();

	thisMP->d3_Position = d3Center;

	return(thisMP);
}
// ----------------------------------------------------------------------------

