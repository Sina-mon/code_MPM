#ifndef MATERIALPOINT_CPDI_CC_H
#define MATERIALPOINT_CPDI_CC_H

#include "Definitions.h"
#include "MaterialPoint_BC.h"

class MaterialPoint_CPDI_CC : public MaterialPoint_BC
{
	public:
		glm::dvec3 d3_Corner[4] = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 0.0)};
	protected:
	private:
};

#endif
