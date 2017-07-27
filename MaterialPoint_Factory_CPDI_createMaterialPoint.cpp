#include "MaterialPoint_Factory_CPDI_CC.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint_BC *> MaterialPoint_Factory_CPDI_CC::createMaterialPoint(glm::dvec3 d3Center, double dOffset)
{
	// formulation from, http://www.baumanneduard.ch/Splitting%20a%20cube%20in%20tetrahedras2.htm
	std::vector<MaterialPoint_BC *> vMaterialPoints;

	dOffset *= 1.0;

	glm::dvec3 d3Cube[8];
	{
		d3Cube[0] = d3Center + glm::dvec3(-0.5*dOffset,-0.5*dOffset,-0.5*dOffset);
		d3Cube[1] = d3Center + glm::dvec3(-0.5*dOffset,-0.5*dOffset,+0.5*dOffset);
		d3Cube[2] = d3Center + glm::dvec3(+0.5*dOffset,-0.5*dOffset,+0.5*dOffset);
		d3Cube[3] = d3Center + glm::dvec3(+0.5*dOffset,-0.5*dOffset,-0.5*dOffset);
		d3Cube[4] = d3Center + glm::dvec3(-0.5*dOffset,+0.5*dOffset,-0.5*dOffset);
		d3Cube[5] = d3Center + glm::dvec3(-0.5*dOffset,+0.5*dOffset,+0.5*dOffset);
		d3Cube[6] = d3Center + glm::dvec3(+0.5*dOffset,+0.5*dOffset,+0.5*dOffset);
		d3Cube[7] = d3Center + glm::dvec3(+0.5*dOffset,+0.5*dOffset,-0.5*dOffset);
	}

	{
		MaterialPoint_CPDI_CC *thisMP = new MaterialPoint_CPDI_CC();
		thisMP->a_Corner[0].d3_Position = d3Cube[0];
		thisMP->a_Corner[1].d3_Position = d3Cube[1];
		thisMP->a_Corner[2].d3_Position = d3Cube[3];
		thisMP->a_Corner[3].d3_Position = d3Cube[7];
		thisMP->d3_Position = 0.25*(thisMP->a_Corner[0].d3_Position+thisMP->a_Corner[1].d3_Position+thisMP->a_Corner[2].d3_Position+thisMP->a_Corner[3].d3_Position);
		vMaterialPoints.push_back(thisMP);
	}
	{
		MaterialPoint_CPDI_CC *thisMP = new MaterialPoint_CPDI_CC();
		thisMP->a_Corner[0].d3_Position = d3Cube[0];
		thisMP->a_Corner[1].d3_Position = d3Cube[4];
		thisMP->a_Corner[2].d3_Position = d3Cube[1];
		thisMP->a_Corner[3].d3_Position = d3Cube[7];
		thisMP->d3_Position = 0.25*(thisMP->a_Corner[0].d3_Position+thisMP->a_Corner[1].d3_Position+thisMP->a_Corner[2].d3_Position+thisMP->a_Corner[3].d3_Position);
		vMaterialPoints.push_back(thisMP);
	}
	{
		MaterialPoint_CPDI_CC *thisMP = new MaterialPoint_CPDI_CC();
		thisMP->a_Corner[0].d3_Position = d3Cube[1];
		thisMP->a_Corner[1].d3_Position = d3Cube[2];
		thisMP->a_Corner[2].d3_Position = d3Cube[3];
		thisMP->a_Corner[3].d3_Position = d3Cube[7];
		thisMP->d3_Position = 0.25*(thisMP->a_Corner[0].d3_Position+thisMP->a_Corner[1].d3_Position+thisMP->a_Corner[2].d3_Position+thisMP->a_Corner[3].d3_Position);
		vMaterialPoints.push_back(thisMP);
	}
	{
		MaterialPoint_CPDI_CC *thisMP = new MaterialPoint_CPDI_CC();
		thisMP->a_Corner[0].d3_Position = d3Cube[1];
		thisMP->a_Corner[1].d3_Position = d3Cube[6];
		thisMP->a_Corner[2].d3_Position = d3Cube[2];
		thisMP->a_Corner[3].d3_Position = d3Cube[7];
		thisMP->d3_Position = 0.25*(thisMP->a_Corner[0].d3_Position+thisMP->a_Corner[1].d3_Position+thisMP->a_Corner[2].d3_Position+thisMP->a_Corner[3].d3_Position);
		vMaterialPoints.push_back(thisMP);
	}
	{
		MaterialPoint_CPDI_CC *thisMP = new MaterialPoint_CPDI_CC();
		thisMP->a_Corner[0].d3_Position = d3Cube[1];
		thisMP->a_Corner[1].d3_Position = d3Cube[4];
		thisMP->a_Corner[2].d3_Position = d3Cube[5];
		thisMP->a_Corner[3].d3_Position = d3Cube[7];
		thisMP->d3_Position = 0.25*(thisMP->a_Corner[0].d3_Position+thisMP->a_Corner[1].d3_Position+thisMP->a_Corner[2].d3_Position+thisMP->a_Corner[3].d3_Position);
		vMaterialPoints.push_back(thisMP);
	}
	{
		MaterialPoint_CPDI_CC *thisMP = new MaterialPoint_CPDI_CC();
		thisMP->a_Corner[0].d3_Position = d3Cube[1];
		thisMP->a_Corner[1].d3_Position = d3Cube[5];
		thisMP->a_Corner[2].d3_Position = d3Cube[6];
		thisMP->a_Corner[3].d3_Position = d3Cube[7];
		thisMP->d3_Position = 0.25*(thisMP->a_Corner[0].d3_Position+thisMP->a_Corner[1].d3_Position+thisMP->a_Corner[2].d3_Position+thisMP->a_Corner[3].d3_Position);
		vMaterialPoints.push_back(thisMP);
	}

	return(vMaterialPoints);
}
// ----------------------------------------------------------------------------
