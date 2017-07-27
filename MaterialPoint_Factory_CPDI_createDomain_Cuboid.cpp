#include "MaterialPoint_Factory_CPDI_CC.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint_BC *> MaterialPoint_Factory_CPDI_CC::createDomain_Cuboid(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset)
{
	std::vector<MaterialPoint_BC *> allMaterialPoint;

	// if the z-dimension is set to zero, create only one layer of material points
	// otherwise create even number of layers that are symmetric relative to the x-y plane
	bool bSingleLayer = false;
	double dz_Start = 0.5*dOffset;
	if(d3Dimension.z <= dOffset)
	{
		bSingleLayer = true;
		dz_Start = 0.0;
	}

	for(double dx = 0.5*dOffset; dx <= 0.5*d3Dimension.x; dx += dOffset)
	{//create a quarter
		for(double dy = 0.5*dOffset; dy <= 0.5*d3Dimension.y; dy += dOffset)
		{
			for(double dz = dz_Start; dz <= 0.5*d3Dimension.z; dz += dOffset)
			{
				std::vector<MaterialPoint_BC *> vMP;
//
//				vMaterialPoints = createMaterialPoint_T4(d3Center + glm::dvec3(+dx, +dy, +dz), dOffset);
//				allMaterialPoint.insert(allMaterialPoint.end(), vMaterialPoints.begin(), vMaterialPoints.end());

//				MaterialPoint_BC *thisMaterialPoint;
//				allMaterialPoint.push_back(thisMaterialPoint);

				vMP = this->createMaterialPoint(d3Center + glm::dvec3(+dx, +dy, +dz), dOffset);
				allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

				vMP = this->createMaterialPoint(d3Center + glm::dvec3(-dx, +dy, +dz), dOffset);
				allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

				vMP = this->createMaterialPoint(d3Center + glm::dvec3(-dx, -dy, +dz), dOffset);
				allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

				vMP = this->createMaterialPoint(d3Center + glm::dvec3(+dx, -dy, +dz), dOffset);
				allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

				if(bSingleLayer == false)
				{
					vMP = this->createMaterialPoint(d3Center + glm::dvec3(+dx, +dy, -dz), dOffset);
					allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

					vMP = this->createMaterialPoint(d3Center + glm::dvec3(-dx, +dy, -dz), dOffset);
					allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

					vMP = this->createMaterialPoint(d3Center + glm::dvec3(-dx, -dy, -dz), dOffset);
					allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

					vMP = this->createMaterialPoint(d3Center + glm::dvec3(+dx, -dy, -dz), dOffset);
					allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());
				}
			}
		}
	}

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

