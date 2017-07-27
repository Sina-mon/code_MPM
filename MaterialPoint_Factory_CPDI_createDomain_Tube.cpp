#include "MaterialPoint_Factory_CPDI_CC.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint_BC *> MaterialPoint_Factory_CPDI_CC::createDomain_Tube(glm::dvec3 d3Center, glm::dvec3 d3Rotation, double dRadius_Outer, double dRadius_Inner, double dLength, double dOffset)
{
	std::vector<MaterialPoint_BC *> allMaterialPoint;

	// if the z-dimension is set to zero, create only one layer of material points
	// otherwise create even number of layers that are symmetric relative to the x-y plane
	bool bSingleLayer = false;
	double dz_Start = 0.5*dOffset;
	if(dLength <= dOffset)
	{
		bSingleLayer = true;
		dz_Start = 0.0;
	}

	for(double dx = 0.5*dOffset; dx <= dRadius_Outer; dx += dOffset)
	{//create a quarter
		for(double dy = 0.5*dOffset; dy <= dRadius_Outer; dy += dOffset)
		{
			for(double dz = dz_Start; dz <= 0.5*dLength; dz += dOffset)
			{
				double dRadialDistance = glm::sqrt(dx*dx + dy*dy);
				if(dRadius_Inner*dRadius_Inner < dx*dx + dy*dy && dx*dx + dy*dy < dRadius_Outer*dRadius_Outer)
				{
					std::vector<MaterialPoint_BC *> vMP;

					vMP = this->createCell_Hexahedron(glm::dvec3(+dx, +dy, +dz), dOffset);
					allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

					vMP = this->createCell_Hexahedron(glm::dvec3(-dx, +dy, +dz), dOffset);
					allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

					vMP = this->createCell_Hexahedron(glm::dvec3(-dx, -dy, +dz), dOffset);
					allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

					vMP = this->createCell_Hexahedron(glm::dvec3(+dx, -dy, +dz), dOffset);
					allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

					if(bSingleLayer == false)
					{
						vMP = this->createCell_Hexahedron(glm::dvec3(+dx, +dy, -dz), dOffset);
						allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

						vMP = this->createCell_Hexahedron(glm::dvec3(-dx, +dy, -dz), dOffset);
						allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

						vMP = this->createCell_Hexahedron(glm::dvec3(-dx, -dy, -dz), dOffset);
						allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());

						vMP = this->createCell_Hexahedron(glm::dvec3(+dx, -dy, -dz), dOffset);
						allMaterialPoint.insert(allMaterialPoint.end(), vMP.begin(), vMP.end());
					}
				}
			}
		}
	}

	for(unsigned int index = 0; index < allMaterialPoint.size(); index++)
	{
		MaterialPoint_CPDI_CC *thisMaterialPoint = (MaterialPoint_CPDI_CC *)allMaterialPoint[index];

		glm::dmat4 m4Transformation_Position = glm::translate(d3Center);
		glm::dmat4 m4Transformation_RotationX = glm::rotate(d3Rotation.x, glm::dvec3(1.0, 0.0, 0.0));
		glm::dmat4 m4Transformation_RotationY = glm::rotate(d3Rotation.y, glm::dvec3(0.0, 1.0, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(d3Rotation.z, glm::dvec3(0.0, 0.0, 1.0));

		glm::dmat4 m4Transformation_Combined;
		m4Transformation_Combined *= m4Transformation_Position;
		m4Transformation_Combined *= m4Transformation_RotationZ;
		m4Transformation_Combined *= m4Transformation_RotationY;
		m4Transformation_Combined *= m4Transformation_RotationX;

		glm::dvec4 d4Position = glm::dvec4(thisMaterialPoint->d3_Position, 1.0);

		thisMaterialPoint->d3_Position = glm::dvec3(m4Transformation_Combined * d4Position);

		for(int index_Corner = 0; index_Corner < 4; index_Corner++)
		{
			glm::dvec4 d4Position = glm::dvec4(thisMaterialPoint->a_Corner[index_Corner].d3_Position, 1.0);

			thisMaterialPoint->a_Corner[index_Corner].d3_Position = glm::dvec3(m4Transformation_Combined * d4Position);
		}
	}

	return(allMaterialPoint);
}
// ----------------------------------------------------------------------------

