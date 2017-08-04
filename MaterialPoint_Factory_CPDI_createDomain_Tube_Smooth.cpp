#include "MaterialPoint_Factory_CPDI_CC.h"

// ----------------------------------------------------------------------------
std::vector<MaterialPoint_BC *> MaterialPoint_Factory_CPDI_CC::createDomain_Tube_Smooth(glm::dvec3 d3Center, glm::dvec3 d3Rotation, double dRadius_Outer, double dRadius_Inner, double dLength, double dOffset)
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

	double dRadius_Average = 0.5*(dRadius_Outer + dRadius_Inner);
	double dAngle = (1.0*_PI)/90.0;
	double dDepth = dOffset;
	double dRadius_Offset = (dRadius_Outer - dRadius_Inner)/4.0;

	for(double dRadius_Inner_Cell = dRadius_Inner; dRadius_Inner_Cell < dRadius_Outer-0.5*dRadius_Offset; dRadius_Inner_Cell += dRadius_Offset)
	{
		double dRadius_Outer_Cell = dRadius_Inner_Cell + dRadius_Offset;
		for(double dz = dz_Start; dz <= 0.5*dLength; dz += dOffset)
		{
			for(double dTheta = -0.5*_PI; dTheta < 0.5*_PI-dAngle; dTheta += dAngle)
			{
				std::vector<MaterialPoint_BC *> vMaterialPoints;

				glm::dvec3 d3Center_Layer = glm::dvec3(0.0,0.0,0.0);
				d3Center_Layer.z = dz;

				std::vector<glm::dvec3> vCorner;
				vCorner.resize(8);
				{
					vCorner[0] = d3Center_Layer + glm::dvec3(dRadius_Inner_Cell*glm::cos(dTheta+dAngle)	,dRadius_Inner_Cell*glm::sin(dTheta+dAngle)	,-0.5*dDepth);
					vCorner[1] = d3Center_Layer + glm::dvec3(dRadius_Inner_Cell*glm::cos(dTheta+dAngle)	,dRadius_Inner_Cell*glm::sin(dTheta+dAngle)	,+0.5*dDepth);
					vCorner[2] = d3Center_Layer + glm::dvec3(dRadius_Inner_Cell*glm::cos(dTheta)		,dRadius_Inner_Cell*glm::sin(dTheta)		,+0.5*dDepth);
					vCorner[3] = d3Center_Layer + glm::dvec3(dRadius_Inner_Cell*glm::cos(dTheta)		,dRadius_Inner_Cell*glm::sin(dTheta)		,-0.5*dDepth);
					vCorner[4] = d3Center_Layer + glm::dvec3(dRadius_Outer_Cell*glm::cos(dTheta+dAngle)	,dRadius_Outer_Cell*glm::sin(dTheta+dAngle)	,-0.5*dDepth);
					vCorner[5] = d3Center_Layer + glm::dvec3(dRadius_Outer_Cell*glm::cos(dTheta+dAngle)	,dRadius_Outer_Cell*glm::sin(dTheta+dAngle)	,+0.5*dDepth);
					vCorner[6] = d3Center_Layer + glm::dvec3(dRadius_Outer_Cell*glm::cos(dTheta)		,dRadius_Outer_Cell*glm::sin(dTheta)		,+0.5*dDepth);
					vCorner[7] = d3Center_Layer + glm::dvec3(dRadius_Outer_Cell*glm::cos(dTheta)		,dRadius_Outer_Cell*glm::sin(dTheta)		,-0.5*dDepth);
				}
				vMaterialPoints = this->createCell_Hexahedron(vCorner);
				allMaterialPoint.insert(allMaterialPoint.end(), vMaterialPoints.begin(), vMaterialPoints.end());

				if(bSingleLayer == false)
				{
					d3Center_Layer.z = -dz;
					std::vector<glm::dvec3> vCorner;
					vCorner.resize(8);
					{
						vCorner[0] = d3Center_Layer + glm::dvec3(dRadius_Inner_Cell*glm::cos(dTheta+dAngle)	,dRadius_Inner_Cell*glm::sin(dTheta+dAngle)	,-0.5*dDepth);
						vCorner[1] = d3Center_Layer + glm::dvec3(dRadius_Inner_Cell*glm::cos(dTheta+dAngle)	,dRadius_Inner_Cell*glm::sin(dTheta+dAngle)	,+0.5*dDepth);
						vCorner[2] = d3Center_Layer + glm::dvec3(dRadius_Inner_Cell*glm::cos(dTheta)		,dRadius_Inner_Cell*glm::sin(dTheta)		,+0.5*dDepth);
						vCorner[3] = d3Center_Layer + glm::dvec3(dRadius_Inner_Cell*glm::cos(dTheta)		,dRadius_Inner_Cell*glm::sin(dTheta)		,-0.5*dDepth);
						vCorner[4] = d3Center_Layer + glm::dvec3(dRadius_Outer_Cell*glm::cos(dTheta+dAngle)	,dRadius_Outer_Cell*glm::sin(dTheta+dAngle)	,-0.5*dDepth);
						vCorner[5] = d3Center_Layer + glm::dvec3(dRadius_Outer_Cell*glm::cos(dTheta+dAngle)	,dRadius_Outer_Cell*glm::sin(dTheta+dAngle)	,+0.5*dDepth);
						vCorner[6] = d3Center_Layer + glm::dvec3(dRadius_Outer_Cell*glm::cos(dTheta)		,dRadius_Outer_Cell*glm::sin(dTheta)		,+0.5*dDepth);
						vCorner[7] = d3Center_Layer + glm::dvec3(dRadius_Outer_Cell*glm::cos(dTheta)		,dRadius_Outer_Cell*glm::sin(dTheta)		,-0.5*dDepth);
					}
					vMaterialPoints = this->createCell_Hexahedron(vCorner);
					allMaterialPoint.insert(allMaterialPoint.end(), vMaterialPoints.begin(), vMaterialPoints.end());
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

