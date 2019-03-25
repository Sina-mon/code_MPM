#include "Canvas2D_CC.h"


// --------------------------------------------------------
Canvas2D_CC::Canvas2D_CC(glm::dvec2 d2Size, double dOffset)
{
	std::cout << "max_size: " << v_Voxels.max_size() << "\n";
	d_Offset = dOffset;
	d2_Size = d2Size;
	u2_Size = glm::floor(d2_Size / dOffset);

	v_Voxels.clear();
	v_Voxels.resize(u2_Size.x * u2_Size.y);

	glm::uvec2 u2Index = glm::uvec2(0,0);
	for(u2Index.y = 0; u2Index.y < u2_Size.y; u2Index.y++)
	{
		for(u2Index.x = 0; u2Index.x < u2_Size.x; u2Index.x++)
		{
			unsigned long int uIndex = (unsigned long int)getIndex(u2Index,u2_Size);

			Voxel_ST *newVoxel = new Voxel_ST;

			newVoxel->b_Active = false;
			newVoxel->u_ID = uIndex;
			newVoxel->d2_Position = glm::dvec2(u2Index) * dOffset + 0.5*glm::dvec2(dOffset,dOffset);

			v_Voxels[uIndex] = newVoxel;
		}
	}
	std::cout << "Canvas voxels: " << Script(v_Voxels.size()) << std::endl;
}
// --------------------------------------------------------
Canvas2D_CC::~Canvas2D_CC()
{
	for(unsigned long int iIndex = 0; iIndex < v_Voxels.size(); iIndex++)
		delete v_Voxels[iIndex];
}
// --------------------------------------------------------
void Canvas2D_CC::drawCircle(glm::dvec2 d2Center, double dRadius)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dvec2 d2Position = v_Voxels[index_Voxel]->d2_Position;

		double dDistance = glm::length(d2Position - d2Center);

		if(dDistance > dRadius)
			continue;

		// if processes reaches this point
		v_Voxels[index_Voxel]->b_Active = true;
	}
}
// --------------------------------------------------------
void Canvas2D_CC::drawRing(glm::dvec2 d2Center, double dRadius_Outer, double dRadius_Inner)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dvec2 d2Position = v_Voxels[index_Voxel]->d2_Position;

		double dDistance = glm::length(d2Position - d2Center);

		if(dDistance < dRadius_Inner + 0.5*d_Offset || dRadius_Outer - 0.5*d_Offset < dDistance)
			continue;

		bool bSurface_Local = false;
		if(glm::abs(dDistance - dRadius_Outer) < d_Offset)
			bSurface_Local = true;
		if(glm::abs(dDistance - dRadius_Inner) < d_Offset)
			bSurface_Local = true;

		if(bSurface_Local == true)
		{
			if(v_Voxels[index_Voxel]->b_Active == true && v_Voxels[index_Voxel]->b_Surface == true)
				bSurface_Local = true;
			if(v_Voxels[index_Voxel]->b_Active == true && v_Voxels[index_Voxel]->b_Surface == false)
				bSurface_Local = false;
		}

		v_Voxels[index_Voxel]->b_Surface = bSurface_Local;
		v_Voxels[index_Voxel]->b_Active = true;
	}
}
// --------------------------------------------------------
void Canvas2D_CC::drawRectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dmat4 m4Transformation_Position = glm::translate(glm::dvec3(-d2Center, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(-dRotation, glm::dvec3(0.0, 0.0, 1.0));

		glm::dmat4 m4Transformation_Combined = glm::dmat4(1.0); // sina, apparently, glm::mat4() and glm::mat4(1.0) build identity matrices
		m4Transformation_Combined *= m4Transformation_RotationZ;
//		m4Transformation_Combined *= m4Transformation_RotationY;
//		m4Transformation_Combined *= m4Transformation_RotationX;
		m4Transformation_Combined *= m4Transformation_Position;

		glm::dvec2 d2Position_Local = glm::dvec2(m4Transformation_Combined * glm::dvec4(v_Voxels[index_Voxel]->d2_Position, 0.0, 1.0));

		if(d2Position_Local.x < -0.5*d2Size.x || +0.5*d2Size.x < d2Position_Local.x)
			continue;
		if(d2Position_Local.y < -0.5*d2Size.y || +0.5*d2Size.y < d2Position_Local.y)
			continue;

		v_Voxels[index_Voxel]->b_Active = true;
	}
}
// --------------------------------------------------------
void Canvas2D_CC::drawPrism(glm::dvec2 d2Center, glm::dvec2 d2Size, double dExtra, double dRotation)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dmat4 m4Transformation_Position = glm::translate(glm::dvec3(-d2Center, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(-dRotation, glm::dvec3(0.0, 0.0, 1.0));

		glm::dmat4 m4Transformation_Combined = glm::dmat4(1.0); // sina, apparently, glm::mat4() and glm::mat4(1.0) build identity matrices
		m4Transformation_Combined *= m4Transformation_RotationZ;
//		m4Transformation_Combined *= m4Transformation_RotationY;
//		m4Transformation_Combined *= m4Transformation_RotationX;
		m4Transformation_Combined *= m4Transformation_Position;

		glm::dvec2 d2Position_Local = glm::dvec2(m4Transformation_Combined * glm::dvec4(v_Voxels[index_Voxel]->d2_Position, 0.0, 1.0));

		if(d2Position_Local.x < -0.5*d2Size.x || +0.5*d2Size.x < d2Position_Local.x)
			continue;
		if(d2Position_Local.y < -0.5*d2Size.y+0.5*dExtra*(2.0*d2Position_Local.x/d2Size.x) || +0.5*d2Size.y-0.5*dExtra*(2.0*d2Position_Local.x/d2Size.x) < d2Position_Local.y)
			continue;

		v_Voxels[index_Voxel]->b_Active = true;
	}
}
// --------------------------------------------------------
void Canvas2D_CC::cutCircle(glm::dvec2 d2Center, double dRadius)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dvec2 d2Position = v_Voxels[index_Voxel]->d2_Position;

		double dDistance = glm::length(d2Position - d2Center);

		if(dDistance > dRadius)
			continue;

		v_Voxels[index_Voxel]->b_Active = false;
	}
}
// --------------------------------------------------------
void Canvas2D_CC::cutRectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dmat4 m4Transformation_Position = glm::translate(glm::dvec3(-d2Center, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(-dRotation, glm::dvec3(0.0, 0.0, 1.0));

		glm::dmat4 m4Transformation_Combined = glm::dmat4(1.0); // sina, apparently, glm::mat4() and glm::mat4(1.0) build identity matrices
		m4Transformation_Combined *= m4Transformation_RotationZ;
//		m4Transformation_Combined *= m4Transformation_RotationY;
//		m4Transformation_Combined *= m4Transformation_RotationX;
		m4Transformation_Combined *= m4Transformation_Position;

		glm::dvec2 d2Position_Local = glm::dvec2(m4Transformation_Combined * glm::dvec4(v_Voxels[index_Voxel]->d2_Position, 0.0, 1.0));

		if(d2Position_Local.x < -0.5*d2Size.x || +0.5*d2Size.x < d2Position_Local.x)
			continue;
		if(d2Position_Local.y < -0.5*d2Size.y || +0.5*d2Size.y < d2Position_Local.y)
			continue;

		v_Voxels[index_Voxel]->b_Active = false;
	}
}
// --------------------------------------------------------
void Canvas2D_CC::setESORectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation, bool bFlag)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dmat4 m4Transformation_Position = glm::translate(glm::dvec3(-d2Center, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(-dRotation, glm::dvec3(0.0, 0.0, 1.0));

		glm::dmat4 m4Transformation_Combined = glm::dmat4(1.0); // sina, apparently, glm::mat4() and glm::mat4(1.0) build identity matrices
		m4Transformation_Combined *= m4Transformation_RotationZ;
//		m4Transformation_Combined *= m4Transformation_RotationY;
//		m4Transformation_Combined *= m4Transformation_RotationX;
		m4Transformation_Combined *= m4Transformation_Position;

		glm::dvec2 d2Position_Local = glm::dvec2(m4Transformation_Combined * glm::dvec4(v_Voxels[index_Voxel]->d2_Position, 0.0, 1.0));

		if(d2Position_Local.x < -0.5*d2Size.x || +0.5*d2Size.x < d2Position_Local.x)
			continue;
//		if(d2Position_Local.y < -0.5*d2Size.y || +0.5*d2Size.y < d2Position_Local.y)
		if(glm::abs(d2Position_Local.y) > 0.5*d2Size.y)
			continue;

		v_Voxels[index_Voxel]->b_ESO = bFlag;
	}
}
// --------------------------------------------------------
void Canvas2D_CC::setLoadRectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation, bool bFlag)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dmat4 m4Transformation_Position = glm::translate(glm::dvec3(-d2Center, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(-dRotation, glm::dvec3(0.0, 0.0, 1.0));

		glm::dmat4 m4Transformation_Combined = glm::dmat4(1.0); // sina, apparently, glm::mat4() and glm::mat4(1.0) build identity matrices
		m4Transformation_Combined *= m4Transformation_RotationZ;
//		m4Transformation_Combined *= m4Transformation_RotationY;
//		m4Transformation_Combined *= m4Transformation_RotationX;
		m4Transformation_Combined *= m4Transformation_Position;

		glm::dvec2 d2Position_Local = glm::dvec2(m4Transformation_Combined * glm::dvec4(v_Voxels[index_Voxel]->d2_Position, 0.0, 1.0));

		if(d2Position_Local.x < -0.5*d2Size.x || +0.5*d2Size.x < d2Position_Local.x)
			continue;
//		if(d2Position_Local.y < -0.5*d2Size.y || +0.5*d2Size.y < d2Position_Local.y)
		if(glm::abs(d2Position_Local.y) > 0.5*d2Size.y)
			continue;

		v_Voxels[index_Voxel]->b_Load = bFlag;
	}
}
// --------------------------------------------------------
void Canvas2D_CC::setSupportRectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation, bool bFlag)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dmat4 m4Transformation_Position = glm::translate(glm::dvec3(-d2Center, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(-dRotation, glm::dvec3(0.0, 0.0, 1.0));

		glm::dmat4 m4Transformation_Combined = glm::dmat4(1.0); // sina, apparently, glm::mat4() and glm::mat4(1.0) build identity matrices
		m4Transformation_Combined *= m4Transformation_RotationZ;
//		m4Transformation_Combined *= m4Transformation_RotationY;
//		m4Transformation_Combined *= m4Transformation_RotationX;
		m4Transformation_Combined *= m4Transformation_Position;

		glm::dvec2 d2Position_Local = glm::dvec2(m4Transformation_Combined * glm::dvec4(v_Voxels[index_Voxel]->d2_Position, 0.0, 1.0));

		if(d2Position_Local.x < -0.5*d2Size.x || +0.5*d2Size.x < d2Position_Local.x)
			continue;
//		if(d2Position_Local.y < -0.5*d2Size.y || +0.5*d2Size.y < d2Position_Local.y)
		if(glm::abs(d2Position_Local.y) > 0.5*d2Size.y)
			continue;

		v_Voxels[index_Voxel]->b_Support = bFlag;
	}
}
// --------------------------------------------------------
std::vector<Voxel_ST *> Canvas2D_CC::getVoxels_Active(bool bState = true)
{
	std::vector<Voxel_ST *> vResult;
	vResult.clear();

	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
	//std::cout << "here" << std::endl;
		if(v_Voxels[index_Voxel]->b_Active == bState)
			vResult.push_back(v_Voxels[index_Voxel]);
	}

	return(vResult);
}
// --------------------------------------------------------
std::vector<Voxel_ST *> Canvas2D_CC::getVoxels_ESO(bool bState = true)
{
	std::vector<Voxel_ST *> vResult;
	vResult.clear();

	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
	//std::cout << "here" << std::endl;
		if(v_Voxels[index_Voxel]->b_ESO == bState)
			vResult.push_back(v_Voxels[index_Voxel]);
	}

	return(vResult);
}
// --------------------------------------------------------
unsigned long int Canvas2D_CC::getCount_Active(bool bState = true)
{
	unsigned long int iResult = 0;

	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		if(v_Voxels[index_Voxel]->b_Active == bState)
			iResult++;
	}

	return(iResult);
}
// --------------------------------------------------------
unsigned long int Canvas2D_CC::getCount_ESO(bool bState = true)
{
	unsigned long int iResult = 0;

	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		if(v_Voxels[index_Voxel]->b_ESO == bState)
			iResult++;
	}

	return(iResult);
}
// --------------------------------------------------------
unsigned long int Canvas2D_CC::getCount_Redundant(bool bState = true)
{
	unsigned long int iResult = 0;

	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		if(v_Voxels[index_Voxel]->b_Redundant == bState)
			iResult++;
	}

	return(iResult);
}
// --------------------------------------------------------
void Canvas2D_CC::filterObjective_Smooth(double dRadius)
{
	std::vector<Voxel_ST *> vVoxels_Active = this->getVoxels_Active(true);
	unsigned long int iVoxels_Active = vVoxels_Active.size();
	std::vector<Voxel_ST> vVoxels_Active_Clone;
	vVoxels_Active_Clone.resize(iVoxels_Active);
	for(unsigned long int index = 0; index < iVoxels_Active; index++)
	{
		vVoxels_Active_Clone[index].u_ID = vVoxels_Active[index]->u_ID;
		vVoxels_Active_Clone[index].d2_Position = vVoxels_Active[index]->d2_Position;
		vVoxels_Active_Clone[index].d_Objective = 0.0;

		int iVoxels_Neighbour = 0;
		double dIntegral_Weight = 0.0;
		for(unsigned long int index_neighbour = 0; index_neighbour < iVoxels_Active; index_neighbour++)
		{
			double dDistance = glm::length(vVoxels_Active[index]->d2_Position - vVoxels_Active[index_neighbour]->d2_Position);

			if(vVoxels_Active[index_neighbour]->b_ESO == false)
				continue;

			if(dDistance < dRadius)
			{
//				iVoxels_Neighbour++;
//				vVoxels_Active_Clone[index].d_Objective += vVoxels_Active[index_neighbour]->d_Objective;

				double dWeight = (dRadius - dDistance)/dRadius;
				dIntegral_Weight += dWeight;
				vVoxels_Active_Clone[index].d_Objective += dWeight * vVoxels_Active[index_neighbour]->d_Objective;
			}
		}
//		vVoxels_Active_Clone[index].d_Objective /= iVoxels_Neighbour;
		vVoxels_Active_Clone[index].d_Objective /= dIntegral_Weight;
	}
	for(unsigned long int index = 0; index < iVoxels_Active; index++)
	{
		vVoxels_Active[index]->d_Objective = vVoxels_Active_Clone[index].d_Objective;
	}
}
// --------------------------------------------------------
double Canvas2D_CC::getObjective_Sum(void)
{
	double dResult = 0.0;

	std::vector<Voxel_ST *> vVoxels_Active = this->getVoxels_Active(true);
	unsigned long int iVoxels_Active = vVoxels_Active.size();


	for(unsigned long int index = 0; index < iVoxels_Active; index++)
	{
		dResult += vVoxels_Active[index]->d_Objective;
	}

	return(dResult);
}
