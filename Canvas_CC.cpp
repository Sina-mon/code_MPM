#include "Canvas_CC.h"


// --------------------------------------------------------
Canvas_CC::Canvas_CC(glm::dvec3 d3Size, double dOffset)
{
	d3_Size = d3Size;
	i3_Size = glm::ivec3(glm::floor(d3_Size / dOffset));

//	v_Voxels.resize(i3_Size.x * i3_Size.y * i3_Size.z);

	v_Voxels.clear();

	glm::ivec3 i3Index = glm::ivec3(0,0,0);
	for(i3Index.x = 0; i3Index.x < i3_Size.x; i3Index.x++)
	{
		for(i3Index.y = 0; i3Index.y < i3_Size.y; i3Index.y++)
		{
			for(i3Index.z = 0; i3Index.z < i3_Size.z; i3Index.z++)
			{
				unsigned int iIndex = getIndex(i3Index,i3_Size);

				Voxel_ST newVoxel;
				newVoxel.b_Active = false;
				newVoxel.i_ID = iIndex;
				newVoxel.d3_Position = glm::dvec3(i3Index) * dOffset;

				v_Voxels.push_back(newVoxel);
//				v_Voxels[iIndex].b_Active = false;
//				v_Voxels[iIndex].i_ID = iIndex;
//				v_Voxels[iIndex].d3_Position = d3Position;
			}
		}
	}
	std::cout << "Canvas voxels: " << Script(v_Voxels.size()) << std::endl;
}
// --------------------------------------------------------
Canvas_CC::~Canvas_CC()
{
	//dtor
}
// --------------------------------------------------------
void Canvas_CC::drawRing(glm::dvec3 d3Center, double dRadius_Outer, double dRadius_Inner)
{
	d3Center.z = 0.0;
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dvec3 d3Position = v_Voxels[index_Voxel].d3_Position;

		double dDistance = glm::length(d3Position - d3Center);

		if(dRadius_Inner <= dDistance &&  dDistance <= dRadius_Outer)
		{
			v_Voxels[index_Voxel].b_Active = true;
		}
	}
}
// --------------------------------------------------------
void Canvas_CC::drawRectangle(glm::dvec3 d3Center, glm::dvec3 d3Size, glm::dvec3 d3Rotation)
{
	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
		glm::dmat4 m4Transformation_Position = glm::translate(-d3Center);
		glm::dmat4 m4Transformation_RotationX = glm::rotate(-d3Rotation.x, glm::dvec3(1.0, 0.0, 0.0));
		glm::dmat4 m4Transformation_RotationY = glm::rotate(-d3Rotation.y, glm::dvec3(0.0, 1.0, 0.0));
		glm::dmat4 m4Transformation_RotationZ = glm::rotate(-d3Rotation.z, glm::dvec3(0.0, 0.0, 1.0));

		glm::dmat4 m4Transformation_Combined = glm::dmat4(1.0); // sina, apparently, glm::mat4() and glm::mat4(1.0) build identity matrices
		m4Transformation_Combined *= m4Transformation_RotationZ;
		m4Transformation_Combined *= m4Transformation_RotationY;
		m4Transformation_Combined *= m4Transformation_RotationX;
		m4Transformation_Combined *= m4Transformation_Position;

		glm::dvec3 d3Position_Local = glm::dvec3(m4Transformation_Combined * glm::dvec4(v_Voxels[index_Voxel].d3_Position, 1.0));

		if(d3Position_Local.x < -0.5*d3Size.x || +0.5*d3Size.x < d3Position_Local.x)
			continue;
		if(d3Position_Local.y < -0.5*d3Size.y || +0.5*d3Size.y < d3Position_Local.y)
			continue;
		if(d3Position_Local.z < -0.5*d3Size.z || +0.5*d3Size.z < d3Position_Local.z)
			continue;

		v_Voxels[index_Voxel].b_Active = true;
	}
}
// --------------------------------------------------------
std::vector<Voxel_ST> Canvas_CC::getVoxels(bool bState = true)
{
	std::vector<Voxel_ST> vResult;
	vResult.clear();

	for(unsigned int index_Voxel = 0; index_Voxel < v_Voxels.size(); index_Voxel++)
	{
	//std::cout << "here" << std::endl;
		if(v_Voxels[index_Voxel].b_Active == bState)
			vResult.push_back(v_Voxels[index_Voxel]);
	}

	return(vResult);
}
// --------------------------------------------------------
