#ifndef CANVAS_CC_H
#define CANVAS_CC_H

#include <vector>
#include <iostream>

#include "Definitions.h"
// --------------------------------------------------------
struct Voxel_ST
{
	bool b_Active = false;

	unsigned int i_ID = 0;

	glm::dvec3 d3_Position = glm::dvec3(0.0, 0.0, 0.0);
};
// --------------------------------------------------------
class Canvas_CC
{
	public:
		Canvas_CC(glm::dvec3 d3Size, double dOffset);
		virtual ~Canvas_CC();

		glm::dvec3 d3_Size;
		glm::dvec3 i3_Size;
		std::vector<Voxel_ST> v_Voxels;

		void drawRing(glm::dvec3 d3Center, double dRadius_Outer, double dRadius_Inner);
		void drawRectangle(glm::dvec3 d3Center, glm::dvec3 d3Size, glm::dvec3 d3Rotation);

		std::vector<Voxel_ST> getVoxels(bool bState);

		static unsigned int getIndex(glm::ivec3 i3Index, glm::ivec3 i3Size)
		{
			int index = i3Index.x + i3Index.y*(i3Size.x) + i3Index.z*(i3Size.x*i3Size.y);

			int iNodes = i3Size.x * i3Size.y * i3Size.z;
			if(0 <= index && index < iNodes)
				return(index);
			else
			{
				std::cout << "Error, in Canvas::getIndex: out of bounds." << std::endl;
				return(0);
			}
		}
	protected:
	private:
};

#endif
// --------------------------------------------------------
