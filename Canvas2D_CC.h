#ifndef CANVAS2D_CC_H
#define CANVAS2D_CC_H

#include <vector>
#include <iostream>

#include "Definitions.h"
// --------------------------------------------------------
struct Voxel_ST
{
	bool b_Active = false;
	bool b_Surface = false;
	bool b_Load = false;
	bool b_ESO = false;// should/not be included for optimization
	bool b_Redundant = false;// marked as redundant but still active

	double d_Objective = 0.0;
	double d_ESO_Opacity = 1.0;

	unsigned long int u_ID = 0;

	glm::dvec2 d2_Position = glm::dvec2(0.0, 0.0);

	bool operator() (Voxel_ST V1, Voxel_ST V2) { return (V1.d_Objective < V2.d_Objective);}
};
// --------------------------------------------------------
class Canvas2D_CC
{
	public:
		Canvas2D_CC(glm::dvec2 d2Size, double dOffset);
		virtual ~Canvas2D_CC();

		double d_Offset;
		glm::dvec2 d2_Size;
		glm::uvec2 u2_Size;
		std::vector<Voxel_ST *> v_Voxels;

		void drawCircle(glm::dvec2 d2Center, double dRadius);
		void drawRing(glm::dvec2 d2Center, double dRadius_Outer, double dRadius_Inner);
		void drawRectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation);
		void drawPrism(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRatio, double dRotation);

		void cutCircle(glm::dvec2 d2Center, double dRadius);
		void cutRectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation);

		void setESORectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation, bool bFlag);
		void setLoadRectangle(glm::dvec2 d2Center, glm::dvec2 d2Size, double dRotation, bool bFlag);

		std::vector<Voxel_ST> getVoxels(bool bState);

		unsigned long int getCount_Active(bool bState);
		unsigned long int getCount_ESO(bool bState);
		unsigned long int getCount_Redundant(bool bState);

		static unsigned long int getIndex(glm::uvec2 u2Index, glm::uvec2 u2Size)
		{
			unsigned long int index = (unsigned long int)u2Index.x + (unsigned long int)u2Index.y*(u2Size.x);// + (unsigned long int)i3Index.z*(i3Size.x*i3Size.y);

			unsigned long int uNodes = (unsigned long int)u2Size.x * u2Size.y;// * i3Size.z;
			if(0 <= index && index < uNodes)
				return(index);
			else
			{
				std::cout << "Error, in Canvas::getIndex: out of bounds." << std::endl;
				std::cout << "index: " << index << std::endl;
				std::cout << "uNodes: " << uNodes << std::endl;
				return(0);
			}
		}
	protected:
	private:
};

#endif
// --------------------------------------------------------
