#ifndef MATERIALPOINT_CPDI_CC_H
#define MATERIALPOINT_CPDI_CC_H

#include "Definitions.h"
#include "MaterialPoint_BC.h"

struct Corner_ST
{
	glm::dvec3 d3_Position = glm::dvec3(0.0,0.0,0.0);
//	std::array<int, 8> a_AGP_Index = {{-1,-1,-1,-1,-1,-1,-1,-1}};
//	std::array<double, 8> a_AGP_ShapeValue = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
//	std::array<glm::dvec3, 8> a_AGP_ShapeGradient = {{glm::dvec3(0.0,0.0,0.0),glm::dvec3(0.0,0.0,0.0),glm::dvec3(0.0,0.0,0.0),glm::dvec3(0.0,0.0,0.0),glm::dvec3(0.0,0.0,0.0),glm::dvec3(0.0,0.0,0.0),glm::dvec3(0.0,0.0,0.0),glm::dvec3(0.0,0.0,0.0)}};
};

class MaterialPoint_CPDI_CC : public MaterialPoint_BC
{
	public:
//		glm::dvec3	d3_Corner[4] = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 0.0)};

		std::array<Corner_ST, 4> a_Corner;

		std::vector<AGPstruct> v_AGP;

//		std::array<std::array<int, 8>, 4> a_Corner_AGP_Index = {{ {{0,0,0,0,0,0,0,0}}, {{0,0,0,0,0,0,0,0}}, {{0,0,0,0,0,0,0,0}}, {{0,0,0,0,0,0,0,0}} }};
//		std::array<std::array<int, 8>, 4> a_Corner_AGP_ShapeValue = {{ {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}, {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}, {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}, {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}} }};
//		std::array<std::array<glm::dvec3, 8>, 4> a_Corner_AGP_ShapeValue = {{ {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}, {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}, {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}, {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}} }};
	protected:
	private:
};

#endif
