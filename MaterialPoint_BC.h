#ifndef MATERIALPOINT_BC_H
#define MATERIALPOINT_BC_H

#include "Definitions.h"

enum MaterialType
{
	_ELASTIC,
	_NEOHOOKEAN,
	_VISCOELASTIC,
	_PLASTIC,
	_VONMISESHARDENING,
	_GASS,
};

struct AGPstruct
{
	unsigned int index = 0;
	double dShapeValue = 0.0;
	glm::dvec3 d3ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
};

class MaterialPoint_BC
{
	public:
		MaterialPoint_BC() {;}
		virtual ~MaterialPoint_BC() {;}

		bool b_DisplacementControl = false;
		float f_DisplacementControl_Multiplier = 1.0;
		bool b_Mark_Force = false;
		bool b_Mark_Stress = false;
		bool b_Mark_Energy_Strain = false;
		bool b_Mark_Energy_Plastic = false;
		bool b_Surface = false;

		unsigned int i_ID = 0;
		unsigned int i_Body = 0;
		unsigned int i_MaterialType = _ELASTIC;

		double d_Mass = 0.0;
		double d_Volume = 0.0;
		double d_Volume_Initial = 0.0;

		double d_ElasticModulus = 0.0;
		double d_PoissonRatio = 0.0;
		double d_YieldStress = 0.0;
		double d_Viscosity = 0.0;

		glm::dvec3 d3_Position = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Force_External = glm::dvec3(0.0, 0.0, 0.0);

		glm::dmat3 d33_DeformationGradient = glm::dmat3(1.0);

		double d6_Strain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6_Stress[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6_Strain_Plastic[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		double d_Energy_Strain = 0.0;
		double d_Energy_Plastic = 0.0;

		double d_BackStress_Isotropic = 0.0;
		double d_Hardening_Isotropic_C0 = 0.0;
		double d_Hardening_Isotropic_C1 = 0.0;

		std::vector<AGPstruct> v_AGP;
	protected:
	private:
};

#endif
