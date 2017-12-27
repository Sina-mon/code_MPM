#ifndef MATERIAL_BC_H
#define MATERIAL_BC_H

#include "Definitions.h"

#define	__ELASTIC 			0
#define	__NEOHOOKEAN		1
#define	__VISCOELASTIC		2
#define	__PLASTIC			3
#define	__VONMISESHARDENING	4
#define	__GASS				5

class Material_BC
{
	public:
		Material_BC() {;}
		virtual ~Material_BC() {;}

		unsigned int i_ID = 0;
		unsigned int i_MaterialType = __ELASTIC;

		double d_Density = 0.0;

		double d_ElasticModulus = 0.0;
		double d_PoissonRatio = 0.0;

		double d_YieldStress = 0.0;

		double d_Viscosity = 0.0;

		double d_Hardening_Isotropic_C0 = 0.0;
		double d_Hardening_Isotropic_C1 = 0.0;
	protected:
	private:
};

#endif
