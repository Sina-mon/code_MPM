/*
make a convenience function to set all outpur variables (stressincrement, plastic strain increment, etc to zero)
*/

#ifndef CONSTITUTIVERELATION_H
#define CONSTITUTIVERELATION_H

#include <iostream>

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp" // linux

#include "Material_BC.h"

class ConstitutiveRelation
{
	public:
		ConstitutiveRelation(){;}
		virtual ~ConstitutiveRelation(){;}

		double d6StressIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6PlasticStrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double dBackstress_IsotropicIncrement = 0.0;
		double d_I1 = 0.0;
		double d_I2 = 0.0;
		double d_I3 = 0.0;
		double d_J2 = 0.0;
		double d_J3 = 0.0;

		void calculateIncrement_Elastic(Material_BC *pMaterial, double d6StrainIncrement[6]);
		void calculateIncrement_Plastic(Material_BC *pMaterial, double d6StressCurrent[6], double d6StrainIncrement[6]);
		void calculateIncrement_RambergOsgood(Material_BC *pMaterial, double dBackstress_Isotropic, double d6StressCurrent[6], double d6StrainIncrement[6]);
		void calculateIncrement_VonMisesHardening(Material_BC *pMaterial, double dBackstress_Isotropic, double d6StressCurrent[6], double d6StrainIncrement[6]);
		void calculateIncrement_VMHS(Material_BC *pMaterial, double dBackstress_Isotropic, double d6StressCurrent[6], double d6StrainCurrent[6], double d6StrainIncrement[6]);

		void calculateIncrement_Elastic(double dE, double dNu, double d6StrainIncrement[6]);
		void calculateIncrement_NeoHookean_6D(double dE, double dNu, double d6Stress_Current[6], glm::dmat3 d33DeformationGradient);
		void calculateIncrement_ViscoElastic_6D(double dE, double dNu, double dEta, double d6Stress_Current[6], double d6Strain_New[6], double d6Strain_Rate_New[6]);
		void calculateIncrement_PerfectlyPlastic_6D(double dE, double dNu, double dYield, double d6StressCurrent[6], double d6StrainIncrement[6]);
		void calculateIncrement_VonMisesHardening_6D(double dE, double dNu, double dYield, double dBackstress_Isotropic, double dHardening_Isotropic_C0, double dHardening_Isotropic_C1, double d6StressCurrent[6], double d6StrainIncrement[6]);
		void calculateIncrement_IdealGass(double dHeatCapacityRatio, double dSpecificHeat, double dShearViscosity, double dBulkViscosity, double d6Stress_Current[6], double dDensity, double dTemperature, double d6Strain_Rate[6]);

		double getState_I1(double d6State[6]);
		double getState_I2(double d6State[6]);
		double getState_I3(double d6State[6]);
		double getState_J2(double d6State[6]);
		double getState_J3(double d6State[6]);
		glm::dvec3 getPrincipal(double d6State[6]);

		void calculateState_I1(double d6State[6]);
		void calculateState_I2(double d6State[6]);
		void calculateState_I3(double d6State[6]);
		void calculateState_J2(double d6State[6]);
		void calculateState_J3(double d6State[6]);

		glm::dmat3 getExp(glm::dmat3 d33A);// https://en.wikipedia.org/wiki/Matrix_exponential
		glm::dmat3 getLog(glm::dmat3 d33A);// https://en.wikipedia.org/wiki/Logarithm_of_a_matrix
	protected:

	private:
};

#endif // CONSTITUTIVERELATION_H
