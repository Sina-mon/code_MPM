#include "ConstitutiveRelation.h"

// from, https://en.wikipedia.org/wiki/Cauchy_stress_tensor#Principal_stresses_and_stress_invariants
// from, http://www.continuummechanics.org/principalstrain.html
// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateState_I1(double d6State[6])
{
	double dResult = 0.0;

	dResult = d6State[0] + d6State[1] + d6State[2];

	this->d_I1 = dResult;
}
// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateState_I2(double d6State[6])
{
	double dResult = 0.0;

	dResult += d6State[0]*d6State[1];
	dResult += d6State[1]*d6State[2];
	dResult += d6State[0]*d6State[2];
	dResult -= d6State[3]*d6State[3];
	dResult -= d6State[4]*d6State[4];
	dResult -= d6State[5]*d6State[5];

	this->d_I2 = dResult;
}
// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateState_I3(double d6State[6])
{
	glm::dmat3 d33State;
	d33State[0][0] = d6State[0];	d33State[0][1] = d6State[3];	d33State[0][2] = d6State[5];
	d33State[1][0] = d6State[3];	d33State[1][1] = d6State[1];	d33State[1][2] = d6State[4];
	d33State[2][0] = d6State[5];	d33State[2][1] = d6State[4];	d33State[2][2] = d6State[2];

	double dResult = 0.0;

//	dResult += d6State[0]*d6State[1]*d6State[2];
//	dResult += 2.0*d6State[3]*d6State[4]*d6State[5];
//	dResult -= d6State[3]*d6State[3]*d6State[2];
//	dResult -= d6State[4]*d6State[4]*d6State[0];
//	dResult -= d6State[5]*d6State[5]*d6State[1];

	dResult = glm::determinant(d33State);

	this->d_I3 = dResult;
}
// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateState_J2(double d6State[6])
{

	double dStress_Hydrostatic = 1.0/3.0 * (d6State[0] + d6State[1] + d6State[2]);

	double d6Stress_Deviatoric[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	for(int index = 0; index < 6; index++)
		d6Stress_Deviatoric[index] = d6State[index];

	d6Stress_Deviatoric[0] -= dStress_Hydrostatic;
	d6Stress_Deviatoric[1] -= dStress_Hydrostatic;
	d6Stress_Deviatoric[2] -= dStress_Hydrostatic;

	double dJ2 = 0.0;
	dJ2 += 1.0/6.0 * (d6State[0] - d6State[1])*(d6State[0] - d6State[1]);
	dJ2 += 1.0/6.0 * (d6State[1] - d6State[2])*(d6State[1] - d6State[2]);
	dJ2 += 1.0/6.0 * (d6State[2] - d6State[0])*(d6State[2] - d6State[0]);
	dJ2 += d6State[3]*d6State[3];
	dJ2 += d6State[4]*d6State[4];
	dJ2 += d6State[5]*d6State[5];

	this->d_J2 = dJ2;
}
// ----------------------------------------------------------------------------
void ConstitutiveRelation::calculateState_J3(double d6State[6])
{
	double dH = 1.0/3.0*(d6State[0] + d6State[1] + d6State[2]);

	glm::dmat3 d33Deviatoric;
	d33Deviatoric[0][0] = d6State[0] - dH;	d33Deviatoric[0][1] = d6State[3];		d33Deviatoric[0][2] = d6State[5];
	d33Deviatoric[1][0] = d6State[3];		d33Deviatoric[1][1] = d6State[1] - dH;	d33Deviatoric[1][2] = d6State[4];
	d33Deviatoric[2][0] = d6State[5];		d33Deviatoric[2][1] = d6State[4];		d33Deviatoric[2][2] = d6State[2] - dH;

	double dResult = 0.0;

	dResult = glm::determinant(d33Deviatoric);

	this->d_J3 = dResult;
}
// ----------------------------------------------------------------------------
double ConstitutiveRelation::getState_I1(double d6State[6])
{
	double dResult = 0.0;

	dResult = d6State[0] + d6State[1] + d6State[2];

	return(dResult);
}
// ----------------------------------------------------------------------------
double ConstitutiveRelation::getState_I2(double d6State[6])
{
	double dResult = 0.0;

	dResult += d6State[0]*d6State[1];
	dResult += d6State[1]*d6State[2];
	dResult += d6State[0]*d6State[2];
	dResult -= d6State[3]*d6State[3];
	dResult -= d6State[4]*d6State[4];
	dResult -= d6State[5]*d6State[5];

	return(dResult);
}
// ----------------------------------------------------------------------------
double ConstitutiveRelation::getState_I3(double d6State[6])
{
	glm::dmat3 d33State;
	d33State[0][0] = d6State[0];	d33State[0][1] = d6State[3];	d33State[0][2] = d6State[5];
	d33State[1][0] = d6State[3];	d33State[1][1] = d6State[1];	d33State[1][2] = d6State[4];
	d33State[2][0] = d6State[5];	d33State[2][1] = d6State[4];	d33State[2][2] = d6State[2];

	double dResult = 0.0;

//	dResult += d6State[0]*d6State[1]*d6State[2];
//	dResult += 2.0*d6State[3]*d6State[4]*d6State[5];
//	dResult -= d6State[3]*d6State[3]*d6State[2];
//	dResult -= d6State[4]*d6State[4]*d6State[0];
//	dResult -= d6State[5]*d6State[5]*d6State[1];

	dResult = glm::determinant(d33State);

	return(dResult);
}
// ----------------------------------------------------------------------------
double ConstitutiveRelation::getState_J2(double d6State[6])
{
	double dStress_Hydrostatic = 1.0/3.0 * (d6State[0] + d6State[1] + d6State[2]);

	double d6Stress_Deviatoric[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	for(int index = 0; index < 6; index++)
		d6Stress_Deviatoric[index] = d6State[index];

	d6Stress_Deviatoric[0] -= dStress_Hydrostatic;
	d6Stress_Deviatoric[1] -= dStress_Hydrostatic;
	d6Stress_Deviatoric[2] -= dStress_Hydrostatic;

	double dResult = 0.0;
	dResult += 1.0/6.0 * (d6State[0] - d6State[1])*(d6State[0] - d6State[1]);
	dResult += 1.0/6.0 * (d6State[1] - d6State[2])*(d6State[1] - d6State[2]);
	dResult += 1.0/6.0 * (d6State[2] - d6State[0])*(d6State[2] - d6State[0]);
	dResult += d6State[3]*d6State[3];
	dResult += d6State[4]*d6State[4];
	dResult += d6State[5]*d6State[5];

	return(dResult);
}
// ----------------------------------------------------------------------------
double ConstitutiveRelation::getState_J3(double d6State[6])
{
	double dH = 1.0/3.0*(d6State[0] + d6State[1] + d6State[2]);

	glm::dmat3 d33Deviatoric;
	d33Deviatoric[0][0] = d6State[0] - dH;	d33Deviatoric[0][1] = d6State[3];		d33Deviatoric[0][2] = d6State[5];
	d33Deviatoric[1][0] = d6State[3];		d33Deviatoric[1][1] = d6State[1] - dH;	d33Deviatoric[1][2] = d6State[4];
	d33Deviatoric[2][0] = d6State[5];		d33Deviatoric[2][1] = d6State[4];		d33Deviatoric[2][2] = d6State[2] - dH;

	double dResult = 0.0;

	dResult = glm::determinant(d33Deviatoric);

	return(dResult);
}
// ----------------------------------------------------------------------------
glm::dvec3 ConstitutiveRelation::getPrincipal(double d6State[6])
{
	// results for principal states, in no specific order
	glm::dvec3 d3Result = glm::dvec3(0.0,0.0,0.0);

	double dI1 = this->getState_I1(d6State);
	double dI2 = this->getState_I2(d6State);
	double dI3 = this->getState_I3(d6State);

	double dQ = (3.0*dI2 - dI1*dI1)/9.0;
	double dR = (2.0*dI1*dI1 - 9.0*dI1*dI2 + 27.0*dI3)/54.0;

	if(dQ > 0)
	{
		std::cout << "Arithmatic error, in ConstitutiveRelation::getPrincipal_max" << std::endl;
		return(d3Result);
	}

	double dT = dR / glm::sqrt(-dQ*dQ*dQ);
	double dTheta = glm::acos(dTheta);

	if(dTheta < 0.0 || 1.0 < dTheta)
	{
		std::cout << "Arithmatic error, in ConstitutiveRelation::getPrincipal_max" << std::endl;
		return(d3Result);
	}

	d3Result[0] = 2.0*glm::sqrt(-dQ) * glm::cos((dTheta+0.0*_PI)/3.0) + 1.0/3.0*dI1;
	d3Result[1] = 2.0*glm::sqrt(-dQ) * glm::cos((dTheta+2.0*_PI)/3.0) + 1.0/3.0*dI1;
	d3Result[2] = 2.0*glm::sqrt(-dQ) * glm::cos((dTheta+4.0*_PI)/3.0) + 1.0/3.0*dI1;

	return(d3Result);
}
// ----------------------------------------------------------------------------
