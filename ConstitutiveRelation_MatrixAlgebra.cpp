#include "ConstitutiveRelation.h"

double getFactorial(int i)
{
	double dResult = 1.0;
	for(int iCount = 1; iCount <= i; iCount++)
	{
		dResult *= iCount;
	}
	return(dResult);
}
glm::dmat3 getPower(glm::dmat3 d33A, int iPow)
{
	glm::dmat3 d33Result = glm::dmat3(1.0);

	for(int iCount = 1; iCount <= iPow; iCount++)
	{
		d33Result *= d33A;
	}

	return(d33Result);

}
// ----------------------------------------------------------------------------
glm::dmat3 ConstitutiveRelation::getExp(glm::dmat3 d33A)
{// https://en.wikipedia.org/wiki/Logarithm_of_a_matrix
	glm::dmat3 d33Result = glm::dmat3(0.0);

	int iTerm_Max = 8;

	for(int iTerm = 0; iTerm <= iTerm_Max; iTerm++)
	{
		d33Result += getPower(d33A, iTerm) / getFactorial(iTerm);
	}

	return(d33Result);
}
// ----------------------------------------------------------------------------
glm::dmat3 ConstitutiveRelation::getLog(glm::dmat3 d33A)
{// https://en.wikipedia.org/wiki/Logarithm_of_a_matrix
	glm::dmat3 d33Result = glm::dmat3(0.0);

	int iTerm_Max = 8;

	glm::dmat3 d33AminusI = d33A - glm::dmat3(1.0);

	for(int iTerm = 1; iTerm <= iTerm_Max; iTerm++)
	{
		d33Result += getPower(d33AminusI, iTerm) * (1.0/iTerm) * glm::pow(-1.0, iTerm+1);
	}

	return(d33Result);
}
