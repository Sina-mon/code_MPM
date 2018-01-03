#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
PhysicsEngine::~PhysicsEngine()
{
	//delete all objects created by Factory classes
	for(unsigned int index = 0; index < allGridPoint.size(); index++)
        delete allGridPoint[index];

	for(unsigned int index = 0; index < allMaterialPoint.size(); index++)
        delete allMaterialPoint[index];

	for(unsigned int index = 0; index < v_allMaterial.size(); index++)
        delete v_allMaterial[index];

//	for(int iThread = 0; iThread < _MAX_N_THREADS; iThread++)
//	{
//		for(unsigned int index = 0; index < allGridPoint_Thread[iThread].size(); index++)
//		{
//			delete allGridPoint_Thread[iThread][index];
//		}
//	}

	for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
	{
		for(unsigned int index = 0; index < allGridPoint_Body[index_Body].size(); index++)
		{
			delete allGridPoint_Body[index_Body][index];
		}
	}

	// destroy GridPoint locks
	for(int index = 0; index < v_GridPoint_Lock.size(); index++)
	{
		omp_destroy_lock(v_GridPoint_Lock[index]);
		delete v_GridPoint_Lock[index];
	}
}
// ----------------------------------------------------------------------------
void PhysicsEngine::reportConsole(std::string sDescription)
{
	std::string strConsole = "";
	strConsole += sDescription;
	strConsole += "\ttime: " + Script(d_Time,4);
	strConsole += "\tRuntime: " + Script(d_Runtime_Total,3);
//	strConsole += "\tCPS: " + Script((d_Time/d_TimeIncrement_Maximum)/d_Runtime_Total,3);

	if(v_MarkedMaterialPoints_Displacement_Monitor.size() > 0)
	{ // displacement
        strConsole += "\tPosition_y: " + Script(v_MarkedMaterialPoints_Displacement_Monitor[0]->d3_Position.y,4);
	}

	{// force
		glm::dvec3 d3Force = glm::dvec3(0.0,0.0,0.0);
		for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
		{// calculate debug values
			GridPoint *thisGP = allGridPoint[index_GP];

			d3Force += thisGP->d3_Force_Temp;
		}
		strConsole += "\tForce_y: " + Script(d3Force.y,4);
	}

	if(v_MarkedMaterialPoints_Momentum.size() > 0)
	{// momentum
		glm::dvec3 d3Momentum = glm::dvec3(0.0,0.0,0.0);
		for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Momentum.size(); index_MP++)
		{// calculate debug values
			d3Momentum += v_MarkedMaterialPoints_Momentum[index_MP]->d3_Velocity * v_MarkedMaterialPoints_Momentum[index_MP]->d_Mass;
		}
		strConsole += "\t momentum_y: " + Script(d3Momentum.y,4);
	}

	if(v_MarkedMaterialPoints_Principal_Monitor.size() > 0)
	{// principal
		ConstitutiveRelation CR;

		double dStrain_max = 0.0;
		double dStrainRate_max = 0.0;
		double dStress_max = 0.0;

		glm::dvec3 d3Strain		= glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3StrainRate	= glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3Stress		= glm::dvec3(0.0,0.0,0.0);
		for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Principal_Monitor.size(); index_MP++)
		{// calculate debug values
			MaterialPoint_BC *thisMP = v_MarkedMaterialPoints_Principal_Monitor[index_MP];

			d3Strain		= CR.getPrincipal(thisMP->d6_Strain);
			d3StrainRate	= CR.getPrincipal(thisMP->d6_Strain_Rate);
			d3Stress		= CR.getPrincipal(thisMP->d6_Stress);

			d3Strain		= glm::abs(d3Strain);
			d3StrainRate	= glm::abs(d3StrainRate);
			d3Stress		= glm::abs(d3Stress);

			double dStrain_max_temp = glm::max(glm::max(d3Strain.x, d3Strain.y), d3Strain.z);
			double dStrainRate_max_temp = glm::max(glm::max(d3StrainRate.x, d3StrainRate.y), d3StrainRate.z);
			double dStress_max_temp = glm::max(glm::max(d3Stress.x, d3Stress.y), d3Stress.z);

			if(dStrain_max < dStrain_max_temp)
				dStrain_max = dStrain_max_temp;
			if(dStrainRate_max < dStrainRate_max_temp)
				dStrainRate_max = dStrainRate_max_temp;
			if(dStress_max < dStress_max_temp)
				dStress_max = dStress_max_temp;
		}

		strConsole += "\tStrain_p_max: " + Script(dStrain_max,4);
		strConsole += "\tStrainRate_p_max: " + Script(dStrainRate_max,4);
		strConsole += "\tStress_p_max: " + Script(dStress_max,4);
	}

	if(v_MarkedMaterialPoints_Stress_Monitor.size() > 0)
	{// stress-strain
		glm::dvec3 d3Stress = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3Strain = glm::dvec3(0.0,0.0,0.0);
		for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Stress_Monitor.size(); index_MP++)
		{// calculate debug values
			MaterialPoint_BC *thisMP = v_MarkedMaterialPoints_Stress_Monitor[index_MP];

			d3Stress += glm::dvec3(thisMP->d6_Stress[0], thisMP->d6_Stress[1], thisMP->d6_Stress[2]);
			d3Strain += glm::dvec3(thisMP->d6_Strain[0], thisMP->d6_Strain[1], thisMP->d6_Strain[2]);
		}
		d3Stress /= v_MarkedMaterialPoints_Stress_Monitor.size();
		d3Strain /= v_MarkedMaterialPoints_Stress_Monitor.size();

		strConsole += "\tStrain_y: " + Script(d3Strain.y,4);
		strConsole += "\tStress_y: " + Script(d3Stress.y,4);
	}

	if(v_MarkedMaterialPoints_Monitor_Energy.size() > 0)
	{// energy log
		double dKineticEnergy = 0.0;
		double dEnergy_Strain = 0.0;
		double dEnergy_Plastic = 0.0;
		for(int index = 0; index < v_MarkedMaterialPoints_Monitor_Energy.size(); index++)
		{
			MaterialPoint_BC *thisMP = v_MarkedMaterialPoints_Monitor_Energy[index];

			dKineticEnergy += 0.5*thisMP->d_Mass * glm::pow(glm::length(thisMP->d3_Velocity),2.0);
			dEnergy_Strain += thisMP->d_Energy_Strain;
			dEnergy_Plastic += thisMP->d_Energy_Plastic;
		}
		strConsole += "\tKinetic Energy: " + Script(dKineticEnergy,4);
		strConsole += "\tStrain Energy: " + Script(dEnergy_Strain,4);
		strConsole += "\tPlastic Energy: " + Script(dEnergy_Plastic,4);
	}



	strConsole += "\n";

	if(false)
	{// runtimes
		double dRuntime_Total = 0.0;
		for(int index = 0; index < a_Runtime.size(); index++)
		{
			strConsole += "\tRuntime_Block[" + Script(index) + "]: " + Script(a_Runtime[index],3) + "(" + Script(a_Runtime[index]/d_Runtime_Total,2,std::ios::fixed) + "\%)";
			strConsole += "\n";
			dRuntime_Total += a_Runtime[index];
		}
		strConsole += "\tRuntime_Total: " + Script(dRuntime_Total,3);
		strConsole += "\n";
	}

	std::cout << strConsole;

	std::string strFileName = str_Log_FileName;//;_STR_LOGFILE;
	std::ofstream OutputFile(strFileName.c_str(), std::ios_base::app);

	OutputFile << strConsole;

	OutputFile.close();
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
