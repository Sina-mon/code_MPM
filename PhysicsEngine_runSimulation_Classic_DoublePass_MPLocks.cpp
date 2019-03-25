#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_ResetGrid(void)
{
	#pragma omp for
	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{
		GridPoint *thisGP = allGridPoint[index_GP];

		if(thisGP->b_Active == false)
			continue;

		thisGP->b_Active = false;
		thisGP->d_Mass = 0.0;
		//thisGP->d3_MassGradient = glm::dvec3(0.0,0.0,0.0);
		thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
		thisGP->d3_Momentum = glm::dvec3(0.0, 0.0, 0.0);
		thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
		thisGP->d3_Force_Temp = glm::dvec3(0.0, 0.0, 0.0);
	}
}
int PhysicsEngine::runSimulation_FindAGPs(int iThread)
{
	#pragma omp for
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{
		MaterialPoint_BC *thisMP = allMaterialPoint[index_MP];
		thisMP->v_AGP.clear();

		mpm_GP_Mediator_Thread[iThread].findAdjacentGridPoints(thisMP->d3_Position);

		for(int index_AGP = 0; index_AGP < mpm_GP_Mediator_Thread[iThread].v_adjacentGridPoints.size(); index_AGP++)
		{
			GridPoint *thisAGP = allGridPoint[mpm_GP_Mediator_Thread[iThread].v_adjacentGridPoints[index_AGP]];
			thisAGP->b_Active = true;
			mpm_GP_Mediator_Thread[iThread].calculateBases_Classic(thisMP->d3_Position, thisAGP->d3_Position);

			// shape value and shape gradient value
			AGPstruct thisAGPstruct;
			thisAGPstruct.index = mpm_GP_Mediator_Thread[iThread].v_adjacentGridPoints[index_AGP];
			thisAGPstruct.dShapeValue = mpm_GP_Mediator_Thread[iThread].d_ShapeValue;
			thisAGPstruct.d3ShapeGradient = mpm_GP_Mediator_Thread[iThread].d3_ShapeGradient;
			thisMP->v_AGP.push_back(thisAGPstruct);
		}
	}
}
int PhysicsEngine::runSimulation_M2G(int nThread)
{
	#pragma omp for
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{
		MaterialPoint_BC *thisMP = allMaterialPoint[index_MP];

//				if(thisMP->b_DisplacementControl == true)
//					continue;

		for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
		{
			GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

			// shape value and shape gradient value
			double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
			glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

			if(nThread > 1)	omp_set_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
			{
				// mass
				thisAGP->d_Mass += dShapeValue * thisMP->d_Mass;

				// momentum
					thisAGP->d3_Momentum += dShapeValue * (thisMP->d_Mass * thisMP->d3_Velocity);
				// velocity, cannot be accumulated at the same time as mass
//						if(thisAGP->d_Mass > d_Mass_Minimum)
//							thisAGP->d3_Velocity += dShapeValue * (thisMP->d_Mass * thisMP->d3_Velocity) / thisAGP->d_Mass;

				// internal forces
				double dVolume = thisMP->d_Volume;
				thisAGP->d3_Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]);
				thisAGP->d3_Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]);
				thisAGP->d3_Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]);

				// external forces
				thisAGP->d3_Force += dShapeValue*thisMP->d3_Force_External;
			}
			if(nThread > 1)	omp_unset_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
		}
	}
}
int PhysicsEngine::runSimulation_IntegrateGrid(double dTimeIncrement)
{
	#pragma omp for
	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{
		GridPoint *thisGP = allGridPoint[index_GP];

		if(thisGP->b_Active == false)
			continue;

		if(glm::length(thisGP->d3_Momentum) > 1.0e-24)
			thisGP->d3_Force -= d_DampingCoefficient * glm::length(glm::dot(thisGP->d3_Force,glm::normalize(thisGP->d3_Momentum))) * glm::normalize(thisGP->d3_Momentum);

//		if(thisGP->d_Mass > d_Mass_Minimum)
//			thisGP->d3_Velocity += thisGP->d3_Force / thisGP->d_Mass * dTimeIncrement;

		if(thisGP->b3_Fixed.x == true)
		{
			glm::dvec3 d3Momentum_Trial = thisGP->d3_Momentum + thisGP->d3_Force * dTimeIncrement;

			if(glm::abs(d3Momentum_Trial.y*dTimeIncrement) > 1.0e-24 && thisGP->d3_Force.x > 0.0)
			{// friction forces
				glm::dvec3 d3Force_Friction = glm::dvec3(0.0,0.0,0.0);
				d3Force_Friction.y = thisGP->d_FrictionCoefficient * glm::sign(-d3Momentum_Trial.y) * glm::abs(thisGP->d3_Force.x);

				if(glm::abs(d3Force_Friction.y * dTimeIncrement) > glm::abs(d3Momentum_Trial.y))
					d3Force_Friction.y *= glm::abs(d3Momentum_Trial.y / dTimeIncrement) / glm::abs(d3Force_Friction.y);

				thisGP->d3_Force.y += d3Force_Friction.y;
			}

			thisGP->d3_Velocity.x = 0.0;
			thisGP->d3_Momentum.x = 0.0;
			thisGP->d3_Force.x = 0.0;
		}
		if(thisGP->b3_Fixed.y == true)
		{
			thisGP->d3_Velocity.y = 0.0;
			thisGP->d3_Momentum.y = 0.0;
			thisGP->d3_Force_Temp.y += thisGP->d3_Force.y;
			thisGP->d3_Force.y = 0.0;
		}
		if(thisGP->b3_Fixed.z == true)
		{
			thisGP->d3_Velocity.z = 0.0;
			thisGP->d3_Momentum.z = 0.0;
			thisGP->d3_Force.z = 0.0;
		}

		thisGP->d3_Momentum += thisGP->d3_Force * dTimeIncrement;

	}
}
int PhysicsEngine::runSimulation_DisplacementControl(void)
{
	#pragma omp for
	for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement_Control.size(); index_MP++)
	{
		MaterialPoint_BC *thisMP = v_MarkedMaterialPoints_Displacement_Control[index_MP];

		for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
		{
			GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

			if(thisMP->b3_DisplacementControl.x == true)
			{
				thisAGP->d3_Velocity.x = thisMP->d3_Velocity.x;
				thisAGP->d3_Momentum.x = thisAGP->d_Mass * thisMP->d3_Velocity.x;
				thisAGP->d3_Force.x = 0.0;
			}
			if(thisMP->b3_DisplacementControl.y == true)
			{
				thisAGP->d3_Velocity.y = thisMP->d3_Velocity.y;
				thisAGP->d3_Momentum.y = thisAGP->d_Mass * thisMP->d3_Velocity.y;
				thisAGP->d3_Force.y = 0.0;
			}
			if(thisMP->b3_DisplacementControl.z == true)
			{
				thisAGP->d3_Velocity.z = thisMP->d3_Velocity.z;
				thisAGP->d3_Momentum.z = thisAGP->d_Mass * thisMP->d3_Velocity.z;
				thisAGP->d3_Force.z = 0.0;
			}
		}
	}
}
int PhysicsEngine::runSimulation_G2P_P2_SmallStrain (double dTimeIncrement)
{
	// grid to material, pass 2
	#pragma omp for
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{
		MaterialPoint_BC *thisMP = allMaterialPoint[index_MP];

//				if(thisMP->b_DisplacementControl == true)
//					continue;

		glm::dmat3 d33VelocityGradient = glm::dmat3(0.0);

		for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
		{
			GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

			// shape value and shape gradient value
			double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
			glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

			// position
			if(thisAGP->d_Mass > d_Mass_Minimum)
				thisMP->d3_Position += dShapeValue * (thisAGP->d3_Momentum/thisAGP->d_Mass) * dTimeIncrement;

			// velocity gradient, to be used to calculate strains
			d33VelocityGradient += glm::outerProduct(thisAGP->d3_Velocity, d3ShapeGradient);// this glm function does the pre-transposition that we want
		}

		ConstitutiveRelation CR_Thread;

		glm::dmat3 d33F_Delta = glm::dmat3(1) + d33VelocityGradient * dTimeIncrement;

		glm::dmat3 d33Strain_Increment = glm::dmat3(0.0);
		{// infinitesimal strain increment
			d33Strain_Increment = 0.5*(glm::transpose(d33F_Delta) + d33F_Delta) - glm::dmat3(1);
		}
		double d6Strain_Increment[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		{
			d6Strain_Increment[0] = d33Strain_Increment[0][0];
			d6Strain_Increment[1] = d33Strain_Increment[1][1];
			d6Strain_Increment[2] = d33Strain_Increment[2][2];
			d6Strain_Increment[3] = d33Strain_Increment[0][1] + d33Strain_Increment[1][0];
			d6Strain_Increment[4] = d33Strain_Increment[1][2] + d33Strain_Increment[2][1];
			d6Strain_Increment[5] = d33Strain_Increment[2][0] + d33Strain_Increment[0][2];
		}

		double d6StrainRate[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		for(int index = 0; index < 6; index++)
			d6StrainRate[index] = d6Strain_Increment[index] / dTimeIncrement;

		if(thisMP->p_Material == NULL)
		{
			std::cout << "Error in PhysicsEngine::runSimulation_CPDI_MultiBody_SinglePass_MPLocks, material not specified" << std::endl;
			continue;
		}

		if(thisMP->p_Material->i_MaterialType == __ELASTIC)
			CR_Thread.calculateIncrement_Elastic(thisMP->p_Material, d6Strain_Increment);
		else if(thisMP->p_Material->i_MaterialType == __PLASTIC)
			CR_Thread.calculateIncrement_Plastic(thisMP->p_Material, thisMP->d6_Stress, d6Strain_Increment);
		else if(thisMP->p_Material->i_MaterialType == __VONMISESHARDENING)
			CR_Thread.calculateIncrement_VonMisesHardening(thisMP->p_Material, thisMP->d_BackStress_Isotropic, thisMP->d6_Stress, d6Strain_Increment);
		else
			CR_Thread.calculateIncrement_Plastic(thisMP->p_Material, thisMP->d6_Stress, d6Strain_Increment);

		// total deformation gradient
		thisMP->d33_DeformationGradient *= d33F_Delta;
		// volume
		double dDet = glm::determinant(thisMP->d33_DeformationGradient);
		thisMP->d_Volume = dDet * thisMP->d_Volume_Initial;
		// strain
		for(int index = 0; index < 6; index++)
			thisMP->d6_Strain[index] += d6Strain_Increment[index];
		// strain rate
//		for(int index = 0; index < 6; index++)
//			thisMP->d6_Strain_Rate[index] = d6Strain_Increment[index] / dTimeIncrement;
		// stress
		for(int index = 0; index < 6; index++)
			thisMP->d6_Stress[index] += CR_Thread.d6StressIncrement[index];
		// plastic strain
		for(int index = 0; index < 6; index++)
			thisMP->d6_Strain_Plastic[index] += CR_Thread.d6PlasticStrainIncrement[index];
		// total strain energy
		double dStrainEnergy_Increment = 0.0;
		for(int index = 0; index < 6; index++)
			dStrainEnergy_Increment += (thisMP->d6_Stress[index]-0.5*CR_Thread.d6StressIncrement[index]) * d6Strain_Increment[index] * thisMP->d_Volume; // sina, be careful, the "-0.5*CR_Thread.d6StressIncrement[index]" is important for getting correct results
		thisMP->d_Energy_Strain += dStrainEnergy_Increment;

		thisMP->d_Energy_Strain_TimeIntegral += thisMP->d_Energy_Strain * dTimeIncrement;

		if(thisMP->d_Energy_Strain > thisMP->d_Energy_Strain_MaxHistory)
			thisMP->d_Energy_Strain_MaxHistory = thisMP->d_Energy_Strain;
//		double strainEnergy = 0.0;
//		for(int index = 0; index < 6; index++)
//			strainEnergy += 0.5*thisMP->d6_Stress[index] * thisMP->d6_Strain[index] * thisMP->d_Volume;
//		thisMP->d_Energy_Strain = strainEnergy;
		// plastic strain energy
		double dPlasticEnergy_Increment = 0.0;
		for(int index = 0; index < 6; index++)
			dPlasticEnergy_Increment += (thisMP->d6_Stress[index]-0.5*CR_Thread.d6StressIncrement[index]) * CR_Thread.d6PlasticStrainIncrement[index] * thisMP->d_Volume;
		thisMP->d_Energy_Plastic += dPlasticEnergy_Increment;
		// hardening variables
		thisMP->d_BackStress_Isotropic += CR_Thread.dBackstress_IsotropicIncrement;
	}
}
int PhysicsEngine::runSimulation_G2P_P2_LargeStrain (double dTimeIncrement)
{
	// grid to material, pass 2
	#pragma omp for
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{
		MaterialPoint_BC *thisMP = allMaterialPoint[index_MP];

//				if(thisMP->b_DisplacementControl == true)
//					continue;

		glm::dmat3 d33VelocityGradient = glm::dmat3(0.0);

		for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
		{
			GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

			// shape value and shape gradient value
			double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
			glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

			// position
			if(thisAGP->d_Mass > d_Mass_Minimum)
				thisMP->d3_Position += dShapeValue * (thisAGP->d3_Momentum/thisAGP->d_Mass) * dTimeIncrement;

			// velocity gradient, to be used to calculate strains
			d33VelocityGradient += glm::outerProduct(thisAGP->d3_Velocity, d3ShapeGradient);// this glm function does the pre-transposition that we want
		}

		ConstitutiveRelation CR_Thread;

		// reset stresses for large strain formulation
		for(int index = 0; index < 6; index++)
			thisMP->d6_Stress[index] = 0.0;

		glm::dmat3 d33F_Delta = glm::dmat3(1) + d33VelocityGradient * dTimeIncrement;
		glm::dmat3 d33F_Total = d33F_Delta * thisMP->d33_DeformationGradient;

		// total strain
		double d6Strain_Total[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		glm::dmat3 d33Strain_Total = glm::dmat3(0.0);
		{
			d33Strain_Total = 0.5*CR_Thread.getLog(thisMP->d33_DeformationGradient * glm::transpose(thisMP->d33_DeformationGradient));

			d6Strain_Total[0] = d33Strain_Total[0][0];
			d6Strain_Total[1] = d33Strain_Total[1][1];
			d6Strain_Total[2] = d33Strain_Total[2][2];
			d6Strain_Total[3] = d33Strain_Total[0][1] + d33Strain_Total[1][0];
			d6Strain_Total[4] = d33Strain_Total[1][2] + d33Strain_Total[2][1];
			d6Strain_Total[5] = d33Strain_Total[2][0] + d33Strain_Total[0][2];
		}

		glm::dmat3 d33Strain_Increment = glm::dmat3(0.0);
		glm::dmat3 d33Strain_Elastic_Trial = glm::dmat3(0.0);
		{// large-strain plasticity Chrisfield Box 12.2
			// initialize
			glm::dmat3 d33F_0_Elastic = thisMP->d33_DeformationGradient_Elastic;
			// trial elastic
			glm::dmat3 d33F_Elastic_Trial = d33F_Delta * d33F_0_Elastic;

			// logarithmic strain tensor trial
			d33Strain_Elastic_Trial = 0.5*CR_Thread.getLog(d33F_Elastic_Trial * glm::transpose(d33F_Elastic_Trial));
			d33Strain_Increment		= 0.5*CR_Thread.getLog(d33F_Delta * glm::transpose(d33F_Delta));

			// green lagrange strain tensor trial
			//d33Strain_Elastic_Trial = 0.5*(glm::transpose(d33F_Elastic_Trial) * d33F_Elastic_Trial - glm::dmat3(1));
			//d33Strain_Increment = 0.5*(glm::transpose(d33F_Delta) * d33F_Delta - glm::dmat3(1));

			// infinitesimal strain increment
			//d33Strain_Elastic_Trial = 0.5*(glm::transpose(d33F_Elastic_Trial) + d33F_Elastic_Trial) - glm::dmat3(1);
			//d33Strain_Increment = 0.5*(glm::transpose(d33F_Delta) + d33F_Delta) - glm::dmat3(1);
		}
		double d6Strain_Increment[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double d6Strain_Elastic_Trial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		{
			d6Strain_Elastic_Trial[0] = d33Strain_Elastic_Trial[0][0];
			d6Strain_Elastic_Trial[1] = d33Strain_Elastic_Trial[1][1];
			d6Strain_Elastic_Trial[2] = d33Strain_Elastic_Trial[2][2];
			d6Strain_Elastic_Trial[3] = d33Strain_Elastic_Trial[0][1] + d33Strain_Elastic_Trial[1][0];
			d6Strain_Elastic_Trial[4] = d33Strain_Elastic_Trial[1][2] + d33Strain_Elastic_Trial[2][1];
			d6Strain_Elastic_Trial[5] = d33Strain_Elastic_Trial[2][0] + d33Strain_Elastic_Trial[0][2];

			d6Strain_Increment[0] = d33Strain_Increment[0][0];
			d6Strain_Increment[1] = d33Strain_Increment[1][1];
			d6Strain_Increment[2] = d33Strain_Increment[2][2];
			d6Strain_Increment[3] = d33Strain_Increment[0][1] + d33Strain_Increment[1][0];
			d6Strain_Increment[4] = d33Strain_Increment[1][2] + d33Strain_Increment[2][1];
			d6Strain_Increment[5] = d33Strain_Increment[2][0] + d33Strain_Increment[0][2];
		}

		double d6StrainRate[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		for(int index = 0; index < 6; index++)
			d6StrainRate[index] = d6Strain_Increment[index] / dTimeIncrement;

		if(thisMP->p_Material == NULL)
		{
			std::cout << "Error in PhysicsEngine::runSimulation_CPDI_MultiBody_SinglePass_MPLocks, material not specified" << std::endl;
			continue;
		}

		if(thisMP->p_Material->i_MaterialType == __ELASTIC)
			CR_Thread.calculateIncrement_Elastic(thisMP->p_Material, d6Strain_Elastic_Trial);
		else if(thisMP->p_Material->i_MaterialType == __PLASTIC)
			CR_Thread.calculateIncrement_Plastic(thisMP->p_Material, thisMP->d6_Stress, d6Strain_Elastic_Trial);
		else if(thisMP->p_Material->i_MaterialType == __VONMISESHARDENING)
			CR_Thread.calculateIncrement_VonMisesHardening(thisMP->p_Material, thisMP->d_BackStress_Isotropic, thisMP->d6_Stress, d6Strain_Elastic_Trial);
		else
			CR_Thread.calculateIncrement_Plastic(thisMP->p_Material, thisMP->d6_Stress, d6Strain_Elastic_Trial);

		// update MP variables
		// elastic deformation gradient
		double d6Strain_Elastic[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		for(int index = 0; index < 6; index++)
			d6Strain_Elastic[index] = d6Strain_Elastic_Trial[index] - CR_Thread.d6PlasticStrainIncrement[index];
		glm::dmat3 d33Strain_Elastic = glm::dmat3(0.0);
		{
			d33Strain_Elastic[0][0] = d6Strain_Elastic[0];
			d33Strain_Elastic[1][1] = d6Strain_Elastic[1];
			d33Strain_Elastic[2][2] = d6Strain_Elastic[2];

			d33Strain_Elastic[0][1] = 0.5*d6Strain_Elastic[3];
			d33Strain_Elastic[1][0] = 0.5*d6Strain_Elastic[3];

			d33Strain_Elastic[1][2] = 0.5*d6Strain_Elastic[4];
			d33Strain_Elastic[2][1] = 0.5*d6Strain_Elastic[4];

			d33Strain_Elastic[2][0] = 0.5*d6Strain_Elastic[5];
			d33Strain_Elastic[0][2] = 0.5*d6Strain_Elastic[5];
		}
		thisMP->d33_DeformationGradient_Elastic = CR_Thread.getExp(d33Strain_Elastic);
		// total deformation gradient
		thisMP->d33_DeformationGradient *= d33F_Delta;
		// volume
		double dDet_F = glm::determinant(thisMP->d33_DeformationGradient);
		thisMP->d_Volume = dDet_F * thisMP->d_Volume_Initial;
		// strain
		for(int index = 0; index < 6; index++)
			thisMP->d6_Strain[index] = d6Strain_Total[index];//+= d6StrainIncrement[index];
		// strain rate
//		for(int index = 0; index < 6; index++)
//			thisMP->d6_Strain_Rate[index] = d6Strain_Increment[index] / dTimeIncrement;
		// stress
		for(int index = 0; index < 6; index++)
			thisMP->d6_Stress[index] = (1.0/dDet_F) * CR_Thread.d6StressIncrement[index];
		// plastic strain
		for(int index = 0; index < 6; index++)
			thisMP->d6_Strain_Plastic[index] += CR_Thread.d6PlasticStrainIncrement[index];
		// total strain energy
		for(int index = 0; index < 6; index++)
			thisMP->d_Energy_Strain += thisMP->d6_Stress[index]*d6Strain_Increment[index] * thisMP->d_Volume;
		// plastic strain energy
		for(int index = 0; index < 6; index++)
			thisMP->d_Energy_Plastic += thisMP->d6_Stress[index]*CR_Thread.d6PlasticStrainIncrement[index] * thisMP->d_Volume;
		// hardening variables
		thisMP->d_BackStress_Isotropic += CR_Thread.dBackstress_IsotropicIncrement;
	}
}
// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_Classic_DoublePass_MPLocks(double dTimeIncrement_Total)
{
	omp_set_num_threads(_MAX_N_THREADS);

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	double dRuntime_MP = 0.0;
	double dRuntime_Block = 0.0;

	double dTimeIncrement_Accumulated = 0.0;
	while(dTimeIncrement_Accumulated < dTimeIncrement_Total)
	{
		// default time-increment value
		double dTimeIncrement = d_TimeIncrement_Maximum;
		// unless reaching the end, where the remaining increment is less than the maximum
		if(d_TimeIncrement_Maximum > (dTimeIncrement_Total-dTimeIncrement_Accumulated))
			dTimeIncrement = dTimeIncrement_Total - dTimeIncrement_Accumulated;
		// calculate increments accumulated
		dTimeIncrement_Accumulated += dTimeIncrement;

		clockCurrent_Total = clock();
		dRuntime_MP = omp_get_wtime();
		#pragma omp parallel
		{
			int iThread_Count = omp_get_num_threads();
			int	iThread_This = omp_get_thread_num();

			#pragma omp barrier
			// reset grid points ---------------------------------------------- reset grid points
			this->runSimulation_ResetGrid();
			#pragma omp barrier
			// Find AGP's and calculate shape values and gradients
			this->runSimulation_FindAGPs(iThread_This);

			#pragma omp barrier
			// material point to grid: mass, momentum and force
			this->runSimulation_M2G(iThread_Count);

			#pragma omp barrier
			// displacement controlled material points ------------------------ displacement control
			this->runSimulation_DisplacementControl();

			#pragma omp barrier
			// update grid momentum and apply boundary conditions ------------- update GP momentum
			this->runSimulation_IntegrateGrid(dTimeIncrement);

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// grid to material, pass 1
			// only update particle velocities
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint_BC *thisMP = allMaterialPoint[index_MP];

//				if(thisMP->b_DisplacementControl == true)
//					continue;

				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

					// shape value and shape gradient value
					double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

					// velocity
					if(thisAGP->d_Mass > d_Mass_Minimum)
						thisMP->d3_Velocity += dShapeValue * (thisAGP->d3_Force/thisAGP->d_Mass) * dTimeIncrement;
				}
			}
			a_Runtime[6] += omp_get_wtime() - dRuntime_Block;

/*			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// reset grid velocites ---------------------------------------------- reset grid points
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == true)
				{
					thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);

					for(int index_Thread = 0; index_Thread < iThread_Count; index_Thread++)
					{
						GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

						thisGP_Thread->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
					}
				}
			}
			a_Runtime[0] += omp_get_wtime() - dRuntime_Block;
*/
			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// map particle velocity back to grid -------------------------
			// mass in NOT mapped here --------------------------------------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint_BC *thisMP = allMaterialPoint[index_MP];

//				if(thisMP->b_DisplacementControl == true)
//					continue;

				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

					// shape value and shape gradient value
					double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

					// velocity
					if(iThread_Count > 1)	omp_set_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
					{
						// velocity
						if(thisAGP->d_Mass > d_Mass_Minimum)
							thisAGP->d3_Velocity += dShapeValue * (thisMP->d_Mass * thisMP->d3_Velocity) / thisAGP->d_Mass;
					}
					if(iThread_Count > 1)	omp_unset_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
				}
			}

			#pragma omp barrier
			// displacement controlled material points ------------------------ displacement control
			this->runSimulation_DisplacementControl();

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// re-apply boundary conditions -----------------------------------
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == false)
					continue;

				if(thisGP->b3_Fixed.x == true)
				{
					thisGP->d3_Velocity.x = 0.0;
					thisGP->d3_Momentum.x = 0.0;
					thisGP->d3_Force.x = 0.0;
				}
				if(thisGP->b3_Fixed.y == true)
				{
					thisGP->d3_Velocity.y = 0.0;
					thisGP->d3_Momentum.y = 0.0;
					thisGP->d3_Force.y = 0.0;
				}
				if(thisGP->b3_Fixed.z == true)
				{
					thisGP->d3_Velocity.z = 0.0;
					thisGP->d3_Momentum.z = 0.0;
					thisGP->d3_Force.z = 0.0;
				}
			}
			a_Runtime[4] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			// grid to particle: pass 2
			this->runSimulation_G2P_P2_SmallStrain(dTimeIncrement);
			//this->runSimulation_G2P_P2_LargeStrain(dTimeIncrement);
/*
			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// grid to material, pass 2
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint_BC *thisMP = allMaterialPoint[index_MP];

//				if(thisMP->b_DisplacementControl == true)
//					continue;

				glm::dmat3 d33VelocityGradient = glm::dmat3(0.0);

				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

					// shape value and shape gradient value
					double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

					// position
					if(thisAGP->d_Mass > d_Mass_Minimum)
						thisMP->d3_Position += dShapeValue * (thisAGP->d3_Momentum/thisAGP->d_Mass) * dTimeIncrement;

					// velocity gradient, to be used to calculate strains
					d33VelocityGradient += glm::outerProduct(thisAGP->d3_Velocity, d3ShapeGradient);// this glm function does the pre-transposition that we want
				}

				ConstitutiveRelation CR_Thread;

				// reset stresses for large strain formulation
				for(int index = 0; index < 6; index++)
					thisMP->d6_Stress[index] = 0.0;

				thisMP->d33_DeformationGradient += (d33VelocityGradient * thisMP->d33_DeformationGradient) * dTimeIncrement;
				double dDet = glm::determinant(thisMP->d33_DeformationGradient);
				thisMP->d_Volume = dDet * thisMP->d_Volume_Initial;

				glm::dmat3 d33DeformationGradientIncrement = glm::dmat3(1.0) + d33VelocityGradient * dTimeIncrement;
				// total strain
				double d6Strain_Total[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				glm::dmat3 d33Strain_Total = glm::dmat3(0.0);
				{
					d33Strain_Total = 0.5*CR_Thread.getLog(thisMP->d33_DeformationGradient * glm::transpose(thisMP->d33_DeformationGradient));

					d6Strain_Total[0] = d33Strain_Total[0][0];
					d6Strain_Total[1] = d33Strain_Total[1][1];
					d6Strain_Total[2] = d33Strain_Total[2][2];
					d6Strain_Total[3] = d33Strain_Total[0][1] + d33Strain_Total[1][0];
					d6Strain_Total[4] = d33Strain_Total[1][2] + d33Strain_Total[2][1];
					d6Strain_Total[5] = d33Strain_Total[2][0] + d33Strain_Total[0][2];
				}

				glm::dmat3 d33Strain_Increment = glm::dmat3(0.0);
				glm::dmat3 d33Strain_Elastic_Trial = glm::dmat3(0.0);
				{// large-strain plasticity
					// initialize
					glm::dmat3 d33F_0_Total = thisMP->d33_DeformationGradient;
					glm::dmat3 d33F_0_Elastic = thisMP->d33_DeformationGradient_Elastic;
					// deformation gradient increment and new total value
					glm::dmat3 d33F_Delta = glm::dmat3(1) + d33VelocityGradient * dTimeIncrement;
					glm::dmat3 d33F_Total = d33F_Delta * d33F_0_Total;
					// trial elastic
					glm::dmat3 d33F_Elastic_Trial = d33F_Delta * d33F_0_Elastic;

					// logarithmic strain tensor trial
					d33Strain_Elastic_Trial = 0.5*CR_Thread.getLog(d33F_Elastic_Trial * glm::transpose(d33F_Elastic_Trial));
					d33Strain_Increment		= 0.5*CR_Thread.getLog(d33F_Delta * glm::transpose(d33F_Delta));

					// green lagrange strain tensor trial
					//d33Strain_Elastic_Trial = 0.5*(glm::transpose(d33F_Elastic_Trial) * d33F_Elastic_Trial - glm::dmat3(1));
					//d33Strain_Increment = 0.5*(glm::transpose(d33F_Delta) * d33F_Delta - glm::dmat3(1));

					// infinitesimal strain increment
					//d33Strain_Elastic_Trial = 0.5*(glm::transpose(d33F_Elastic_Trial) + d33F_Elastic_Trial) - glm::dmat3(1);
					//d33Strain_Increment = 0.5*(glm::transpose(d33F_Delta) + d33F_Delta) - glm::dmat3(1);
				}
				double d6Strain_Increment[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double d6Strain_Elastic_Trial[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				{
					d6Strain_Elastic_Trial[0] = d33Strain_Elastic_Trial[0][0];
					d6Strain_Elastic_Trial[1] = d33Strain_Elastic_Trial[1][1];
					d6Strain_Elastic_Trial[2] = d33Strain_Elastic_Trial[2][2];
					d6Strain_Elastic_Trial[3] = d33Strain_Elastic_Trial[0][1] + d33Strain_Elastic_Trial[1][0];
					d6Strain_Elastic_Trial[4] = d33Strain_Elastic_Trial[1][2] + d33Strain_Elastic_Trial[2][1];
					d6Strain_Elastic_Trial[5] = d33Strain_Elastic_Trial[2][0] + d33Strain_Elastic_Trial[0][2];

					d6Strain_Increment[0] = d33Strain_Increment[0][0];
					d6Strain_Increment[1] = d33Strain_Increment[1][1];
					d6Strain_Increment[2] = d33Strain_Increment[2][2];
					d6Strain_Increment[3] = d33Strain_Increment[0][1] + d33Strain_Increment[1][0];
					d6Strain_Increment[4] = d33Strain_Increment[1][2] + d33Strain_Increment[2][1];
					d6Strain_Increment[5] = d33Strain_Increment[2][0] + d33Strain_Increment[0][2];
				}

				double d6StrainRate[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				for(int index = 0; index < 6; index++)
					d6StrainRate[index] = d6Strain_Increment[index] / dTimeIncrement;

				if(thisMP->p_Material == NULL)
				{
					std::cout << "Error in PhysicsEngine::runSimulation_CPDI_MultiBody_SinglePass_MPLocks, material not specified" << std::endl;
					continue;
				}

				if(thisMP->p_Material->i_MaterialType == __ELASTIC)
					CR_Thread.calculateIncrement_Elastic(thisMP->p_Material, d6Strain_Elastic_Trial);
				else if(thisMP->p_Material->i_MaterialType == __PLASTIC)
					CR_Thread.calculateIncrement_Plastic(thisMP->p_Material, thisMP->d6_Stress, d6Strain_Elastic_Trial);
				else if(thisMP->p_Material->i_MaterialType == __VONMISESHARDENING)
					CR_Thread.calculateIncrement_VonMisesHardening(thisMP->p_Material, thisMP->d_BackStress_Isotropic, thisMP->d6_Stress, d6Strain_Elastic_Trial);
				else
					CR_Thread.calculateIncrement_Plastic(thisMP->p_Material, thisMP->d6_Stress, d6Strain_Elastic_Trial);

				// large strain plasticity
				// correct elastic strain
				double d6Strain_Elastic[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				for(int index = 0; index < 6; index++)// elastic strain
					d6Strain_Elastic[index] = d6Strain_Elastic_Trial[index] - CR_Thread.d6PlasticStrainIncrement[index];
				glm::dmat3 d33Strain_Elastic = glm::dmat3(0.0);
				{
					d33Strain_Elastic[0][0] = d6Strain_Elastic[0];
					d33Strain_Elastic[1][1] = d6Strain_Elastic[1];
					d33Strain_Elastic[2][2] = d6Strain_Elastic[2];

					d33Strain_Elastic[0][1] = 0.5*d6Strain_Elastic[3];
					d33Strain_Elastic[1][0] = 0.5*d6Strain_Elastic[3];

					d33Strain_Elastic[1][2] = 0.5*d6Strain_Elastic[4];
					d33Strain_Elastic[2][1] = 0.5*d6Strain_Elastic[4];

					d33Strain_Elastic[2][0] = 0.5*d6Strain_Elastic[5];
					d33Strain_Elastic[0][2] = 0.5*d6Strain_Elastic[5];
				}
				thisMP->d33_DeformationGradient_Elastic = CR_Thread.getExp(d33Strain_Elastic);

				// update MP variables
				for(int index = 0; index < 6; index++)// strain
					thisMP->d6_Strain[index] = d6Strain_Total[index];//+= d6StrainIncrement[index];
//				for(int index = 0; index < 6; index++)// strain rate
//					thisMP->d6_Strain_Rate[index] = d6StrainIncrement[index] / dTimeIncrement;
				for(int index = 0; index < 6; index++)// stress
					thisMP->d6_Stress[index] += CR_Thread.d6StressIncrement[index];
				for(int index = 0; index < 6; index++)// plastic strain
					thisMP->d6_Strain_Plastic[index] += CR_Thread.d6PlasticStrainIncrement[index];
				// energies
				for(int index = 0; index < 6; index++)
					thisMP->d_Energy_Strain += thisMP->d6_Stress[index]*d6Strain_Increment[index] * thisMP->d_Volume;
				for(int index = 0; index < 6; index++)
					thisMP->d_Energy_Plastic += thisMP->d6_Stress[index]*CR_Thread.d6PlasticStrainIncrement[index] * thisMP->d_Volume;

				thisMP->d_BackStress_Isotropic += CR_Thread.dBackstress_IsotropicIncrement;
			}
			a_Runtime[6] += omp_get_wtime() - dRuntime_Block;
*/
			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint_BC *thisMP = v_MarkedMaterialPoints_Displacement_Control[index_MP];

				if(m_TimeLine.v_Time.size() != 0)
				{
					if(thisMP->b3_DisplacementControl.x == true)
					{
						thisMP->d3_Velocity.x = thisMP->f_DisplacementControl_Multiplier * m_TimeLine.getVelocity(d_Time).x;
					}
					if(thisMP->b3_DisplacementControl.y == true)
					{
						thisMP->d3_Velocity.y = thisMP->f_DisplacementControl_Multiplier * m_TimeLine.getVelocity(d_Time).y;
					}
					if(thisMP->b3_DisplacementControl.z == true)
					{
						thisMP->d3_Velocity.z = thisMP->f_DisplacementControl_Multiplier * m_TimeLine.getVelocity(d_Time).z;
					}
				}
//					thisMP->d3_Velocity = glm::dvec3(thisMP->f_DisplacementControl_Multiplier) * m_TimeLine.getVelocity(d_Time);

				//thisMP->d3_Position += thisMP->d3_Velocity * dTimeIncrement;
			}
			a_Runtime[7] += omp_get_wtime() - dRuntime_Block;

		}

		d_Runtime_Total += omp_get_wtime() - dRuntime_MP;
		//report to console ---------------------------------------------------
		if(d_Time - d_TimeConsole_Last > d_TimeConsole_Interval)
		{
			d_TimeConsole_Last = d_Time;
			this->reportConsole();
		}

		d_Time += dTimeIncrement;
		i_TimeCycle++;
	}

	if(d_Time < d_TimeEnd)
		return(0);
	else
	{
		d_TimeConsole_Last = d_Time;
		this->reportConsole();
		std::cout << "Analysis complete ..." << std::endl;
		return(1);
	}
}
// ----------------------------------------------------------------------------
