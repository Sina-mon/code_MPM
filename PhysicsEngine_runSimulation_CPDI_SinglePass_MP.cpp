#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_CPDI_SinglePass_MP(double dTimeIncrement_Total)
{
	omp_set_num_threads(_MAX_N_THREADS);

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	double dRuntime_MP = 0.0;
	double dRuntime_Block = 0.0;
	double dDebug_ContactCutoff = -1000.0;

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
			dRuntime_Block = omp_get_wtime();
			// reset grid points ---------------------------------------------- reset grid points
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == false)
					continue;

				thisGP->b_Active = false;
				thisGP->d_Mass = 0.0;
				thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force_Temp = glm::dvec3(0.0, 0.0, 0.0);

				for(int index_Thread = 0; index_Thread < iThread_Count; index_Thread++)
				{
					GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

					thisGP_Thread->b_Active = false;
					thisGP_Thread->d_Mass = 0.0;
					thisGP_Thread->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
					thisGP_Thread->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
					thisGP_Thread->d3_Force_Temp = glm::dvec3(0.0, 0.0, 0.0);
				}
			}
			a_Runtime[0] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// Find AGP's and calculate shape values and gradients
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint_CPDI.size(); index_MP++)
			{
				MaterialPoint_CPDI_CC *thisMP = allMaterialPoint_CPDI[index_MP];
				thisMP->v_AGP.clear();

				mpm_GP_Mediator_Thread[iThread_This].findAdjacentGridPoints_CPDI(thisMP);
				for(unsigned int index_AGP = 0; index_AGP < mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints[index_AGP]];
					thisAGP->b_Active = true;
					mpm_GP_Mediator_Thread[iThread_This].calculateBases_CPDI(thisMP, thisAGP->d3_Position);

					// shape value and shape gradient value
					AGPstruct thisAGPstruct;
					thisAGPstruct.index = mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints[index_AGP];
					thisAGPstruct.dShapeValue = mpm_GP_Mediator_Thread[iThread_This].d_ShapeValue;
					thisAGPstruct.d3ShapeGradient = mpm_GP_Mediator_Thread[iThread_This].d3_ShapeGradient;
					thisMP->v_AGP.push_back(thisAGPstruct);
				}
			}
			a_Runtime[1] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// material point to grid, mass only
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint_CPDI.size(); index_MP++)
			{
				MaterialPoint_CPDI_CC *thisMP = allMaterialPoint_CPDI[index_MP];

				if(thisMP->b_DisplacementControl == true)
					continue;

				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP_Thread = allGridPoint_Thread[iThread_This][thisMP->v_AGP[index_AGP].index];

					// shape value and shape gradient value
					double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

					// mass
					thisAGP_Thread->d_Mass += dShapeValue * thisMP->d_Mass;
				}
			}
			#pragma omp barrier
			// accumulate GP-layered values ----------------------------------- GP-layered
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];
				if(thisGP->b_Active == true)
				{
					for(int index_Thread = 0; index_Thread < iThread_Count; index_Thread++)
					{
						GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

						thisGP->d_Mass += thisGP_Thread->d_Mass;
					}
				}
			}
			a_Runtime[2] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// material point to grid, velocity and force
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint_CPDI.size(); index_MP++)
			{
				MaterialPoint_CPDI_CC *thisMP = allMaterialPoint_CPDI[index_MP];

				if(thisMP->b_DisplacementControl == true)
					continue;

				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];
					GridPoint *thisAGP_Thread = allGridPoint_Thread[iThread_This][thisMP->v_AGP[index_AGP].index];

					// shape value and shape gradient value
					double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

					// velocity
					if(thisAGP->d_Mass > d_Mass_Minimum)
//						thisAGP->d3_Velocity += dShapeValue * (thisMP->d_Mass * thisMP->d3_Velocity) / thisAGP->d_Mass;
						thisAGP_Thread->d3_Velocity += dShapeValue * (thisMP->d_Mass * thisMP->d3_Velocity) / thisAGP->d_Mass;

					// internal forces
					double dVolume = thisMP->d_Volume;
//					thisAGP->d3_Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]);
//					thisAGP->d3_Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]);
//					thisAGP->d3_Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]);
					thisAGP_Thread->d3_Force.x += -dVolume * (d3ShapeGradient.x*thisMP->d6_Stress[0] + d3ShapeGradient.y*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[5]);
					thisAGP_Thread->d3_Force.y += -dVolume * (d3ShapeGradient.y*thisMP->d6_Stress[1] + d3ShapeGradient.x*thisMP->d6_Stress[3] + d3ShapeGradient.z*thisMP->d6_Stress[4]);
					thisAGP_Thread->d3_Force.z += -dVolume * (d3ShapeGradient.z*thisMP->d6_Stress[2] + d3ShapeGradient.x*thisMP->d6_Stress[5] + d3ShapeGradient.y*thisMP->d6_Stress[4]);

					// external forces
//					thisAGP->d3_Force += dShapeValue*thisMP->d3_Force_External;
					thisAGP_Thread->d3_Force += dShapeValue*thisMP->d3_Force_External;
				}
			}
			#pragma omp barrier
			// accumulate GP-layered values ----------------------------------- GP-layered
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];
				if(thisGP->b_Active == true)
				{
					for(int index_Thread = 0; index_Thread < iThread_Count; index_Thread++)
					{
						GridPoint *thisGP_Thread = allGridPoint_Thread[index_Thread][index_GP];

						thisGP->d3_Velocity	+= thisGP_Thread->d3_Velocity;
						thisGP->d3_Force	+= thisGP_Thread->d3_Force;
					}
				}
			}
			a_Runtime[3] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// update grid momentum and apply boundary conditions ------------- update GP momentum and damping
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == false)
					continue;

				if(glm::length(thisGP->d3_Velocity) > 1.0e-9)
					thisGP->d3_Force -= d_DampingCoefficient * glm::length(thisGP->d3_Force) * glm::normalize(thisGP->d3_Velocity);

				if(thisGP->d_Mass > d_Mass_Minimum)
					thisGP->d3_Velocity += thisGP->d3_Force / thisGP->d_Mass * dTimeIncrement;

				if(thisGP->b3_Fixed.x == true)
				{
					thisGP->d3_Velocity.x = 0.0;
					thisGP->d3_Force.x = 0.0;
				}
				if(thisGP->b3_Fixed.y == true)
				{
					thisGP->d3_Velocity.y = 0.0;
					thisGP->d3_Force_Temp.y += thisGP->d3_Force.y;
					thisGP->d3_Force.y = 0.0;
				}
				if(thisGP->b3_Fixed.z == true)
				{
					thisGP->d3_Velocity.z = 0.0;
					thisGP->d3_Force.z = 0.0;
				}
			}
			a_Runtime[4] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_CPDI_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint_CPDI_CC *thisMP = v_MarkedMaterialPoints_CPDI_Displacement_Control[index_MP];

				mpm_GP_Mediator_Thread[iThread_This].findAdjacentGridPoints(thisMP->d3_Position);

				for(int index_AGP = 0; index_AGP < mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints.size(); index_AGP++)
				{
					unsigned int index_GP = mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints[index_AGP];
					GridPoint *thisAGP = allGridPoint[index_GP];
//				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
//				{
//					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

					thisAGP->d3_Velocity = thisMP->d3_Velocity;
					//thisAGP->d3_Force_Temp += thisAGP->d3_Force;
					thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
				}
			}
			a_Runtime[5] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// grid to material ----------------------------------------------- GP to MP
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint_CPDI.size(); index_MP++)
			{
				MaterialPoint_CPDI_CC *thisMP = allMaterialPoint_CPDI[index_MP];

				if(thisMP->b_DisplacementControl == true)
					continue;

				glm::dmat3 d33VelocityGradient = glm::dmat3(0.0);

				//thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

					// shape value and shape gradient value
					double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

					if(thisAGP->d_Mass > d_Mass_Minimum)
						thisMP->d3_Velocity += dShapeValue * (thisAGP->d3_Force/thisAGP->d_Mass) * dTimeIncrement; // and before the loop, remove thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
						//thisMP->d3_Velocity += dShapeValue * thisAGP->d3_Velocity; // and before the loop, add thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);

					// velocity gradient, to be used to calculate strains
					d33VelocityGradient += glm::outerProduct(thisAGP->d3_Velocity, d3ShapeGradient);// this glm function does the pre-transposition that we want
				}

				thisMP->d33_DeformationGradient += (d33VelocityGradient * thisMP->d33_DeformationGradient) * dTimeIncrement;

				glm::dmat3 d33DeformationGradientIncrement = glm::dmat3(1.0) + d33VelocityGradient * dTimeIncrement;

				double d6StrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				d6StrainIncrement[0] = d33DeformationGradientIncrement[0][0] - 1.0;
				d6StrainIncrement[1] = d33DeformationGradientIncrement[1][1] - 1.0;
				d6StrainIncrement[2] = d33DeformationGradientIncrement[2][2] - 1.0;
				d6StrainIncrement[3] = d33DeformationGradientIncrement[0][1] + d33DeformationGradientIncrement[1][0];
				d6StrainIncrement[4] = d33DeformationGradientIncrement[1][2] + d33DeformationGradientIncrement[2][1];
				d6StrainIncrement[5] = d33DeformationGradientIncrement[2][0] + d33DeformationGradientIncrement[0][2];

				double d6StrainRate[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				for(int index = 0; index < 6; index++)
					d6StrainRate[index] = d6StrainIncrement[index] / dTimeIncrement;

				for(int index = 0; index < 6; index++)
					thisMP->d6_Strain[index] += d6StrainIncrement[index];

				// elastic
				double dE = thisMP->d_ElasticModulus;
				// plastic
				double dNu = thisMP->d_PoissonRatio;
				double dYield = thisMP->d_YieldStress;

				double d6StressIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double d6PlasticStrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				ConstitutiveRelation vonMises_Thread;

				if(thisMP->i_MaterialType == _ELASTIC)
					vonMises_Thread.calculateIncrement_Elastic(dE, dNu, d6StrainIncrement);
				else if(thisMP->i_MaterialType == _PLASTIC)
					vonMises_Thread.calculateIncrement_PerfectlyPlastic_6D(dE, dNu, dYield, thisMP->d6_Stress, d6StrainIncrement);
				else if(thisMP->i_MaterialType == _VONMISESHARDENING)
					vonMises_Thread.calculateIncrement_VonMisesHardening_6D(dE, dNu, dYield, thisMP->d_BackStress_Isotropic, thisMP->d_Hardening_Isotropic_C0, thisMP->d_Hardening_Isotropic_C1, thisMP->d6_Stress, d6StrainIncrement);
				else
					vonMises_Thread.calculateIncrement_PerfectlyPlastic_6D(dE, dNu, dYield, thisMP->d6_Stress, d6StrainIncrement);

				// update MP variables
				for(int index = 0; index < 6; index++)
					d6StressIncrement[index] = vonMises_Thread.d6StressIncrement[index];
				for(int index = 0; index < 6; index++)
					d6PlasticStrainIncrement[index] = vonMises_Thread.d6PlasticStrainIncrement[index];
				for(int index = 0; index < 6; index++)
					thisMP->d6_Stress[index] += d6StressIncrement[index];
				for(int index = 0; index < 6; index++)
					thisMP->d6_Strain_Plastic[index] += d6PlasticStrainIncrement[index];

				for(int index = 0; index < 6; index++)
					thisMP->d_Energy_Strain += thisMP->d6_Stress[index]*d6StrainIncrement[index] * thisMP->d_Volume;
				for(int index = 0; index < 6; index++)
					thisMP->d_Energy_Plastic += thisMP->d6_Stress[index]*d6PlasticStrainIncrement[index] * thisMP->d_Volume;

				thisMP->d_BackStress_Isotropic += vonMises_Thread.dBackstress_IsotropicIncrement;
			}
			a_Runtime[6] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			// update corner positions --------------------------------------------
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint_CPDI.size(); index_MP++)
			{
				MaterialPoint_CPDI_CC *thisMP = allMaterialPoint_CPDI[index_MP];

				if(thisMP->b_DisplacementControl == true)
					continue;

				for(unsigned int index_Corner = 0; index_Corner < 4; index_Corner++)
				{
					glm::dvec3 d3Displacement = glm::dvec3(0.0, 0.0, 0.0);

					mpm_GP_Mediator_Thread[iThread_This].findAdjacentGridPoints(thisMP->a_Corner[index_Corner].d3_Position);
					for(unsigned int index_AGP = 0; index_AGP < mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints.size(); index_AGP++)
					{
						GridPoint *thisAGP = allGridPoint[mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints[index_AGP]];

						mpm_GP_Mediator_Thread[iThread_This].calculateBases_Classic(thisMP->a_Corner[index_Corner].d3_Position, thisAGP->d3_Position);

						double dShapeValue = mpm_GP_Mediator_Thread[iThread_This].d_ShapeValue;
						glm::dvec3 d3ShapeGradient = mpm_GP_Mediator_Thread[iThread_This].d3_ShapeGradient;

						d3Displacement += dShapeValue * (thisAGP->d3_Velocity) * dTimeIncrement;
					}

					thisMP->a_Corner[index_Corner].d3_Position += d3Displacement;
				}
				// calculate material point position, as average of corners
				thisMP->d3_Position = glm::dvec3(0.0, 0.0, 0.0);
				for(unsigned int index_Corner = 0; index_Corner < 4; index_Corner++)
					thisMP->d3_Position += 0.25 * thisMP->a_Corner[index_Corner].d3_Position;
				// calculate material point volume
				MaterialPoint_Factory_CPDI_CC MP_Factory;
				thisMP->d_Volume = MP_Factory.getVolume(thisMP);
			}

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_CPDI_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint_CPDI_CC *thisMP = v_MarkedMaterialPoints_CPDI_Displacement_Control[index_MP];

				if(m_TimeLine.v_Time.size() != 0)
					thisMP->d3_Velocity = glm::dvec3(thisMP->f_DisplacementControl_Multiplier) * m_TimeLine.getVelocity(d_Time);

				thisMP->d3_Position += thisMP->d3_Velocity * dTimeIncrement;
				for(unsigned int index_Corner = 0; index_Corner < 4; index_Corner++)
					thisMP->a_Corner[index_Corner].d3_Position += thisMP->d3_Velocity * dTimeIncrement;
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
