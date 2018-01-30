#include "PhysicsEngine.h"

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
				thisGP->d3_MassGradient = glm::dvec3(0.0,0.0,0.0);
				thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Momentum = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force_Temp = glm::dvec3(0.0, 0.0, 0.0);
			}
			a_Runtime[0] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// Find AGP's and calculate shape values and gradients
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
			{
				MaterialPoint_BC *thisMP = allMaterialPoint[index_MP];
				thisMP->v_AGP.clear();

				mpm_GP_Mediator_Thread[iThread_This].findAdjacentGridPoints(thisMP->d3_Position);

				for(int index_AGP = 0; index_AGP < mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints[index_AGP]];
					thisAGP->b_Active = true;
					mpm_GP_Mediator_Thread[iThread_This].calculateBases_Classic(thisMP->d3_Position, thisAGP->d3_Position);

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
			// material point to grid: mass, momentum and force
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

					if(iThread_Count > 1)	omp_set_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
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
					if(iThread_Count > 1)	omp_unset_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
				}
			}

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint_BC *thisMP = v_MarkedMaterialPoints_Displacement_Control[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

//					thisAGP->d3_Velocity = thisMP->d3_Velocity;
//					thisAGP->d3_Momentum = thisAGP->d_Mass * thisMP->d3_Velocity;
//					thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);

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
			a_Runtime[5] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// update grid momentum and apply boundary conditions ------------- update GP momentum
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];

				if(thisGP->b_Active == false)
					continue;

				if(glm::length(thisGP->d3_Momentum) > 1.0e-24)
					thisGP->d3_Force -= d_DampingCoefficient * glm::length(glm::dot(thisGP->d3_Force,glm::normalize(thisGP->d3_Momentum))) * glm::normalize(thisGP->d3_Momentum);

//				if(thisGP->d_Mass > d_Mass_Minimum)
//					thisGP->d3_Velocity += thisGP->d3_Force / thisGP->d_Mass * dTimeIncrement;

				thisGP->d3_Momentum += thisGP->d3_Force * dTimeIncrement;

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
					thisGP->d3_Force_Temp.y += thisGP->d3_Force.y;
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
			dRuntime_Block = omp_get_wtime();
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint_BC *thisMP = v_MarkedMaterialPoints_Displacement_Control[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

//					thisAGP->d3_Velocity = thisMP->d3_Velocity;
//					thisAGP->d3_Momentum = thisAGP->d_Mass * thisMP->d3_Velocity;
//					thisAGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);

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
			a_Runtime[5] += omp_get_wtime() - dRuntime_Block;

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
//					thisGP->d3_Force_Temp.y += thisGP->d3_Force.y;
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

				thisMP->d33_DeformationGradient += (d33VelocityGradient * thisMP->d33_DeformationGradient) * dTimeIncrement;

				double dDet = glm::determinant(thisMP->d33_DeformationGradient);
				thisMP->d_Volume = dDet * thisMP->d_Volume_Initial;

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

				double d6StressIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				double d6PlasticStrainIncrement[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				ConstitutiveRelation CR_Thread;

				if(thisMP->p_Material == NULL)
				{
					std::cout << "Error in PhysicsEngine::runSimulation_CPDI_MultiBody_SinglePass_MPLocks, material not specified" << std::endl;
					continue;
				}

				if(thisMP->p_Material->i_MaterialType == __ELASTIC)
					CR_Thread.calculateIncrement_Elastic(thisMP->p_Material, d6StrainIncrement);
					//CR_Thread.calculateIncrement_Elastic(dE, dNu, d6StrainIncrement);
				else if(thisMP->p_Material->i_MaterialType == __PLASTIC)
					CR_Thread.calculateIncrement_Plastic(thisMP->p_Material, thisMP->d6_Stress, d6StrainIncrement);
				else if(thisMP->p_Material->i_MaterialType == __VONMISESHARDENING)
					CR_Thread.calculateIncrement_VonMisesHardening(thisMP->p_Material, thisMP->d_BackStress_Isotropic, thisMP->d6_Stress, d6StrainIncrement);
				else
					CR_Thread.calculateIncrement_Plastic(thisMP->p_Material, thisMP->d6_Stress, d6StrainIncrement);

				// update MP variables
				for(int index = 0; index < 6; index++)
					thisMP->d6_Strain_Rate[index] = d6StrainIncrement[index] / dTimeIncrement;
				for(int index = 0; index < 6; index++)
					d6StressIncrement[index] = CR_Thread.d6StressIncrement[index];
				for(int index = 0; index < 6; index++)
					d6PlasticStrainIncrement[index] = CR_Thread.d6PlasticStrainIncrement[index];
				for(int index = 0; index < 6; index++)
					thisMP->d6_Stress[index] += d6StressIncrement[index];
				for(int index = 0; index < 6; index++)
					thisMP->d6_Strain_Plastic[index] += d6PlasticStrainIncrement[index];
				// energies
				for(int index = 0; index < 6; index++)
					thisMP->d_Energy_Strain += thisMP->d6_Stress[index]*d6StrainIncrement[index] * thisMP->d_Volume;
				for(int index = 0; index < 6; index++)
					thisMP->d_Energy_Plastic += thisMP->d6_Stress[index]*d6PlasticStrainIncrement[index] * thisMP->d_Volume;

				thisMP->d_BackStress_Isotropic += CR_Thread.dBackstress_IsotropicIncrement;
			}
			a_Runtime[6] += omp_get_wtime() - dRuntime_Block;

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
