#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation_CPDI_MultiBody_SinglePass_MPLocks(double dTimeIncrement_Total)
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
				thisGP->b_Contact = false;
				thisGP->d_Mass = 0.0;
				thisGP->d3_MassGradient = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
				thisGP->d3_Force_Temp = glm::dvec3(0.0, 0.0, 0.0);

				for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
				{
					GridPoint *thisGP = allGridPoint_Body[index_Body][index_GP];

					thisGP->b_Active = false;
					thisGP->b_Contact = false;
					thisGP->d_Mass = 0.0;
					thisGP->d3_MassGradient = glm::dvec3(0.0, 0.0, 0.0);
					thisGP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
					thisGP->d3_Force = glm::dvec3(0.0, 0.0, 0.0);
					thisGP->d3_Force_Temp = glm::dvec3(0.0, 0.0, 0.0);
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
					GridPoint *thisAGP;
//					if(thisMP->i_Body < _MAX_N_BODIES)
						thisAGP = allGridPoint_Body[thisMP->i_Body][thisMP->v_AGP[index_AGP].index];
//					else
//						thisAGP = allGridPoint_Body[0][thisMP->v_AGP[index_AGP].index];

//					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

					// shape value and shape gradient value
					double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

					// mass
					if(iThread_Count > 1)	omp_set_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
					{
						thisAGP->d_Mass += dShapeValue * thisMP->d_Mass;
						thisAGP->d3_MassGradient += d3ShapeGradient * thisMP->d_Mass;
					}
					if(iThread_Count > 1)	omp_unset_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
				}
			}
			a_Runtime[2] += omp_get_wtime() - dRuntime_Block;
			#pragma omp barrier
			// accumulate GP-Body values ----------------------------------- GP-Body
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];
				if(thisGP->b_Active == true)
				{
					for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
					{
						GridPoint *thisGP_Body = allGridPoint_Body[index_Body][index_GP];

						thisGP->d_Mass += thisGP_Body->d_Mass;
						thisGP->d3_MassGradient += thisGP_Body->d3_MassGradient;
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
					GridPoint *thisAGP;
//					if(thisMP->i_Body < _MAX_N_BODIES)
						thisAGP = allGridPoint_Body[thisMP->i_Body][thisMP->v_AGP[index_AGP].index];
//					else
//						thisAGP = allGridPoint_Body[0][thisMP->v_AGP[index_AGP].index];

//					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

					// shape value and shape gradient value
					double dShapeValue = thisMP->v_AGP[index_AGP].dShapeValue;
					glm::dvec3 d3ShapeGradient = thisMP->v_AGP[index_AGP].d3ShapeGradient;

					if(iThread_Count > 1)	omp_set_lock(v_GridPoint_Lock[thisMP->v_AGP[index_AGP].index]);
					{
						// velocity
						if(thisAGP->d_Mass > d_Mass_Minimum)
							thisAGP->d3_Velocity += dShapeValue * (thisMP->d_Mass * thisMP->d3_Velocity) / thisAGP->d_Mass;
//						if(allGridPoint[thisMP->v_AGP[index_AGP].index]->d_Mass > d_Mass_Minimum)
//							thisAGP->d3_Velocity += dShapeValue * (thisMP->d_Mass * thisMP->d3_Velocity) / allGridPoint[thisMP->v_AGP[index_AGP].index]->d_Mass;

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
			a_Runtime[3] += omp_get_wtime() - dRuntime_Block;
			#pragma omp barrier
			// accumulate GP-Body values ----------------------------------- GP-Body
			#pragma omp for
			for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
			{
				GridPoint *thisGP = allGridPoint[index_GP];
				if(thisGP->b_Active == true)
				{
					// add all bodies
					for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
					{
						GridPoint *thisGP_Body = allGridPoint_Body[index_Body][index_GP];

//						if(thisGP_Body->d_Mass > d_Mass_Minimum)
							thisGP->d3_Velocity	+= (thisGP_Body->d3_Velocity * thisGP_Body->d_Mass) / thisGP->d_Mass;
//						thisGP->d3_Velocity	+= thisGP_Body->d3_Velocity;
						//thisGP->d3_Momentum	+= thisGP_Body->d3_Momentum;
						thisGP->d3_Force	+= thisGP_Body->d3_Force;
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

				// helper combined variable
				GridPoint combinedGP;
				{
					combinedGP.d3_Velocity	= glm::dvec3(0.0,0.0,0.0);
					combinedGP.d3_Force		= glm::dvec3(0.0,0.0,0.0);
				}
				for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
				{
					GridPoint *thisGP_Body = allGridPoint_Body[index_Body][index_GP];

//					if(thisGP->d_Mass > d_Mass_Minimum)
						combinedGP.d3_Velocity	+= (thisGP_Body->d_Mass * thisGP_Body->d3_Velocity)/thisGP->d_Mass;
					combinedGP.d3_Force	+= thisGP_Body->d3_Force;
				}
				// check for contact
				for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
				{
					GridPoint *thisGP_Body = allGridPoint_Body[index_Body][index_GP];

//					glm::dvec3 d3Normal = glm::dvec3(0.0,0.0,0.0);
//					if(glm::length(thisGP_Body->d3_MassGradient) > d_Mass_Minimum)
					glm::dvec3 d3Normal = glm::normalize(thisGP_Body->d3_MassGradient);
					double dContact = glm::dot(thisGP_Body->d3_Velocity - combinedGP.d3_Velocity, d3Normal);
					if(dContact > 1.0e-12)
					{// if there is contact, adjust the normal velocity component
						thisGP->b_Contact = true;

						thisGP_Body->d3_Velocity = thisGP_Body->d3_Velocity - dContact*d3Normal;

						thisGP_Body->d3_Force = (-thisGP_Body->d_Mass * dContact / dTimeIncrement) * d3Normal;
					}
				}
				// update body momenta and apply boundary conditions
				for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
				{
					GridPoint *thisGP_Body = allGridPoint_Body[index_Body][index_GP];

					if(thisGP_Body->d_Mass > d_Mass_Minimum)
						thisGP_Body->d3_Velocity += thisGP_Body->d3_Force / thisGP_Body->d_Mass * dTimeIncrement;

					if(thisGP->b3_Fixed.x == true)
					{
						thisGP_Body->d3_Velocity.x = 0.0;
						thisGP_Body->d3_Force.x = 0.0;
					}
					if(thisGP->b3_Fixed.y == true)
					{
						thisGP_Body->d3_Velocity.y = 0.0;
						thisGP->d3_Force_Temp.y += thisGP_Body->d3_Force.y;
						thisGP_Body->d3_Force.y = 0.0;
					}
					if(thisGP->b3_Fixed.z == true)
					{
						thisGP_Body->d3_Velocity.z = 0.0;
						thisGP_Body->d3_Force.z = 0.0;
					}
				}


//				if(thisGP->d_Mass > d_Mass_Minimum)
//					thisGP->d3_Velocity += thisGP->d3_Force / thisGP->d_Mass * dTimeIncrement;
//
//				if(thisGP->b3_Fixed.x == true)
//				{
//					thisGP->d3_Velocity.x = 0.0;
//					thisGP->d3_Force.x = 0.0;
//				}
//				if(thisGP->b3_Fixed.y == true)
//				{
//					thisGP->d3_Velocity.y = 0.0;
//					thisGP->d3_Force_Temp.y += thisGP->d3_Force.y;
//					thisGP->d3_Force.y = 0.0;
//				}
//				if(thisGP->b3_Fixed.z == true)
//				{
//					thisGP->d3_Velocity.z = 0.0;
//					thisGP->d3_Force.z = 0.0;
//				}

			}
			a_Runtime[4] += omp_get_wtime() - dRuntime_Block;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// displacement controlled material points ------------------------ displacement control
			#pragma omp for
			for(unsigned int index_MP = 0; index_MP < v_MarkedMaterialPoints_CPDI_Displacement_Control.size(); index_MP++)
			{
				MaterialPoint_CPDI_CC *thisMP = v_MarkedMaterialPoints_CPDI_Displacement_Control[index_MP];

				for(unsigned int index_AGP = 0; index_AGP < thisMP->v_AGP.size(); index_AGP++)
				{
					GridPoint *thisAGP;
						thisAGP = allGridPoint_Body[thisMP->i_Body][thisMP->v_AGP[index_AGP].index];

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
					GridPoint *thisAGP;
						thisAGP = allGridPoint_Body[thisMP->i_Body][thisMP->v_AGP[index_AGP].index];

//					GridPoint *thisAGP = allGridPoint[thisMP->v_AGP[index_AGP].index];

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

				for(int index = 0; index < 6; index++)
					d6StressIncrement[index] = vonMises_Thread.d6StressIncrement[index];

				for(int index = 0; index < 6; index++)
					d6PlasticStrainIncrement[index] = vonMises_Thread.d6PlasticStrainIncrement[index];

				for(int index = 0; index < 6; index++)
					thisMP->d6_Stress[index] += d6StressIncrement[index];

				for(int index = 0; index < 6; index++)
					thisMP->d6_Strain_Plastic[index] += d6PlasticStrainIncrement[index];

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
						GridPoint *thisAGP;
							thisAGP = allGridPoint_Body[thisMP->i_Body][mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints[index_AGP]];

//						GridPoint *thisAGP = allGridPoint[mpm_GP_Mediator_Thread[iThread_This].v_adjacentGridPoints[index_AGP]];

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
			// sina, need to upadate this from classic to tetrahedra corners
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
