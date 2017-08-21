#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_CPDI_HalfRing_Xiang_PlainStress(void)
{
	MaterialPoint_Factory_CPDI_CC	MP_Factory;
	GridPoint_Factory				GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	glm::dvec3 d3_Length_Grid = glm::dvec3(0.040, 0.060, 0.002);
	glm::ivec3 i3_Cells = 4*glm::ivec3(40, 60, 2);
	glm::dvec3 d3_Length_Cell = d3_Length_Grid / glm::dvec3(i3_Cells);
	glm::ivec3 i3_Nodes = i3_Cells + glm::ivec3(1, 1, 1);
	for(int indexThread = 0; indexThread < _MAX_N_THREADS; indexThread++)
	{// initialize GP mediator

		d3_Length_World = d3_Length_Grid;

		mpm_GP_Mediator_Thread[indexThread].d3_Length_Grid = d3_Length_Grid;
		mpm_GP_Mediator_Thread[indexThread].d3_Length_Cell = d3_Length_Cell;
		mpm_GP_Mediator_Thread[indexThread].i3_Cells = i3_Cells;
		mpm_GP_Mediator_Thread[indexThread].i3_Node_Count = i3_Nodes;
	}
	allGridPoint = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);

	for(int iThread = 0; iThread < _MAX_N_THREADS; iThread++)
	{
		allGridPoint_Thread[iThread] = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);
	}

	// contact kernel grid ---------------------------------------------------- contact grid
	{// initialize GP kernel mediator
		d3_Length_Grid_Kernel = d3_Length_Grid;
		i3_Cells_Kernel = 1*i3_Cells;

		d3_Length_Cell_Kernel = d3_Length_Grid_Kernel / glm::dvec3(i3_Cells_Kernel);

		i3_Nodes_Kernel = i3_Cells_Kernel + glm::ivec3(1, 1, 1);
	}
	v_GridPoint_Kernel = GP_Factory.createGrid(d3_Length_Grid_Kernel, i3_Cells_Kernel);

	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{// grid point boundary conditions
		GridPoint *thisGridPoint = allGridPoint[index_GP];

		double dx = thisGridPoint->d3_Position[0];
		double dy = thisGridPoint->d3_Position[1];
		double dz = thisGridPoint->d3_Position[2];
		double dTolerance = 0.1*d3_Length_Cell.x;

		//fixed grid points
		thisGridPoint->b3_Fixed = glm::bvec3(false, false, false);

		if(fabs(dx - 0.0) < 1.5*d3_Length_Cell.x)
		{
			thisGridPoint->b3_Fixed.x = true;
		}
		if(fabs(dx - d3_Length_Grid.x) < dTolerance)
		{
		}
		if(fabs(dy - 0.0) < 0.5*d3_Length_Cell.y)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dy - d3_Length_Grid.y) < dTolerance)
		{
		}
		if(fabs(dz - 0.0) < 1.5*d3_Length_Cell.z)
		{
			thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - 0.0) < 2.0*d3_Length_Grid.z)
		{
			//thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - d3_Length_Grid.z) < dTolerance)
		{
			//thisGridPoint->b3_Fixed.z = true;
		}
	}

	// locks on grid points for atomic operations
	v_GridPoint_Lock.resize(allGridPoint.size());
	for(int index = 0; index < v_GridPoint_Lock.size(); index++)
	{
		v_GridPoint_Lock[index] = new omp_lock_t;
		omp_init_lock(v_GridPoint_Lock[index]);
	}


	double dThickness_Ring = 0.00148;// - d_Offset;
	double dDiameter_Inner = 0.0479;
	double dDiameter_Outer = dDiameter_Inner + 2.0*dThickness_Ring;

	double dAngle_Start = -_PI/2.0;
	double dAngle_End = _PI/2.0;
	double dRadius_Inner = 0.5*dDiameter_Inner;
	double dRadius_Outer = 0.5*dDiameter_Outer;
	double dLength_Ring = 0.001;//1.0*d_Offset;

	int iDivision_Angular = 360;
	int iDivision_Radial = 32;
	int iDivision_Longitudinal = 1;
	glm::dvec3 d3Center_Ring = glm::dvec3(0.0,0.0,0.5)*d3_Length_Grid;
	d3Center_Ring.x += 0.5*d3_Length_Cell.x;
	d3Center_Ring.y = 0.5*dDiameter_Outer + 3.5*d3_Length_Cell.y;
	d3Center_Ring.z = d3_Length_Cell.z + 0.5*dLength_Ring;
	if(true)
	{// ring material points -------------------------------------------------- tube MP
		double dGravity = 0.0;

//		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Tube_Smooth(d3Center_Ring, glm::dvec3(0.0,0.0,0.0), dRadius_Outer, dRadius_Inner, dLength_Ring, dLength_Ring/1.0);
		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Arc(d3Center_Ring, glm::dvec3(0.0,0.0,0.0), dAngle_Start, dAngle_End, dRadius_Outer, dRadius_Inner, dLength_Ring, iDivision_Angular, iDivision_Radial, iDivision_Longitudinal);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

			thisMP->i_MaterialType = _VONMISESHARDENING;
//			thisMP->i_MaterialType = _PLASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = MP_Factory.getVolume((MaterialPoint_CPDI_CC *)thisMP);
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 78000.0 * thisMP->d_Volume;
			d_Mass_Minimum = 0.001 * dMass;
			thisMP->d_Mass = dMass;

			thisMP->d_ElasticModulus = 210.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.3;
			thisMP->d_YieldStress = 310.0e6;

			thisMP->d_Hardening_Isotropic_C0 = 4.0;
			thisMP->d_Hardening_Isotropic_C1 = 150.0e6;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			thisMP->d3_Force_External = thisMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint_CPDI_CC *thisMP = (MaterialPoint_CPDI_CC *)thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint_CPDI.push_back(thisMP);
			// moment log
//			v_MarkedMaterialPoints_Momentum.push_back(thisMP);
			// mark for stress monitor
//			v_MarkedMaterialPoints_Stress_Monitor.push_back(thisMP);
		}
	}
	if(true)
	{// top platen material points -------------------------------------------- platen MP
		glm::dvec3 d3Dimension = glm::dvec3(0.9*d3_Length_World.x,2.0*d3_Length_Cell.y,dLength_Ring);
		glm::dvec3 d3Center = d3Center_Ring;
		d3Center.x = 0.5*d3Dimension.x;
		d3Center.y = d3Center_Ring.y + 0.5*dDiameter_Outer + 2.0*d3_Length_Cell.y + 0.5*d3Dimension.y;

		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center, d3Dimension, d3_Length_Cell.y);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

			thisMP->i_MaterialType = _ELASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = MP_Factory.getVolume((MaterialPoint_CPDI_CC *)thisMP);
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 7800.0 * thisMP->d_Volume;
			thisMP->d_Mass = dMass;

			thisMP->d_ElasticModulus = 210.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.33;
			thisMP->d_YieldStress = 200.0e6;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
//			thisMP->d3_Momentum = thisMP->d3_Mass * thisMP->d3_Velocity;
			thisMP->d3_Force_External = thisMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint_CPDI_CC *thisMP = (MaterialPoint_CPDI_CC *)thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint_CPDI.push_back(thisMP);
			// displacement control
			if(true)
			{
				thisMP->b_DisplacementControl = true;
				thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_CPDI_Displacement_Control.push_back(thisMP);
				v_MarkedMaterialPoints_CPDI_Displacement_Monitor.push_back(thisMP);
			}
		}
	}
	if(true)
	{// bottom platen material points ----------------------------------------- platen MP
		glm::dvec3 d3Dimension = glm::dvec3(0.9*d3_Length_World.x,2.0*d3_Length_Cell.y,dLength_Ring);
		glm::dvec3 d3Center = d3Center_Ring;
		d3Center.x = 0.5*d3Dimension.x;
		d3Center.y = 0.5*d3Dimension.y;

		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center, d3Dimension, d3_Length_Cell.y);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

			thisMP->i_MaterialType = _ELASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = MP_Factory.getVolume((MaterialPoint_CPDI_CC *)thisMP);
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 78000.0 * thisMP->d_Volume;
			thisMP->d_Mass = dMass;

			thisMP->d_ElasticModulus = 210.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.33;
			thisMP->d_YieldStress = 200.0e6;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			thisMP->d3_Force_External = thisMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint_CPDI_CC *thisMP = (MaterialPoint_CPDI_CC *)thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint_CPDI.push_back(thisMP);
			// stress monitor
//			v_MarkedMaterialPoints_Stress_Monitor.push_back(thisMP);
		}
	}

	double dPlatenSpeed = -0.1;
	double dTime_On  = 0.2e-3;
	double dTime_Off = 0.8e-3;
	if(true)
	{// timeline events -------------------------------------------------------
	    double dTime_Line = 0.0;

   		m_TimeLine.addTimePoint(0.0, glm::dvec3(0.0, dPlatenSpeed, 0.0));
   		m_TimeLine.addTimePoint(0.1, glm::dvec3(0.0, dPlatenSpeed, 0.0));

//   		m_TimeLine.addTimePoint(0.0, glm::dvec3(0.0, dPlatenSpeed, 0.0));
//   		m_TimeLine.addTimePoint(0.0025, glm::dvec3(0.0, dPlatenSpeed, 0.0));
//
//   		m_TimeLine.addTimePoint(0.0025+1.0e-12, glm::dvec3(0.0, 0.0, 0.0));
//   		m_TimeLine.addTimePoint(0.1, glm::dvec3(0.0, 0.0, 0.0));

//   	m_TimeLine.addTimePoint(dTime_Line, glm::dvec3(0.0, dPlatenSpeed, 0.0));
//		for(dTime_Line = 0.0; dTime_Line < 0.02; dTime_Line += (dTime_On+dTime_Off))
//		{
//			m_TimeLine.addTimePoint(dTime_Line-1.0e-12,				glm::dvec3(0.0, 0.0, 0.0));
//			m_TimeLine.addTimePoint(dTime_Line,						glm::dvec3(0.0, dPlatenSpeed, 0.0));
//
//			m_TimeLine.addTimePoint(dTime_Line+dTime_On-1.0e-12,	glm::dvec3(0.0, dPlatenSpeed, 0.0));
//			m_TimeLine.addTimePoint(dTime_Line+dTime_On,        	glm::dvec3(0.0, 0.0, 0.0));
//		}
//		m_TimeLine.addTimePoint(1.0e6,          glm::dvec3(0.0, 0.0, 0.0));
	}

	double dMass_Domain = 0.0;
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint_CPDI.size(); index_MP++)
	{// calculate debug values
		dMass_Domain += allMaterialPoint_CPDI[index_MP]->d_Mass;
	}

	a_Runtime.fill(0.0);
	d_DampingCoefficient = 0.0;

	d_TimeIncrement_Maximum = 5.0e-8;
	d_TimeEnd = 0.8*dDiameter_Outer / glm::abs(dPlatenSpeed);
	d_TimeConsole_Interval = 1.0e-4;

	std::string sDescription = "";
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time (&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer,80,"%d-%m-%Y %H:%M:%S",timeinfo);
		std::string strTime(buffer);

		sDescription += "-------------------------------------------------------------\n";
		sDescription += "CPDI formulation, Plain Stress, Half Ring Xiang (2017) ------\n";
		sDescription += "Process started on: " + strTime + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Number of threads: " + Script(_MAX_N_THREADS) + "\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Material Point count: " + Script(allMaterialPoint_CPDI.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3_Cells.x) + "," + Script(i3_Cells.y) + "," + Script(i3_Cells.z) + ")\n";
//		sDescription += "Kernel Resolution: (" + Script(i3_Cells_Kernel.x) + "," + Script(i3_Cells_Kernel.y) + "," + Script(i3_Cells_Kernel.z) + ")\n";
		sDescription += "Division (Angular): " + Script(iDivision_Angular) + "\n";
		sDescription += "Division (Radial): " + Script(iDivision_Radial) + "\n";
		sDescription += "Division (Longitudinal): " + Script(iDivision_Longitudinal) + "\n";
		sDescription += "Tube length: " + Script(dLength_Ring,3) + "\n";
		sDescription += "Timeline Speed: " + Script(m_TimeLine.getVelocity(1.0e-4).y, 3) + " m/s" + "\n";
		sDescription += "Yield: " + Script(allMaterialPoint_CPDI[0]->d_YieldStress, 3) + " N/m^2" + "\n";
		sDescription += "Modulus: " + Script(allMaterialPoint_CPDI[0]->d_ElasticModulus, 3) + " N/m^2" + "\n";
		sDescription += "Hardening 0: " + Script(allMaterialPoint_CPDI[0]->d_Hardening_Isotropic_C0, 3) + "\n";
		sDescription += "Hardening 1: " + Script(allMaterialPoint_CPDI[0]->d_Hardening_Isotropic_C1, 3) + "\n";

		sDescription += "Global Damping: " + Script(d_DampingCoefficient, 3) + "\n";
		sDescription += "Non-slip contact\n";
//		sDescription += "Elastic-Perfectly plastic\n";
	}

	d_TimeConsole_Last = 0.0;
	this->reportConsole(sDescription);
}
// ----------------------------------------------------------------------------
