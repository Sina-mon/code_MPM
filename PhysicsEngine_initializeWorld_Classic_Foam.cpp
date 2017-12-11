#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Classic_Foam(void)
{
	MaterialPoint_Factory_Classic_CC	MP_Factory;
	GridPoint_Factory					GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	glm::dvec3 d3_Length_Grid = glm::dvec3(0.050, 0.050, 0.004/2.0);
	glm::ivec3 i3_Cells = glm::ivec3(2.0*50, 2.0*50, 2);
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

	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{// grid point boundary conditions
		GridPoint *thisGridPoint = allGridPoint[index_GP];

		double dx = thisGridPoint->d3_Position[0];
		double dy = thisGridPoint->d3_Position[1];
		double dz = thisGridPoint->d3_Position[2];
		double dTolerance = 0.1*d3_Length_Cell.x;

		//fixed grid points
		thisGridPoint->b3_Fixed = glm::bvec3(false, false, false);

		if(fabs(dx - 0.0) < 0.5*d3_Length_Cell.x)
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
		if(fabs(dz - 0.0) < 0.5*d3_Length_Cell.z)
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
	double dDiameter_Average = 0.5*(dDiameter_Inner + dDiameter_Outer);

	// for classic
	double dOffset_Core = 0.5*d3_Length_Cell.y;
	double dOffset_Ring = 1.0/4.0*dThickness_Ring;
	double dOffset_Platen = 0.5*d3_Length_Cell.y;
	// for CPDI
	int iDivision_Angular = 4500;
	int iDivision_Radial = 400;
	int iDivision_Longitudinal = 1;

	double dAngle_Start	= -0.5*_PI;
	double dAngle_End	= +0.5*_PI;
	double dRadius_Inner = 0.5*dDiameter_Inner;
	double dRadius_Outer = 0.5*dDiameter_Outer;
	double dLength_Ring = dOffset_Ring;//glm::min(0.5*d3_Length_Grid.z, 0.5*(dThickness_Ring/iDivision_Radial + 0.5*_PI*dDiameter_Average/iDivision_Angular));

	glm::dvec3 d3Dimension_Core				= glm::dvec3(0.8*d3_Length_World.x,0.8*d3_Length_World.y,dOffset_Core);
	glm::dvec3 d3Dimension_Platen_Top		= glm::dvec3(0.95*d3_Length_World.x,1.0*d3_Length_Cell.y,d3_Length_Grid.z);

	glm::dvec3 d3Center_Core			= glm::dvec3(0.5*d3Dimension_Core.x, 0.5*d3Dimension_Core.y+0.0*d3_Length_Cell.y,0.5*d3Dimension_Core.z);
	glm::dvec3 d3Center_Ring			= glm::dvec3(d3_Length_Cell.x, 0.5*dDiameter_Outer+2.6*d3_Length_Cell.y,0.5*dLength_Ring);
	glm::dvec3 d3Center_Platen_Top		= glm::dvec3(0.5*d3Dimension_Platen_Top.x, d3Center_Core.y+0.5*d3Dimension_Core.y+2.0*d3_Length_Cell.y,0.0);

	if(true)
	{// core ------------------------------------------------------------------ core
		double dGravity = 0.0;

		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center_Core, d3Dimension_Core, dOffset_Core);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

//			thisMP->i_MaterialType = _VONMISESHARDENING;
			thisMP->i_MaterialType = _PLASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = glm::pow(dOffset_Ring, 3);
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 2760.0 * thisMP->d_Volume;
			d_Mass_Minimum = 0.0001 * dMass;
			thisMP->d_Mass = dMass;

			thisMP->d_ElasticModulus = 70.0e9;
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
			MaterialPoint_Classic_CC *thisMP = (MaterialPoint_Classic_CC *)thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint.push_back(thisMP);
			// moment log
//			v_MarkedMaterialPoints_Momentum.push_back(thisMP);
			// mark for stress monitor
//			v_MarkedMaterialPoints_Stress_Monitor.push_back(thisMP);
			// mark for energy monitor
//			v_MarkedMaterialPoints_Monitor_Energy.push_back(thisMP);
		}
	}
	if(true)
	{// random holes
		std::vector<glm::dvec4> vHoles;
		for(int iHole = 0; iHole < 400; iHole++)
		{
			double dx = _RANDOM((d3Center_Core.x-0.5*d3Dimension_Core.x), (d3Center_Core.x+0.5*d3Dimension_Core.x), 100);
			double dy = _RANDOM((d3Center_Core.y-0.5*d3Dimension_Core.y), (d3Center_Core.y+0.5*d3Dimension_Core.y), 100);
//			double dz = _RANDOM(d3Center_Core.z-0.5*d3Dimension_Core.z, d3Center_Core.z+0.5*d3Dimension_Core.z, 100000);
			double dr = _RANDOM(0.0005, 0.001, 1000);

//			glm::dvec4 d4Hole = glm::dvec4(0.025,0.005,0.0,0.002);
			glm::dvec4 d4Hole = glm::dvec4(dx,dy,d3Center_Core.z,dr);

			vHoles.push_back(d4Hole);
		}
		for(int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_Classic_CC *thisMP = (MaterialPoint_Classic_CC *)allMaterialPoint[index_MP];

			for(int iHole = 0; iHole < vHoles.size(); iHole++)
			{
				glm::dvec3 d3Hole_Coordinate = glm::dvec3(vHoles[iHole]);
				double dRadius = vHoles[iHole].w;

				if(glm::length(thisMP->d3_Position - d3Hole_Coordinate) <= dRadius)
				{
//					delete allMaterialPoint[index_MP];
					thisMP->b_Surface = true;
					allMaterialPoint.erase(allMaterialPoint.begin()+index_MP);
					index_MP--;
//					if (index_MP != allMaterialPoint.size() - 1)// Beware of move assignment to self
//						allMaterialPoint[index_MP] = std::move(allMaterialPoint.back());
//					allMaterialPoint.pop_back();

					break;
				}
			}
		}
	}
	if(true)
	{// top platen material points -------------------------------------------- platen MP
		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center_Platen_Top, d3Dimension_Platen_Top, dOffset_Platen);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

			thisMP->i_MaterialType = _ELASTIC;
			thisMP->i_ID = 1;

			thisMP->d_Volume_Initial = glm::pow(dOffset_Platen, 3);
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			double dMass = 7800.0 * thisMP->d_Volume;
			thisMP->d_Mass = dMass;

			thisMP->d_ElasticModulus = 210.0e9;
			thisMP->d_Viscosity = 0.0;
			thisMP->d_PoissonRatio = 0.3;
			thisMP->d_YieldStress = 200.0e6;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			thisMP->d3_Force_External = thisMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint_Classic_CC *thisMP = (MaterialPoint_Classic_CC *)thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint.push_back(thisMP);
			// displacement control
//			if(false)
			//if(thisMP->d3_Position.y > d3Center_Platen_Top.y + 0.25*d3Dimension_Platen_Top.y)
			{
				thisMP->b_DisplacementControl = true;
				thisMP->f_DisplacementControl_Multiplier = -1.0;
				thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_Displacement_Control.push_back(thisMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(thisMP);
			}
		}
	}

	double dPlatenSpeed = +10.0;

	d_TimeIncrement_Maximum = 10.0e-8;
	d_TimeEnd = 0.5*d3Dimension_Core.y / glm::abs(dPlatenSpeed);
	d_TimeConsole_Interval = 0.1e-3 / glm::abs(dPlatenSpeed);

	// timeline events -------------------------------------------------------
	m_TimeLine.addTimePoint(0.0,	glm::dvec3(0.0, +dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(10.0,	glm::dvec3(0.0, +dPlatenSpeed, 0.0));

	double dMass_Domain = 0.0;
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint_CPDI.size(); index_MP++)
	{// calculate debug values
		dMass_Domain += allMaterialPoint[index_MP]->d_Mass;
	}

	a_Runtime.fill(0.0);
	d_DampingCoefficient = 0.0;

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
		sDescription += "Classic formulation, Plain Stress, Half Ring Xiang (2017) ------\n";
		sDescription += "Process started on: " + strTime + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Number of threads: " + Script(_MAX_N_THREADS) + "\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Material Point count: " + Script(allMaterialPoint.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3_Cells.x) + "," + Script(i3_Cells.y) + "," + Script(i3_Cells.z) + ")\n";
//		sDescription += "Kernel Resolution: (" + Script(i3_Cells_Kernel.x) + "," + Script(i3_Cells_Kernel.y) + "," + Script(i3_Cells_Kernel.z) + ")\n";
		sDescription += "Offset_Ring: " + Script(dOffset_Ring,4) + "\n";
		sDescription += "Offset_Ring: " + Script(dOffset_Ring,4) + "\n";
		sDescription += "Offset_Platen: " + Script(dOffset_Platen) + "\n";
//		sDescription += "Division (Angular): " + Script(iDivision_Angular) + " (offset: " + Script(0.5*_PI*dDiameter_Average/iDivision_Angular,4) + ")" + "\n";
//		sDescription += "Division (Radial): " + Script(iDivision_Radial) + " (offset: " + Script(dThickness_Ring/iDivision_Radial,4) + ")" + "\n";
//		sDescription += "Division (Longitudinal): " + Script(iDivision_Longitudinal) + " (offset: " + Script(dLength_Ring/iDivision_Longitudinal,4) + ")" + "\n";
		sDescription += "Tube length: " + Script(dLength_Ring,3) + "\n";
		sDescription += "Timeline Speed: " + Script(m_TimeLine.getVelocity(1.0e-4).y, 3) + " m/s" + "\n";
		sDescription += "Yield: " + Script(allMaterialPoint[0]->d_YieldStress, 3) + " N/m^2" + "\n";
		sDescription += "Modulus: " + Script(allMaterialPoint[0]->d_ElasticModulus, 3) + " N/m^2" + "\n";
		sDescription += "Hardening 0: " + Script(allMaterialPoint[0]->d_Hardening_Isotropic_C0, 3) + "\n";
		sDescription += "Hardening 1: " + Script(allMaterialPoint[0]->d_Hardening_Isotropic_C1, 3) + "\n";

		sDescription += "Global Damping: " + Script(d_DampingCoefficient, 3) + "\n";
		sDescription += "Non-slip contact\n";
//		sDescription += "Elastic-Perfectly plastic\n";
	}

	d_TimeConsole_Last = 0.0;
	this->reportConsole(sDescription);
}
// ----------------------------------------------------------------------------
