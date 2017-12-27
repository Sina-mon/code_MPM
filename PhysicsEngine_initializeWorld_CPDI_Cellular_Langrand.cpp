#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_CPDI_Cellular_Langrand(void)
{
	MaterialPoint_Factory_CPDI_CC	MP_Factory;
	GridPoint_Factory				GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	glm::dvec3 d3_Length_Grid = glm::dvec3(0.040, 0.030, 0.001/20.0);
	glm::ivec3 i3_Cells = glm::ivec3(20.0*40, 20.0*30, 1);
	glm::dvec3 d3_Length_Cell = d3_Length_Grid / glm::dvec3(i3_Cells);
	glm::ivec3 i3_Nodes = i3_Cells + glm::ivec3(1, 1, 1);

	for(int indexThread = 0; indexThread < _MAX_N_THREADS; indexThread++)
	{// initialize GP mediator

		d3_Length_World = d3_Length_Grid;

		mpm_GP_Mediator_Thread[indexThread].d3_Length_Grid = d3_Length_Grid;
		mpm_GP_Mediator_Thread[indexThread].d3_Length_Cell = d3_Length_Cell;
		mpm_GP_Mediator_Thread[indexThread].i3_Cells = i3_Cells;
		mpm_GP_Mediator_Thread[indexThread].i3_Node_Count = i3_Nodes;

		//allGridPoint_Thread[iThread] = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);
	}
	allGridPoint = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);

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
		if(fabs(dy - 0.0) < 1.5*d3_Length_Cell.y)
		{
			thisGridPoint->b3_Fixed.y = true;
//			thisGridPoint->b3_Fixed.z = true;
//			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
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
			thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - d3_Length_Grid.z) < dTolerance)
		{
//			thisGridPoint->b3_Fixed.z = true;
		}
	}

	// multi-body implementation
	for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
	{
		allGridPoint_Body[index_Body] = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);
	}
	// locks on grid points for atomic operations
	v_GridPoint_Lock.resize(allGridPoint.size());
	for(int index = 0; index < v_GridPoint_Lock.size(); index++)
	{
		v_GridPoint_Lock[index] = new omp_lock_t;
		omp_init_lock(v_GridPoint_Lock[index]);
	}

	Material_BC *pInconel = new Material_BC;
	v_allMaterial.push_back(pInconel);
	{
		pInconel->i_ID = 0;
		pInconel->i_MaterialType = __VONMISESHARDENING;

		pInconel->d_Density = 8250.0;

		pInconel->d_ElasticModulus = 197.6e9;
		pInconel->d_PoissonRatio = 0.29;

		pInconel->d_YieldStress = 280.0e6;
		pInconel->d_Hardening_Isotropic_C0 = 15.0;
		pInconel->d_Hardening_Isotropic_C1 = 420.0e6;
	}
	Material_BC *pSteel = new Material_BC;
	v_allMaterial.push_back(pSteel);
	{
		pSteel->i_ID = 0;
		pSteel->i_MaterialType = __ELASTIC;

		pSteel->d_Density = 10.0*7800.0;

		pSteel->d_ElasticModulus = 10.0*210.0e9;
		pSteel->d_PoissonRatio = 0.3;
	}

	double dPlatenSpeed = +5.0;
//	ET03
	double dDiameter_Average = 0.0045;
	double dThickness_Ring = 0.0005;
//	ET04
//	double dDiameter_Average = 0.04987;
//	double dThickness_Ring = 0.00296;

	double dDiameter_Outer = dDiameter_Average + dThickness_Ring;
	double dDiameter_Inner = dDiameter_Average - dThickness_Ring;

	glm::ivec2 i2Array_Count = glm::ivec2(1,1);
	glm::dvec2 d2Array_Offset = glm::dvec2(dDiameter_Outer, dDiameter_Outer);

	int iDivision_Angular = 360;
	int iDivision_Radial = 16;
	int iDivision_Longitudinal = 1;

	double dAngle_Start	= -0.0*_PI;
	double dAngle_End	= +2.0*_PI;
	double dRadius_Inner = 0.5*dDiameter_Inner;
	double dRadius_Outer = 0.5*dDiameter_Outer;
	double dLength_Ring = glm::min(dThickness_Ring/iDivision_Radial, 2.0*_PI*dDiameter_Average/iDivision_Angular);
	if(dLength_Ring > d3_Length_Cell.z)
		dLength_Ring = d3_Length_Cell.z;

	glm::dvec3 d3Dimension_Platen_Bottom	= glm::dvec3(0.8*d3_Length_World.x,.0*d3_Length_Cell.y,d3_Length_Grid.z);
	glm::dvec3 d3Dimension_Platen_Top		= glm::dvec3(i2Array_Count.x*dDiameter_Outer,2.0*d3_Length_Cell.y,dLength_Ring);

	glm::dvec3 d3Center_Platen_Bottom	= 0.5*d3Dimension_Platen_Bottom + glm::dvec3(0.0,0.2*d3_Length_Cell.y,0.0);
	glm::dvec3 d3Center_Array			= glm::dvec3(0.5*dDiameter_Outer + 1.0*d3_Length_Cell.x, 0.5*dDiameter_Outer+d3Dimension_Platen_Bottom.y+1.0*d3_Length_Cell.y,0.5*d3_Length_Grid.z);
	glm::dvec3 d3Center_Platen_Top		= glm::dvec3(0.5*d3_Length_Cell.x+0.5*d3Dimension_Platen_Top.x, d3Center_Array.y + (i2Array_Count.y-0.5)*d2Array_Offset.y+0.5*d3Dimension_Platen_Top.y+0.0*d3_Length_Cell.y,0.5*d3_Length_Grid.z);

	for(int ix = 0; ix < i2Array_Count.x; ix++)
	{
		for(int iy = 0; iy < i2Array_Count.y; iy++)
		{
			if(true)
			{// ring material points -------------------------------------------------- tube MP
				glm::dvec3 d3Center_Ring = d3Center_Array + glm::dvec3(ix*d2Array_Offset.x,iy*d2Array_Offset.y,0.0);

				std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Arc(d3Center_Ring, glm::dvec3(0.0,0.0,0.0), dAngle_Start, dAngle_End, dRadius_Outer, dRadius_Inner, dLength_Ring, iDivision_Angular, iDivision_Radial, iDivision_Longitudinal);
				for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
				{// assign material point initial values
					MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

					thisMP->p_Material = pInconel;

					thisMP->i_Body = 0;

					thisMP->d_Volume_Initial = MP_Factory.getVolume((MaterialPoint_CPDI_CC *)thisMP);
					thisMP->d_Volume = thisMP->d_Volume_Initial;

					thisMP->d_Mass = thisMP->p_Material->d_Density * thisMP->d_Volume;
					d_Mass_Minimum = 0.0 * thisMP->d_Mass;

					thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
					thisMP->d3_Force_External = thisMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);

					// identify surface MPs
//					if(glm::length(thisMP->d3_Position-d3Center_Ring) < (dRadius_Inner+dThickness_Ring/iDivision_Radial))
//						thisMP->b_Surface = true;
//					if(glm::length(thisMP->d3_Position-d3Center_Ring) > (dRadius_Outer-dThickness_Ring/iDivision_Radial))
//						thisMP->b_Surface = true;
				}
				for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
				{// send to MP vectors
					MaterialPoint_CPDI_CC *thisMP = (MaterialPoint_CPDI_CC *)thisMaterialDomain[index_MP];
					// all MPs
					allMaterialPoint_CPDI.push_back(thisMP);
					// moment log
					v_MarkedMaterialPoints_Momentum.push_back(thisMP);
					// mark for principal stress/strain monitor
					v_MarkedMaterialPoints_Principal_Monitor.push_back(thisMP);
					// mark for energy monitor
		//			v_MarkedMaterialPoints_Monitor_Energy.push_back(thisMP);
				}
			}
		}
	}

	if(true)
	{// top platen material points -------------------------------------------- platen MP
		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center_Platen_Top, d3Dimension_Platen_Top, 2.0*dLength_Ring);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

			thisMP->p_Material = pInconel;

			thisMP->i_Body = 1;

			thisMP->d_Volume_Initial = MP_Factory.getVolume((MaterialPoint_CPDI_CC *)thisMP);
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			thisMP->d_Mass = thisMP->p_Material->d_Density * thisMP->d_Volume;

			thisMP->d3_Velocity = glm::dvec3(0.0, -dPlatenSpeed, 0.0);
			thisMP->d3_Force_External = thisMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint_CPDI_CC *thisMP = (MaterialPoint_CPDI_CC *)thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint_CPDI.push_back(thisMP);
			// displacement control
//			if(thisMP->d3_Position.y > d3Center_Platen_Top.y + 0.375*d3Dimension_Platen_Top.y)
			{
				thisMP->b_DisplacementControl = true;
				thisMP->f_DisplacementControl_Multiplier = -1.0;
				thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_CPDI_Displacement_Control.push_back(thisMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(thisMP);
			}
		}
	}
	if(false)
	{// bottom platen material points ----------------------------------------- platen MP
//		glm::dvec3 d3Dimension = glm::dvec3(0.8*d3_Length_World.x,4.0*d3_Length_Cell.y,d3_Length_Grid.z);
//		glm::dvec3 d3Center = d3Center_Ring;
//		d3Center.x = 0.5*d3Dimension.x;
//		d3Center.y = 0.5*d3Dimension.y + 0.2*d3_Length_Cell.y;
//		d3Center.z = 0.5*d3Dimension.z;

		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center_Platen_Bottom, d3Dimension_Platen_Bottom, 1.0*d3_Length_Cell.y);
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
			thisMP->d_PoissonRatio = 0.3;
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
			// displacement control
			if(false)
//			if(thisMP->d3_Position.y < d3Center.y - 0.25*d3Dimension.y)
			{
				thisMP->b_DisplacementControl = true;
				thisMP->f_DisplacementControl_Multiplier = +.0;
				thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_CPDI_Displacement_Control.push_back(thisMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(thisMP);
			}
		}
	}

	d_TimeIncrement_Maximum = 1.0/20.0*5.0e-8;
	d_TimeEnd = (d3Center_Platen_Top.y - 0.2*dDiameter_Outer) / glm::abs(dPlatenSpeed);
//	d_TimeEnd = 0.8*dDiameter_Outer / glm::abs(dPlatenSpeed);
	d_TimeConsole_Interval = 0.2e-4 / glm::abs(dPlatenSpeed);

	// timeline events -------------------------------------------------------
	m_TimeLine.addTimePoint(0.0,					glm::dvec3(0.0, +dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(1.0*d_TimeEnd,			glm::dvec3(0.0, +dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(1.0*d_TimeEnd + 1.0e-6,	glm::dvec3(0.0, -dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(10.,					glm::dvec3(0.0, -dPlatenSpeed, 0.0));

	double dMass_Domain = 0.0;
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint_CPDI.size(); index_MP++)
	{// calculate debug values
		dMass_Domain += allMaterialPoint_CPDI[index_MP]->d_Mass;
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
		sDescription += "CPDI formulation, Cellular, Langrand (2017) ------\n";
		sDescription += "Process started on: " + strTime + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Number of threads: " + Script(_MAX_N_THREADS) + "\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Material Point count: " + Script(allMaterialPoint_CPDI.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3_Cells.x) + "," + Script(i3_Cells.y) + "," + Script(i3_Cells.z) + ")" + "(" + Script(d3_Length_Cell.x,3) + ")\n";
//		sDescription += "Kernel Resolution: (" + Script(i3_Cells_Kernel.x) + "," + Script(i3_Cells_Kernel.y) + "," + Script(i3_Cells_Kernel.z) + ")\n";
		sDescription += "Division (Angular): " + Script(iDivision_Angular) + " (offset: " + Script(0.5*_PI*dDiameter_Average/iDivision_Angular,4) + ")" + "\n";
		sDescription += "Division (Radial): " + Script(iDivision_Radial) + " (offset: " + Script(dThickness_Ring/iDivision_Radial,4) + ")" + "\n";
		sDescription += "Division (Longitudinal): " + Script(iDivision_Longitudinal) + " (offset: " + Script(dLength_Ring/iDivision_Longitudinal,4) + ")" + "\n";
		sDescription += "Tube average_diameter: " + Script(dDiameter_Average,3) + "\n";
		sDescription += "Tube thickness: " + Script(dThickness_Ring,3) + "\n";
		sDescription += "Tube length: " + Script(dLength_Ring,3) + "\n";
		sDescription += "Timeline Speed: " + Script(m_TimeLine.getVelocity(1.0e-4).y, 3) + " m/s" + "\n";
		sDescription += "Yield: " + Script(pInconel->d_YieldStress, 3) + " N/m^2" + "\n";
		sDescription += "Modulus: " + Script(pInconel->d_ElasticModulus, 3) + " N/m^2" + "\n";
		sDescription += "Hardening0: " + Script(pInconel->d_Hardening_Isotropic_C0, 3) + "\n";
		sDescription += "Hardening1: " + Script(pInconel->d_Hardening_Isotropic_C1, 3) + "\n";

		sDescription += "Global Damping: " + Script(d_DampingCoefficient, 3) + "\n";
		sDescription += "Non-slip contact\n";
//		sDescription += "Elastic-Perfectly plastic\n";
	}

	d_TimeConsole_Last = 0.0;
	this->reportConsole(sDescription);
}
// ----------------------------------------------------------------------------
