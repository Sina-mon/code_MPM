#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Classic_Foam(void)
{
	MaterialPoint_Factory_Classic_CC	MP_Factory;
	GridPoint_Factory					GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	glm::dvec3 d3_Length_Grid = glm::dvec3(0.025, 0.025, 0.001/10.0);
	glm::ivec3 i3_Cells = glm::ivec3(10.0*25, 10.0*25, 1);
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
		if(fabs(dx - d3_Length_Grid.x) < 1.5*d3_Length_Cell.x)
		{
			//thisGridPoint->b3_Fixed.x = true;
		}
		if(fabs(dy - 0.0) < 1.5*d3_Length_Cell.y)
		{
//			thisGridPoint->b3_Fixed.y = true;
//			thisGridPoint->b3_Fixed.z = true;
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

	Material_BC *pAluminum = new Material_BC;
	v_allMaterial.push_back(pAluminum);
	{
		pAluminum->i_ID = 0;
		pAluminum->i_MaterialType = __VONMISESHARDENING;

		pAluminum->d_Density = 2760.0;

		pAluminum->d_ElasticModulus = 70.0e9;
		pAluminum->d_PoissonRatio = 0.3;

		pAluminum->d_YieldStress = 190.0e6;
		pAluminum->d_Hardening_Isotropic_C0 = 30.0;
		pAluminum->d_Hardening_Isotropic_C1 = 50.0e6;
	}
	Material_BC *pSteel = new Material_BC;
	v_allMaterial.push_back(pSteel);
	{
		pSteel->i_ID = 0;
		pSteel->i_MaterialType = __ELASTIC;

		pSteel->d_Density = 7800.0;

		pSteel->d_ElasticModulus = 210.0e9;
		pSteel->d_PoissonRatio = 0.3;
	}

	double dPlatenSpeed = +10.0;
	double dThickness_Min = 0.8*0.085e-3;
	double dThickness_Max = 0.8*0.085e-3;
	double dDiameter_Min = 0.00;
	double dDiameter_Max = 0.006;
	double dOffset = dThickness_Min/4.0;

	glm::dvec3 d3Dimension_Platen_Top	= glm::dvec3(0.020,2.0*d3_Length_Cell.y,dOffset);
	glm::dvec3 d3Center_Platen_Top		= glm::dvec3(20.0*d3_Length_Cell.x+0.5*d3Dimension_Platen_Top.x, 0.02 + 0.5*d3Dimension_Platen_Top.y + 2.0*d3_Length_Cell.y, 0.5*d3_Length_Grid.z);

	double dRelativeDensity = 0.0;
	if(false)
	{// sample, non-overlapping circular pores with minimum thickness
		Canvas_CC Canvas(glm::dvec3(0.020, 0.020, dOffset), dOffset);

		Canvas.drawRectangle(0.5*Canvas.d3_Size, Canvas.d3_Size, glm::dvec3(0.0,0.0,0.0));

		std::vector<glm::dvec4> vPores;
		for(int iCount = 0; iCount < 5000; iCount++)
		{
			glm::dvec3 d3Center = glm::dvec3(_RANDOM(0.0,1.0*Canvas.d3_Size.x,1000),_RANDOM(0.0,1.0*Canvas.d3_Size.y,1000), 0.0);
			double dRadius = _RANDOM(4.0*dOffset,0.5*0.003,1000);

			bool bIntersect = false;
			for(int index_Pores = 0; index_Pores < vPores.size(); index_Pores++)
			{
				double dDistance = glm::length(glm::dvec3(vPores[index_Pores]) - d3Center);

				if(dDistance < dRadius + vPores[index_Pores].w + dThickness_Min)
					bIntersect  = true;
			}
			if(bIntersect == false)
			{
				vPores.push_back(glm::dvec4(d3Center, dRadius));
				Canvas.cutCircle(d3Center, dRadius);
			}
		}

		std::vector<Voxel_ST> vVoxels = Canvas.getVoxels(true);
		std::cout << "Canvas voxels: " << Script(Canvas.v_Voxels.size()) << " (" << Script(Canvas.i3_Size.x) << "," << Script(Canvas.i3_Size.y) << "," << Script(Canvas.i3_Size.z) << ")" << std::endl;
		std::cout << "Active voxels: " << Script(vVoxels.size()) << std::endl;
		dRelativeDensity = (double)vVoxels.size()/Canvas.v_Voxels.size();
		std::cout << "Ratio: " << Script(dRelativeDensity,3) << std::endl;

		for(unsigned int index_Voxel = 0; index_Voxel < vVoxels.size(); index_Voxel++)
		{// get canvas voxels and create material points
			MaterialPoint_BC *newMP;
			newMP = MP_Factory.createMaterialPoint(vVoxels[index_Voxel].d3_Position);

			newMP->i_ID = vVoxels[index_Voxel].i_ID;
			newMP->p_Material = pAluminum;

			newMP->d_Volume_Initial = dOffset*dOffset*dOffset;
			newMP->d_Volume = newMP->d_Volume_Initial;

			newMP->d_Mass = newMP->p_Material->d_Density * newMP->d_Volume;

			newMP->d3_Position = vVoxels[index_Voxel].d3_Position + glm::dvec3(2.0*d3_Length_Cell.x-dOffset, 2.0*d3_Length_Cell.y-dOffset, 0.5*d3_Length_Grid.z);
			newMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			newMP->d3_Force_External = newMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);

			allMaterialPoint.push_back(newMP);
			// log
			v_MarkedMaterialPoints_Momentum.push_back((MaterialPoint_Classic_CC *)newMP);
			v_MarkedMaterialPoints_Principal_Monitor.push_back((MaterialPoint_Classic_CC *)newMP);
		}
	}

	if(true)
	{// sample, voronoi cells
		Canvas_CC Canvas(glm::dvec3(0.020, 0.020, dOffset), dOffset);
		Canvas.drawRectangle(0.5*Canvas.d3_Size, Canvas.d3_Size, glm::dvec3(0.0,0.0,0.0));
		Canvas.cutRectangle(0.5*Canvas.d3_Size, Canvas.d3_Size-glm::dvec3(2.0*dThickness_Min,2.0*dThickness_Min,0.0), glm::dvec3(0.0,0.0,0.0));

		std::vector<glm::dvec4> vPores;
		for(int iCount = 0; iCount < 2000; iCount++)
		{
			glm::dvec3 d3Center = glm::dvec3(_RANDOM(0.0,1.0*Canvas.d3_Size.x,1000),_RANDOM(0.0,1.0*Canvas.d3_Size.y,1000), 0.0);
			double dRadius = _RANDOM(0.5*dDiameter_Min,0.5*dDiameter_Max,1000);
//			double dRadius = 0.5*0.003;
//			double dRadius = Canvas.d3_Size.x / glm::sqrt(20);

			bool bIntersect = false;
			for(int index_Pores = 0; index_Pores < vPores.size(); index_Pores++)
			{
				double dDistance = glm::length(glm::dvec3(vPores[index_Pores]) - d3Center);

				if(dDistance < dRadius + vPores[index_Pores].w + dThickness_Min)
					bIntersect  = true;
			}
			if(bIntersect == false)
			{
				vPores.push_back(glm::dvec4(d3Center, dRadius));
//				Canvas.drawRing(d3Center, dRadius, 0.5*dRadius);
//				Canvas.drawRing(d3Center, 4.0*dOffset, 3.0*dOffset);
			}

//			vPores.push_back(glm::dvec4(d3Center, dRadius));

//			Canvas.drawRing(d3Center,4.0*dOffset,2.0*dOffset);
		}

		for(int index_Voxel = 0; index_Voxel < Canvas.v_Voxels.size(); index_Voxel++)
		{
			double vDistance_Min[2] = {1.0,1.0};
			glm::dvec4 vPore_Min[2];
			for(int index_BP = 0; index_BP < vPores.size(); index_BP++)
			{
				double dDistance = glm::length(glm::dvec2(Canvas.v_Voxels[index_Voxel].d3_Position) - glm::dvec2(vPores[index_BP]));

				if(dDistance < vDistance_Min[0])
				{
					vDistance_Min[1] = vDistance_Min[0];
					vDistance_Min[0] = dDistance;

					vPore_Min[1] = vPore_Min[0];
					vPore_Min[0] = vPores[index_BP];
				}
				else if(vDistance_Min[0] < dDistance && dDistance < vDistance_Min[1])
				{
					vDistance_Min[1] = dDistance;

					vPore_Min[1] = vPores[index_BP];
				}
			}

			double dDistance_BP = glm::length(glm::dvec3(vPore_Min[0]) - glm::dvec3(vPore_Min[1]));

			if(glm::abs(vDistance_Min[0] - vDistance_Min[1]) < dThickness_Min)
			{
				Canvas.v_Voxels[index_Voxel].b_Active = true;
			}
		}

		std::vector<Voxel_ST> vVoxels = Canvas.getVoxels(true);
		std::cout << "Canvas voxels: " << Script(Canvas.v_Voxels.size()) << " (" << Script(Canvas.i3_Size.x) << "," << Script(Canvas.i3_Size.y) << "," << Script(Canvas.i3_Size.z) << ")" << std::endl;
		std::cout << "Active voxels: " << Script(vVoxels.size()) << std::endl;
		dRelativeDensity = (double)vVoxels.size()/Canvas.v_Voxels.size();
		std::cout << "Ratio: " << Script(dRelativeDensity,3) << std::endl;

		for(unsigned int index_Voxel = 0; index_Voxel < vVoxels.size(); index_Voxel++)
		{// get canvas voxels and create material points
			MaterialPoint_BC *newMP;
			newMP = MP_Factory.createMaterialPoint(vVoxels[index_Voxel].d3_Position);

			newMP->i_ID = vVoxels[index_Voxel].i_ID;
			newMP->p_Material = pAluminum;

			newMP->d_Volume_Initial = dOffset*dOffset*dOffset;
			newMP->d_Volume = newMP->d_Volume_Initial;

			newMP->d_Mass = newMP->p_Material->d_Density * newMP->d_Volume;

			newMP->d3_Position = vVoxels[index_Voxel].d3_Position + glm::dvec3(20.0*d3_Length_Cell.x+dOffset, 2.0*d3_Length_Cell.y+dOffset, 0.5*d3_Length_Grid.z);
			newMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			newMP->d3_Force_External = newMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);

			allMaterialPoint.push_back(newMP);
			// log
			v_MarkedMaterialPoints_Momentum.push_back((MaterialPoint_Classic_CC *)newMP);
			v_MarkedMaterialPoints_Principal_Monitor.push_back((MaterialPoint_Classic_CC *)newMP);
		}
	}

	if(true)
	{// top platen material points -------------------------------------------- platen MP
		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center_Platen_Top, d3Dimension_Platen_Top, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

			thisMP->p_Material = pSteel;

			thisMP->i_Body = 1;

			thisMP->d_Volume_Initial = dOffset*dOffset*dOffset;
			thisMP->d_Volume = thisMP->d_Volume_Initial;

			thisMP->d_Mass = thisMP->p_Material->d_Density * thisMP->d_Volume;

			thisMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			thisMP->d3_Force_External = thisMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// send to MP vectors
			MaterialPoint_Classic_CC *thisMP = (MaterialPoint_Classic_CC *)thisMaterialDomain[index_MP];
			// all MPs
			allMaterialPoint.push_back(thisMP);
			// displacement control
//			if(thisMP->d3_Position.y > d3Center_Platen_Top.y + 0.375*d3Dimension_Platen_Top.y)
			{
				thisMP->b_DisplacementControl = true;
				thisMP->f_DisplacementControl_Multiplier = -1.0;
				thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_Displacement_Control.push_back(thisMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(thisMP);
			}
		}
	}

	d_TimeIncrement_Maximum = 1.0/10.0*10.0e-8;
	d_TimeEnd = 0.8*(d3Center_Platen_Top.y) / glm::abs(dPlatenSpeed);
	d_TimeConsole_Interval = 0.05e-3 / glm::abs(dPlatenSpeed);

	// timeline events -------------------------------------------------------
	m_TimeLine.addTimePoint(0.0,	glm::dvec3(0.0, +dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(10.,	glm::dvec3(0.0, +dPlatenSpeed, 0.0));

	double dMass_Domain = 0.0;
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{// calculate debug values
		dMass_Domain += allMaterialPoint[index_MP]->d_Mass;
	}

	a_Runtime.fill(0.0);
	d_DampingCoefficient = 0.0;

	std::string sDescription = "";
	{// description
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time (&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer,80,"%d-%m-%Y %H:%M:%S",timeinfo);
		std::string strTime(buffer);

		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Classic formulation, Ring, Fan (2013) -----------------------\n";
		sDescription += "Process started on: " + strTime + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Number of threads: " + Script(_MAX_N_THREADS) + "\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Material Point count: " + Script(allMaterialPoint.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3_Cells.x) + "," + Script(i3_Cells.y) + "," + Script(i3_Cells.z) + ")" + "(" + Script(d3_Length_Cell.x,3) + ")\n";
		sDescription += "dOffset: " + Script(dOffset,4) + "\n";
		sDescription += "Relative density: " + Script(dRelativeDensity,4) + "\n";
		sDescription += "Thickness_min: " + Script(dThickness_Min,4) + "\n";
		sDescription += "Thickness_max: " + Script(dThickness_Max,4) + "\n";
		sDescription += "Diameter_min: " + Script(dDiameter_Min,4) + "\n";
		sDescription += "Diameter_max: " + Script(dDiameter_Max,4) + "\n";
		sDescription += "Sample depth: " + Script(dOffset,3) + "\n";
		sDescription += "Timeline Speed: " + Script(m_TimeLine.getVelocity(1.0e-4).y, 3) + " m/s" + "\n";
		sDescription += "Yield: " + Script(pAluminum->d_YieldStress, 3) + " N/m^2" + "\n";
		sDescription += "Modulus: " + Script(pAluminum->d_ElasticModulus, 3) + " N/m^2" + "\n";
		sDescription += "Hardening0: " + Script(pAluminum->d_Hardening_Isotropic_C0, 3) + "\n";
		sDescription += "Hardening1: " + Script(pAluminum->d_Hardening_Isotropic_C1, 3) + "\n";

		sDescription += "Global Damping: " + Script(d_DampingCoefficient, 3) + "\n";
		sDescription += "Non-slip contact\n";
//		sDescription += "Elastic-Perfectly plastic\n";
	}

	{// save to file
		d_TimeConsole_Last = 0.0;
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time (&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer,80,"%Y%m%d%H%M%S",timeinfo);
		std::string strTime(buffer);

		str_Log_FileName = _STR_LOGFILE + strTime + ".txt";
		str_Snapshot_FileName = _STR_SNAPFILE + strTime + "_";
		this->reportConsole(sDescription);
	}
}
// ----------------------------------------------------------------------------
