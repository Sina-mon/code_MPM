#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Classic_Cellular_Shim_Square(void)
{
	MaterialPoint_Factory_Classic_CC	MP_Factory;
	GridPoint_Factory					GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.028, 0.050, 0.002/2.0);// for 2x2
//	glm::ivec3 i3_Cells = glm::ivec3(2.0*28, 2.0*50, 2);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.02675, 0.050, 0.002/4.0);// for 2x2
//	glm::ivec3 i3_Cells = glm::ivec3(4.0*26.75, 4.0*50, 2);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.026, 0.050, 0.002/8.0);// for 2x2
//	glm::ivec3 i3_Cells = glm::ivec3(8.0*26.0, 8.0*50, 2);

//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.0545, 0.100, 0.002/2.0);// for 4x4
//	glm::ivec3 i3_Cells = glm::ivec3(2.0*54.5, 2.0*100, 2);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.0525, 0.100, 0.002/4.0);// for 4x4
//	glm::ivec3 i3_Cells = glm::ivec3(4.0*52.5, 4.0*100, 2);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.05175, 0.100, 0.002/8.0);// for 4x4
//	glm::ivec3 i3_Cells = glm::ivec3(8.0*51.75, 8.0*100, 2);

//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.1335, 0.150, 0.002/2.0);// for 10x10 with 1-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(2.0*133.5, 2.0*150, 2);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.0665, 0.150, 0.002/2.0);// for 5x10 with 1-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(2.0*66.5, 2.0*150, 2);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.1390, 0.150, 0.002/2.0);// for 10x10 with 2-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(2.0*139.0, 2.0*150, 2);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.5*0.1380, 0.150, 0.002/2.0);// for 5x10 with 2-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(2.0*0.5*138.0, 2.0*150, 2);

//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.13025, 0.150, 0.001/4.0);// for 10x10 with 1-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(4.0*130.25, 4.0*150, 1);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.5*0.130, 0.150, 0.001/4.0);// for 5x10 with 1-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(4.0*0.5*130.0, 4.0*150, 1);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.065, 0.150, 0.001/4.0);// for 5x10 with 1-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(4.0*65.0, 4.0*150, 1);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.1330, 0.140, 0.001/4.0);// for 10x10 with 2-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(4.0*133.0, 4.0*150, 1);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.5*0.1330, 0.140, 0.001/4.0);// for 5x10 with 2-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(4.0*0.5*133.0, 4.0*150, 1);

	// calculated using 5D + 4C + 2C
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.0688, 0.150, 0.001/1.25);// for 5x10 with 1-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(1.25*68.8, 1.25*150, 1);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.066, 0.150, 0.001/2.5);// for 5x10 with 1-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(2.5*66.0, 2.5*150, 1);
//	glm::dvec3 d3_Length_Grid = glm::dvec3(0.065, 0.150, 0.001/5.0);// for 5x10 with 1-cell gap
//	glm::ivec3 i3_Cells = glm::ivec3(5.0*65.0, 5.0*150, 1);
	glm::dvec3 d3_Length_Grid = glm::dvec3(0.0641, 0.130, 0.001/10.0);// for 5x10 with 1-cell gap
	glm::ivec3 i3_Cells = glm::ivec3(10.0*64.1, 10.0*130, 1);

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
			thisGridPoint->b3_Fixed.x = true;
			//thisGridPoint->b3_Fixed.y = true;
			thisGridPoint->d_FrictionCoefficient = 0.3;
		}
		if(fabs(dy - 0.0) < 1.5*d3_Length_Cell.y)
		{
			thisGridPoint->b3_Fixed.y = true;
//			thisGridPoint->b3_Fixed.z = true;
//			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dz - 0.0) < 0.5*d3_Length_Cell.z)
		{
			thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - 0.0) < 2.0*d3_Length_Grid.z)
		{
			thisGridPoint->b3_Fixed.z = true;
		}
	}

	// multi-body implementation
//	for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
//	{
//		allGridPoint_Body[index_Body] = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);
//	}
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
//		pInconel->i_MaterialType = __PLASTIC;
//		pInconel->i_MaterialType = __ELASTIC;

		pInconel->d_Density = 8250.0;

		pInconel->d_ElasticModulus = 197.6e9;
		pInconel->d_PoissonRatio = 0.29;

		pInconel->d_YieldStress = 350.0e6;
		pInconel->d_Hardening_Isotropic_C0 = 15.0;
		pInconel->d_Hardening_Isotropic_C1 = 350.0e6;
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
	Material_BC *pAluminum = new Material_BC;
	v_allMaterial.push_back(pAluminum);
	{
		pAluminum->i_ID = 0;
		pAluminum->i_MaterialType = __VONMISESHARDENING;
//		pAluminum->i_MaterialType = __PLASTIC;

		pAluminum->d_Density = 2700.0;

		pAluminum->d_ElasticModulus = 70.0e9;
		pAluminum->d_PoissonRatio = 0.3;

		pAluminum->d_YieldStress = 70.0e6;
		pAluminum->d_Hardening_Isotropic_C0 = 35.0;
		pAluminum->d_Hardening_Isotropic_C1 = 55.0e6;
	}

	double dThickness_Ring = 0.00071;
	double dDiameter_Outer = 0.01269;

	double dDiameter_Inner = dDiameter_Outer - 2.0*dThickness_Ring;

	double dDiameter_Average = 0.5*(dDiameter_Outer + dDiameter_Inner);

	glm::ivec2 i2Array_Count = glm::ivec2(5,10);// ------------------------------------------------------------------
	glm::dvec2 d2Array_Offset = glm::dvec2(dDiameter_Outer+1.0*d3_Length_Cell.x, dDiameter_Outer+1.0*d3_Length_Cell.y);

//	double dPlatenSpeed = 2.0*+1.0*i2Array_Count.y;
//	double dOffset = dThickness_Ring/16.0;
	double dPlatenSpeed = 20.0;
	double dOffset = d3_Length_Cell.x/4.0;

	double dRadius_Inner = 0.5*dDiameter_Inner;
	double dRadius_Outer = 0.5*dDiameter_Outer;
	double dLength_Ring = dOffset;
	if(dLength_Ring > d3_Length_Cell.z)
		dLength_Ring = d3_Length_Cell.z;

	glm::dvec3 d3Dimension_Platen_Top		= glm::dvec3(0.98*d3_Length_Grid.x,2.0*d3_Length_Cell.y,dOffset);

	// for half model
	glm::dvec2 d2Center_Array			= glm::dvec2(0.5*dDiameter_Outer + 1.0*d3_Length_Cell.x, 0.5*dDiameter_Outer + 2.0*d3_Length_Cell.y);
	glm::dvec3 d3Center_Platen_Top		= glm::dvec3(0.5*d3Dimension_Platen_Top.x + 1.0*d3_Length_Cell.x, d2Center_Array.y + (i2Array_Count.y-0.5)*d2Array_Offset.y+0.5*d3Dimension_Platen_Top.y+1.0*d3_Length_Cell.y,0.5*d3_Length_Grid.z);
	// for full model
//	glm::dvec2 d2Center_Array			= glm::dvec2(0.5*d3_Length_Grid.x - 0.5*i2Array_Count.x*d2Array_Offset.x + 0.5*d2Array_Offset.x, 0.5*dDiameter_Outer + 2.0*d3_Length_Cell.y);
//	glm::dvec3 d3Center_Platen_Top		= glm::dvec3(0.5*d3_Length_Grid.x, d2Center_Array.y + (i2Array_Count.y-0.5)*d2Array_Offset.y+0.5*d3Dimension_Platen_Top.y+1.0*d3_Length_Cell.y,0.5*d3_Length_Grid.z);

	if(true)
	{// sample
		Canvas2D_CC Canvas(glm::dvec2(d3_Length_Grid.x, d3_Length_Grid.y), dOffset);
		for(int ix = 0; ix < i2Array_Count.x; ix++)
		{// sample array
			for(int iy = 0; iy < i2Array_Count.y; iy++)
			{
				glm::dvec2 d2Center_Ring = d2Center_Array + glm::dvec2(ix*d2Array_Offset.x,iy*d2Array_Offset.y);

				Canvas.drawRing(d2Center_Ring, dRadius_Outer, dRadius_Inner);
				/*
				{// braze
					for(float fAngle = 0.0; fAngle < 2.0*_PI; fAngle += 0.5*_PI)
					{
						glm::dvec3 d3Center_Braze = d3Center_Ring + glm::dvec3(0.5*dDiameter_Average*glm::cos(fAngle), 0.5*dDiameter_Average*glm::sin(fAngle), 0.0);
						Canvas.drawRectangle(d3Center_Braze, glm::dvec3(0.0015,dThickness_Ring,d3_Length_Cell.z), glm::dvec3(0.0,0.0,fAngle+0.5*_PI));
					}
				}
				*/
			}
		}
		std::vector<Voxel_ST *> vVoxels = Canvas.getVoxels_Active(true);
		std::cout << "Canvas voxels: " << Script(Canvas.v_Voxels.size()) << " (" << Script(Canvas.u2_Size.x) << "," << Script(Canvas.u2_Size.y) << ")" << std::endl;
		std::cout << "Active voxels: " << Script(vVoxels.size()) << std::endl;

		allMaterialPoint.reserve(vVoxels.size());

		for(unsigned int index_Voxel = 0; index_Voxel < vVoxels.size(); index_Voxel++)
		{// get canvas voxels and create material points
			MaterialPoint_BC *newMP;
			newMP = MP_Factory.createMaterialPoint(glm::dvec3(vVoxels[index_Voxel]->d2_Position, 0.0));

			newMP->i_ID = vVoxels[index_Voxel]->u_ID;
			newMP->p_Material = pAluminum;

			newMP->d_Volume_Initial = dOffset*dOffset*dOffset;
			newMP->d_Volume = newMP->d_Volume_Initial;

			newMP->d_Mass = newMP->p_Material->d_Density * newMP->d_Volume;

			newMP->d3_Position = glm::dvec3(vVoxels[index_Voxel]->d2_Position, 0.0) + glm::dvec3(0.0,0.0, 0.5*d3_Length_Grid.z);
			newMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			newMP->d3_Force_External = newMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);

			allMaterialPoint.push_back(newMP);
			// moment log
//			v_MarkedMaterialPoints_Momentum.push_back((MaterialPoint_Classic_CC *)newMP);
//			v_MarkedMaterialPoints_Monitor_Energy.push_back((MaterialPoint_Classic_CC *)newMP);
			if(index_Voxel % 10000 == 0)
				std::cout << "here2: " << Script(index_Voxel) << std::endl;

//			delete vVoxels[index_Voxel];
		}
	}

	if(true)
	{// top platen material points -------------------------------------------- platen MP
		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(d3Center_Platen_Top, d3Dimension_Platen_Top, dOffset);
		for(unsigned int index_MP = 0; index_MP < thisMaterialDomain.size(); index_MP++)
		{// assign material point initial values
			MaterialPoint_BC *thisMP = thisMaterialDomain[index_MP];

			thisMP->p_Material = pInconel;

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
				thisMP->b3_DisplacementControl = glm::bvec3(true, true, true);
				thisMP->f_DisplacementControl_Multiplier = -1.0;
				thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_Displacement_Control.push_back(thisMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(thisMP);
			}
		}
	}

	d_TimeIncrement_Maximum = 1.0/10.0*5.0e-8;
	d_TimeEnd = 0.45*d3Center_Platen_Top.y / glm::abs(dPlatenSpeed);
	d_TimeConsole_Interval = 0.01*d_TimeEnd;
	d_TimeSnapshot_Interval = 10.0*d_TimeConsole_Interval;

	// timeline events -------------------------------------------------------
	m_TimeLine.addTimePoint(0.0,					glm::dvec3(0.0, 0.0, 0.0));
	m_TimeLine.addTimePoint(5.0e-5,					glm::dvec3(0.0, +dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(1.0*d_TimeEnd,			glm::dvec3(0.0, +dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(1.0*d_TimeEnd + 1.0e-6,	glm::dvec3(0.0, -dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(10.,					glm::dvec3(0.0, -dPlatenSpeed, 0.0));

	double dMass_Domain = 0.0;
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
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
		sDescription += "Classic formulation, Cellular, Shim (1986) --------------\n";
		sDescription += "Process started on: " + strTime + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Number of threads: " + Script(_MAX_N_THREADS) + "\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Material Point count: " + Script(allMaterialPoint.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3_Cells.x) + "," + Script(i3_Cells.y) + "," + Script(i3_Cells.z) + ")" + "(" + Script(d3_Length_Cell.x,3) + ")\n";
		sDescription += "dOffset: " + Script(dOffset,4) + "\n";
//		sDescription += "Division (Angular): " + Script(iDivision_Angular) + " (offset: " + Script(_PI*dDiameter_Average/iDivision_Angular,4) + ")" + "\n";
//		sDescription += "Division (Radial): " + Script(iDivision_Radial) + " (offset: " + Script(dThickness_Ring/iDivision_Radial,4) + ")" + "\n";
//		sDescription += "Division (Longitudinal): " + Script(iDivision_Longitudinal) + " (offset: " + Script(dLength_Ring/iDivision_Longitudinal,4) + ")" + "\n";
		sDescription += "Tube average_diameter: " + Script(dDiameter_Average,3) + "\n";
		sDescription += "Tube thickness: " + Script(dThickness_Ring,3) + "\n";
		sDescription += "Tube length: " + Script(dLength_Ring,3) + "\n";
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
