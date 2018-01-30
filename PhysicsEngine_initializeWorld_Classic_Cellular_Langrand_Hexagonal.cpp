#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Classic_Cellular_Langrand_Hexagonal(void)
{
	MaterialPoint_Factory_Classic_CC	MP_Factory;
	GridPoint_Factory					GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	glm::dvec3 d3Length_Grid = glm::dvec3(0.010, 0.010, 0.002/10.0);
	glm::ivec3 i3Cells = glm::ivec3(10.0*10, 10.0*10, 2);
//	glm::dvec3 d3Length_Grid = glm::dvec3(0.030, 0.030, 0.001/10.0);
//	glm::ivec3 i3Cells = glm::ivec3(10.0*30, 10.0*30, 1);
	glm::dvec3 d3Length_Cell = d3Length_Grid / glm::dvec3(i3Cells);
	glm::ivec3 i3Nodes = i3Cells + glm::ivec3(1, 1, 1);
	this->d3_Length_World = d3Length_Grid;

	for(int indexThread = 0; indexThread < _MAX_N_THREADS; indexThread++)
	{// initialize GP mediator
		mpm_GP_Mediator_Thread[indexThread].d3_Length_Grid = d3Length_Grid;
		mpm_GP_Mediator_Thread[indexThread].d3_Length_Cell = d3Length_Cell;
		mpm_GP_Mediator_Thread[indexThread].i3_Cells = i3Cells;
		mpm_GP_Mediator_Thread[indexThread].i3_Node_Count = i3Nodes;

		allGridPoint = GP_Factory.createGrid(d3Length_Grid, i3Cells);
	}


	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{// grid point boundary conditions
		GridPoint *thisGridPoint = allGridPoint[index_GP];

		double dx = thisGridPoint->d3_Position[0];
		double dy = thisGridPoint->d3_Position[1];
		double dz = thisGridPoint->d3_Position[2];

		//fixed grid points
		thisGridPoint->b3_Fixed = glm::bvec3(false, false, false);

		if(fabs(dx - 0.0) < 1.5*d3Length_Cell.x)
		{
			thisGridPoint->b3_Fixed.x = true;
		}
		if(fabs(dy - 0.0) < 1.5*d3Length_Cell.y)
		{
			thisGridPoint->b3_Fixed.y = true;
//			thisGridPoint->b3_Fixed.z = true;
//			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dz - 0.0) < 2.0*d3Length_Grid.z)
		{
			thisGridPoint->b3_Fixed.z = true;
		}
	}

	// multi-body implementation
	for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
	{
		allGridPoint_Body[index_Body] = GP_Factory.createGrid(d3Length_Grid, i3Cells);
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

	double dThickness_Ring = 0.0005;
	double dDiameter_Outer = 0.005;
	double dDiameter_Inner = dDiameter_Outer - 2.0*dThickness_Ring;
	double dDiameter_Average = 0.5*(dDiameter_Outer + dDiameter_Inner);

	glm::ivec2 i2Array_Count = glm::ivec2(1,2);
	glm::dvec2 d2Array_Offset = glm::dvec2(dDiameter_Outer, dDiameter_Outer*glm::sin(_PI/3.0));

	double dPlatenSpeed = +10.0*i2Array_Count.y;
//	double dOffset = glm::min(dThickness_Ring/16.0, d3Length_Cell.x/8.0);
	double dOffset = d3Length_Cell.x/4.0;
//	double dOffset = dThickness_Ring/16.0;

	double dRadius_Inner = 0.5*dDiameter_Inner;
	double dRadius_Outer = 0.5*dDiameter_Outer;
	double dLength_Ring = glm::min(dOffset,d3Length_Cell.z);

	glm::dvec3 d3Dimension_Platen_Top		= glm::dvec3(i2Array_Count.x*dDiameter_Outer, 8.0*d3Length_Cell.y, dOffset);

	glm::dvec3 d3Center_Array			= glm::dvec3(0.5*dDiameter_Outer, 0.0, 0.5*d3Length_Grid.z);
	glm::dvec3 d3Center_Platen_Top		= glm::dvec3(0.5*d3Dimension_Platen_Top.x, d3Center_Array.y + (i2Array_Count.y-1.0)*d2Array_Offset.y+0.5*dDiameter_Outer+0.5*d3Dimension_Platen_Top.y + 1.5*d3Length_Cell.y, 0.5*d3Length_Grid.z);

	if(true)
	{// sample
		Canvas2D_CC Canvas(glm::dvec2(d3Length_Grid.x, d3Length_Grid.y), dOffset);
		for(int ix = 0; ix < i2Array_Count.x; ix++)
		{// sample array
			for(int iy = 0; iy < i2Array_Count.y; iy++)
			{
				glm::dvec2 d2Center_Ring = glm::dvec2(0.0,0.0);

				if(iy % 2 == 0)
					d2Center_Ring = glm::dvec2(d3Center_Array) + glm::dvec2(ix*d2Array_Offset.x, iy*d2Array_Offset.y);
				else if(iy % 2 == 1)
					d2Center_Ring = glm::dvec2(d3Center_Array) - glm::dvec2(0.5*dDiameter_Outer,0.0) + glm::dvec2(ix*d2Array_Offset.x, iy*d2Array_Offset.y);

				if(0.5*ix + iy > 2)
					d2Center_Ring += glm::dvec2(0.001,0.0);

				std::vector<double> vAngle = {0.0, 1.0*_PI/3.0, 2.0*_PI/3.0, 3.0*_PI/3.0, 4.0*_PI/3.0, 5.0*_PI/3.0};
				if(iy == i2Array_Count.y-1)
				{
					vAngle.clear();
					vAngle = {0.0, 0.5*_PI, 1.0*_PI, 4.0*_PI/3.0, 5.0*_PI/3.0};
				}
				Canvas.drawRing(d2Center_Ring, dRadius_Outer, dRadius_Inner);
				{// braze
//					for(double dAngle = 0.0; dAngle < 2.0*_PI; dAngle += _PI/3.0)
					for(int index_Angle = 0; index_Angle < vAngle.size(); index_Angle++)
					{
						double dAngle = vAngle[index_Angle];
						glm::dvec2 d2Center_Braze = d2Center_Ring + glm::dvec2(0.5*dDiameter_Average*glm::cos(dAngle), 0.5*dDiameter_Average*glm::sin(dAngle));
						//Canvas.drawRectangle(d3Center_Braze, glm::dvec3(0.0015,dThickness_Ring,d3Length_Cell.z), glm::dvec3(0.0,0.0,dAngle+0.5*_PI));
					}
				}
			}
		}
		std::vector<Voxel_ST *> vVoxels = Canvas.getVoxels(true);
		std::cout << "Canvas voxels: " << Script(Canvas.v_Voxels.size()) << " (" << Script(Canvas.u2_Size.x) << "," << Script(Canvas.u2_Size.y) << ")" << std::endl;
		std::cout << "Active voxels: " << Script(vVoxels.size()) << std::endl;

		for(unsigned int index_Voxel = 0; index_Voxel < vVoxels.size(); index_Voxel++)
		{// get canvas voxels and create material points
			MaterialPoint_BC *newMP;
			newMP = MP_Factory.createMaterialPoint(glm::dvec3(vVoxels[index_Voxel]->d2_Position,0.0));

			newMP->i_ID = vVoxels[index_Voxel]->u_ID;
			newMP->b_Surface = vVoxels[index_Voxel]->b_Surface;
			newMP->p_Material = pInconel;

			newMP->d_Volume_Initial = dOffset*dOffset*dOffset;
			newMP->d_Volume = newMP->d_Volume_Initial;

			newMP->d_Mass = newMP->p_Material->d_Density * newMP->d_Volume;

			newMP->d3_Position = glm::dvec3(vVoxels[index_Voxel]->d2_Position,0.0) + glm::dvec3(d3Length_Cell.x+dOffset, d3Length_Cell.y+dOffset, 0.5*d3Length_Grid.z);
			newMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			newMP->d3_Force_External = newMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);

			allMaterialPoint.push_back(newMP);
			// moment log
			v_MarkedMaterialPoints_Momentum.push_back((MaterialPoint_Classic_CC *)newMP);
		}
	}

	if(true)
	{// platen
		Canvas2D_CC Canvas(glm::dvec2(d3Length_Grid.x, d3Length_Grid.y), dOffset);

		Canvas.drawRectangle(glm::dvec2(d3Center_Platen_Top), glm::dvec2(d3Dimension_Platen_Top), 0.0);

		std::vector<Voxel_ST *> vVoxels = Canvas.getVoxels(true);

		for(unsigned int index_Voxel = 0; index_Voxel < vVoxels.size(); index_Voxel++)
		{// get canvas voxels and create material points
			MaterialPoint_BC *newMP;
			newMP = MP_Factory.createMaterialPoint(glm::dvec3(vVoxels[index_Voxel]->d2_Position,0.0));

			newMP->i_ID = vVoxels[index_Voxel]->u_ID;
			newMP->b_Surface = vVoxels[index_Voxel]->b_Surface;
			newMP->p_Material = pInconel;

			newMP->d_Volume_Initial = dOffset*dOffset*dOffset;
			newMP->d_Volume = newMP->d_Volume_Initial;

			newMP->d_Mass = newMP->p_Material->d_Density * newMP->d_Volume;

			newMP->d3_Position = glm::dvec3(vVoxels[index_Voxel]->d2_Position,0.0) + glm::dvec3(d3Length_Cell.x+dOffset, d3Length_Cell.y+dOffset, 0.5*d3Length_Grid.z);
			newMP->d3_Velocity = glm::dvec3(0.0, -1.0, 0.0);
			newMP->d3_Force_External = newMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);

			allMaterialPoint.push_back(newMP);
			// displacement control
//			if(thisMP->d3_Position.y > d3Center_Platen_Top.y + 0.375*d3Dimension_Platen_Top.y)
			{
//				newMP->b_DisplacementControl = true;
				newMP->b3_DisplacementControl = glm::bvec3(false,true,false);
				newMP->f_DisplacementControl_Multiplier = -1.0;
				newMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_Displacement_Control.push_back(newMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(newMP);
			}
		}
	}

	d_TimeIncrement_Maximum = 1.0/10.0*5.0e-8;
	d_TimeEnd = 0.5*d3Center_Platen_Top.y / glm::abs(dPlatenSpeed);
//	d_TimeEnd = 0.8*dDiameter_Outer / glm::abs(dPlatenSpeed);
	d_TimeConsole_Interval = 0.1e-4 / glm::abs(dPlatenSpeed);

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
		sDescription += "Classic formulation, Cellular, Langrand (2017) --------------\n";
		sDescription += "Process started on: " + strTime + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Number of threads: " + Script(_MAX_N_THREADS) + "\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Material Point count: " + Script(allMaterialPoint.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3Cells.x) + "," + Script(i3Cells.y) + "," + Script(i3Cells.z) + ")" + "(" + Script(d3Length_Cell.x,3) + ")\n";
		sDescription += "dOffset: " + Script(dOffset,4) + "\n";
//		sDescription += "Division (Angular): " + Script(iDivision_Angular) + " (offset: " + Script(_PI*dDiameter_Average/iDivision_Angular,4) + ")" + "\n";
//		sDescription += "Division (Radial): " + Script(iDivision_Radial) + " (offset: " + Script(dThickness_Ring/iDivision_Radial,4) + ")" + "\n";
//		sDescription += "Division (Longitudinal): " + Script(iDivision_Longitudinal) + " (offset: " + Script(dLength_Ring/iDivision_Longitudinal,4) + ")" + "\n";
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
