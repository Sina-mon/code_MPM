#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Classic_Bar(void)
{
	MaterialPoint_Factory_Classic_CC	MP_Factory;
	GridPoint_Factory					GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	glm::dvec3 d3Length_Grid = glm::dvec3(0.025, 0.025, 0.025);
	glm::ivec3 i3Cells = glm::ivec3(1.0*25, 1.0*25, 1*25);
	glm::dvec3 d3Length_Cell = d3Length_Grid / glm::dvec3(i3Cells);
	glm::ivec3 i3Nodes = i3Cells + glm::ivec3(1, 1, 1);

	for(int indexThread = 0; indexThread < _MAX_N_THREADS; indexThread++)
	{// initialize GP mediator

		d3_Length_World = d3Length_Grid;

		mpm_GP_Mediator_Thread[indexThread].d3_Length_Grid = d3Length_Grid;
		mpm_GP_Mediator_Thread[indexThread].d3_Length_Cell = d3Length_Cell;
		mpm_GP_Mediator_Thread[indexThread].i3_Cells = i3Cells;
		mpm_GP_Mediator_Thread[indexThread].i3_Node_Count = i3Nodes;

		//allGridPoint_Thread[iThread] = GP_Factory.createGrid(d3_Length_Grid, i3_Cells);
	}
	allGridPoint = GP_Factory.createGrid(d3Length_Grid, i3Cells);

	for(unsigned int index_GP = 0; index_GP < allGridPoint.size(); index_GP++)
	{// grid point boundary conditions
		GridPoint *thisGridPoint = allGridPoint[index_GP];

		double dx = thisGridPoint->d3_Position[0];
		double dy = thisGridPoint->d3_Position[1];
		double dz = thisGridPoint->d3_Position[2];
		double dTolerance = 0.1*d3Length_Cell.x;

		//fixed grid points
		thisGridPoint->b3_Fixed = glm::bvec3(false, false, false);

		if(fabs(dx - 0.0) < 1.5*d3Length_Cell.x)
		{
			//thisGridPoint->b3_Fixed.x = true;
		}
		if(fabs(dx - d3Length_Grid.x) < 1.5*d3Length_Cell.x)
		{
			//thisGridPoint->b3_Fixed.x = true;
		}
		if(fabs(dy - 0.0) < 1.5*d3Length_Cell.y)
		{
			thisGridPoint->b3_Fixed = glm::bvec3(false, true, false);
		}
		if(fabs(dy - d3Length_Grid.y) < dTolerance)
		{
		}
		if(fabs(dz - 0.0) < 0.5*d3Length_Cell.z)
		{
			//thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - 0.0) < 2.0*d3Length_Grid.z)
		{
			//thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - d3Length_Grid.z) < dTolerance)
		{
//			thisGridPoint->b3_Fixed.z = true;
		}
	}

	// multi-body implementation
//	for(int index_Body = 0; index_Body < _MAX_N_BODIES; index_Body++)
//	{
//		allGridPoint_Body[index_Body] = GP_Factory.createGrid(d3Length_Grid, i3Cells);
//	}
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
		Material_BC *thisMaterial = pAluminum;

		thisMaterial->i_ID = 0;
		thisMaterial->i_MaterialType = __VONMISESHARDENING;

		thisMaterial->d_Density = 2700.0;

		thisMaterial->d_ElasticModulus = 70.0e9;
		thisMaterial->d_PoissonRatio = 0.3;

		thisMaterial->d_YieldStress = 70.0e6;
		thisMaterial->d_Hardening_Isotropic_C0 = 35.0;
		thisMaterial->d_Hardening_Isotropic_C1 = 55.0e6;
	}
	Material_BC *pSteel = new Material_BC;
	v_allMaterial.push_back(pSteel);
	{
		Material_BC *thisMaterial = pSteel;

		thisMaterial->i_ID = 0;
		thisMaterial->i_MaterialType = __PLASTIC;
//		thisMaterial->i_MaterialType = __VONMISESHARDENING;

		thisMaterial->d_Density = 7800.0;

		thisMaterial->d_ElasticModulus = 210.0e9;
		thisMaterial->d_PoissonRatio = 0.3;

		thisMaterial->d_YieldStress = 310.0e6;
		thisMaterial->d_Hardening_Isotropic_C0 = 4.0;
		thisMaterial->d_Hardening_Isotropic_C1 = 150.0e6;
	}

	double dOffset = d3Length_Cell.x/8.0;

	double dPlatenSpeed = +1.0;
	double dBar_Length = 0.02;
	double dBar_Width = 0.004;
	double dBar_Thickness = 0.004;//2.0*dOffset;

	glm::dvec2 d2Bar_Center = glm::dvec2(0.5*d3Length_Grid.x, 0.5*dBar_Length);
	glm::dvec2 d2Bar_Dimension = glm::dvec2(dBar_Width, dBar_Length);

	glm::dvec3 d3Dimension_Platen_Top	= glm::dvec3(0.020,2.0*d3Length_Cell.y,dOffset);
	glm::dvec3 d3Center_Platen_Top		= glm::dvec3(20.0*d3Length_Cell.x+0.5*d3Dimension_Platen_Top.x, 0.02 + 0.5*d3Dimension_Platen_Top.y + 2.0*d3Length_Cell.y, 0.5*d3Length_Grid.z);

	Canvas2D_CC Canvas(glm::dvec2(0.020, 0.020), dOffset);
	{// draw setup on canvas
		Canvas.drawRectangle(d2Bar_Center, d2Bar_Dimension, 0.0);

		glm::dvec2 d2Platen_Center = d2Bar_Center+glm::dvec2(0.0, 0.5*dBar_Length);
		glm::dvec2 d2Platen_Size = glm::dvec2(dBar_Width, 4.0*d3Length_Cell.y);
		Canvas.setLoadRectangle(d2Platen_Center, d2Platen_Size, 0.0, true);
	}
	std::vector<Voxel_ST *> vVoxels = Canvas.getVoxels_Active(true);
	std::cout << "Canvas voxels: " << Script(Canvas.v_Voxels.size()) << " (" << Script(Canvas.u2_Size.x) << "," << Script(Canvas.u2_Size.y) << ")" << std::endl;
	std::cout << "Active voxels: " << Script(vVoxels.size()) << std::endl;

	for(unsigned int index_Voxel = 0; index_Voxel < vVoxels.size(); index_Voxel++)
	{// get canvas voxels and create material points
		for(double dz = -0.5*(dBar_Thickness-dOffset); dz <= 0.5*(dBar_Thickness-dOffset); dz += dOffset)
		{
			MaterialPoint_BC *newMP;
			newMP = MP_Factory.createMaterialPoint(glm::dvec3(0.0, 0.0, 0.0));

			newMP->i_ID = vVoxels[index_Voxel]->u_ID;
			newMP->p_Material = pAluminum;

			newMP->d_Volume_Initial = dOffset*dOffset*dOffset;
			newMP->d_Volume = newMP->d_Volume_Initial;

			newMP->d_Mass = newMP->p_Material->d_Density * newMP->d_Volume;

			// position shifts in x, y and z are added here
			glm::dvec3 d3Position_Shift = glm::dvec3(0.0, 1.0*d3Length_Cell.y, 0.5*d3Length_Grid.z);
			newMP->d3_Position = glm::dvec3(vVoxels[index_Voxel]->d2_Position, dz) + d3Position_Shift;
			newMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			newMP->d3_Force_External = newMP->d_Mass * glm::dvec3(0.0, 0.0, 0.0);

			allMaterialPoint.push_back(newMP);

			if(vVoxels[index_Voxel]->b_Load == true)
			{
				newMP->b3_DisplacementControl = glm::bvec3(false,true,false);
				newMP->f_DisplacementControl_Multiplier = 1.0;
				v_MarkedMaterialPoints_Displacement_Control.push_back(newMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(newMP);
			}

			// log
			v_MarkedMaterialPoints_Momentum.push_back((MaterialPoint_Classic_CC *)newMP);
			//v_MarkedMaterialPoints_Principal_Monitor.push_back((MaterialPoint_Classic_CC *)newMP);
			if(glm::abs(newMP->d3_Position.y - 0.5*dBar_Length) < 2.0*d3Length_Cell.y)
			{// stress monitor
				newMP->b_Mark_Stress = true;
				v_MarkedMaterialPoints_Stress_Monitor.push_back(newMP);
			}
		}
	}

	d_TimeIncrement_Maximum = 1.0/1.0*1.0e-8;
	d_TimeEnd = 0.2*(d3Center_Platen_Top.y) / glm::abs(dPlatenSpeed);
	d_TimeConsole_Interval = 1.0e-5 / glm::abs(dPlatenSpeed);

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
		sDescription += "Grid Resolution: (" + Script(i3Cells.x) + "," + Script(i3Cells.y) + "," + Script(i3Cells.z) + ")" + "(" + Script(d3Length_Cell.x,3) + ")\n";
		sDescription += "dOffset: " + Script(dOffset,4) + "\n";
		sDescription += "Bar_Length: " + Script(dBar_Length,4) + "\n";
		sDescription += "Bar_Width: " + Script(dBar_Width,4) + "\n";
		sDescription += "Bar_Thickness: " + Script(dBar_Thickness,4) + "\n";
		sDescription += "Platen Speed: " + Script(dPlatenSpeed, 3) + " m/s" + "\n";
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
