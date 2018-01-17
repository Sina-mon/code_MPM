#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Classic_ESO(Canvas2D_CC *pCanvas)
{
	MaterialPoint_Factory_Classic_CC	MP_Factory;
	GridPoint_Factory					GP_Factory;
	// ------------------------------------------------------------------------
	// grid points ------------------------------------------------------------
	double dOffset = pCanvas->d_Offset;
	glm::dvec3 d3Length_Grid = glm::dvec3(pCanvas->d2_Size, 4.0*pCanvas->d_Offset/1.0);

	glm::ivec3 i3Cells = glm::ivec3(d3Length_Grid/(4.0*pCanvas->d_Offset));
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
			thisGridPoint->b3_Fixed.x = true;
			thisGridPoint->b3_Fixed.y = true;
		}
		if(fabs(dx - d3Length_Grid.x) < dTolerance)
		{
		}
		if(fabs(dy - 0.0) < 1.5*d3Length_Cell.y)
		{
			thisGridPoint->b3_Fixed.y = true;
//			thisGridPoint->b3_Fixed.z = true;
//			thisGridPoint->b3_Fixed = glm::bvec3(true, true, true);
		}
		if(fabs(dy - d3Length_Grid.y) < dTolerance)
		{
		}
		if(fabs(dz - 0.0) < 0.5*d3Length_Cell.z)
		{
			thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - 0.0) < 2.0*d3Length_Grid.z)
		{
			thisGridPoint->b3_Fixed.z = true;
		}
		if(fabs(dz - d3Length_Grid.z) < dTolerance)
		{
//			thisGridPoint->b3_Fixed.z = true;
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
//		pInconel->i_MaterialType = __VONMISESHARDENING;
//		pInconel->i_MaterialType = __PLASTIC;
		pInconel->i_MaterialType = __ELASTIC;

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

	if(true)
	{// sample based on canvas
		std::vector<Voxel_ST> vVoxels = pCanvas->getVoxels(true);

		for(unsigned int index_Voxel = 0; index_Voxel < vVoxels.size(); index_Voxel++)
		{// get canvas voxels and create material points
			MaterialPoint_BC *newMP;
			newMP = MP_Factory.createMaterialPoint(glm::dvec3(vVoxels[index_Voxel].d2_Position,0.0));

			newMP->i_ID = vVoxels[index_Voxel].u_ID;
			newMP->p_Material = pInconel;

			newMP->d_Volume_Initial = dOffset*dOffset*dOffset * vVoxels[index_Voxel].d_ESO_Opacity;
			newMP->d_Volume = newMP->d_Volume_Initial;

			newMP->d_Mass = newMP->p_Material->d_Density * newMP->d_Volume;

			newMP->d3_Position = glm::dvec3(vVoxels[index_Voxel].d2_Position,0.0) + glm::dvec3(0.0,0.0,0.5*d3Length_Grid.z);
			newMP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
			newMP->d3_Force_External = newMP->d_Mass * glm::dvec3(0.0, -00.0, 0.0);

			allMaterialPoint.push_back(newMP);
			// monitor
			newMP->b_Monitor = false;
			if(vVoxels[index_Voxel].b_ESO == true)
			{
				newMP->b_Mark_ESO = true;
				newMP->b_Monitor = true;
				//v_MarkedMaterialPoints_Stress_Monitor.push_back((MaterialPoint_Classic_CC *)newMP);
			}
			// moment log
			v_MarkedMaterialPoints_Momentum.push_back((MaterialPoint_Classic_CC *)newMP);
		}
	}

	double dPlatenSpeed = 0.1;
	if(true)
	{// top platen material points -------------------------------------------- platen MP
		std::vector<MaterialPoint_BC *> thisMaterialDomain = MP_Factory.createDomain_Cuboid(glm::dvec3(0.8*d3Length_Grid.x,0.5*d3Length_Grid.y,0.5*d3Length_Grid.z), glm::dvec3(0.002,0.002,d3Length_Grid.z), 0.5*d3Length_Grid.z);
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
//			if(thisMP->d3_Position.x > 0.9*d3Length_Grid.x)
			{
				thisMP->b_DisplacementControl = true;
				thisMP->f_DisplacementControl_Multiplier = -1.0;
				thisMP->d3_Velocity = glm::dvec3(0.0,0.0,0.0);
				v_MarkedMaterialPoints_Displacement_Control.push_back(thisMP);
				v_MarkedMaterialPoints_Displacement_Monitor.push_back(thisMP);
			}
		}
	}

	d_TimeIncrement_Maximum = 1.0/1.0*5.0e-8;
	d_TimeEnd = 2.0e-4;
	d_TimeConsole_Interval = 1e-5;

	// timeline events -------------------------------------------------------
//	m_TimeLine.addTimePoint(0.0,					glm::dvec3(0.0, 0.0, 0.0));
	m_TimeLine.addTimePoint(0.0,					glm::dvec3(0.0, +dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(1.0e-4,					glm::dvec3(0.0, +dPlatenSpeed, 0.0));
	m_TimeLine.addTimePoint(1.0e-4+1.0e-24,			glm::dvec3(0.0, 0.0, 0.0));
	m_TimeLine.addTimePoint(d_TimeEnd,				glm::dvec3(0.0, 0.0, 0.0));
//	m_TimeLine.addTimePoint(d_TimeEnd,				glm::dvec3(0.0, +dPlatenSpeed, 0.0));

	double dMass_Domain = 0.0;
	for(unsigned int index_MP = 0; index_MP < allMaterialPoint.size(); index_MP++)
	{// calculate debug values
		dMass_Domain += allMaterialPoint[index_MP]->d_Mass;
	}

	a_Runtime.fill(0.0);
	d_DampingCoefficient = 0.1;

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
		sDescription += "Material Point count: " + Script(allMaterialPoint.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Grid Resolution: (" + Script(i3Cells.x) + "," + Script(i3Cells.y) + "," + Script(i3Cells.z) + ")" + "(" + Script(d3Length_Cell.x,3) + ")\n";
		sDescription += "dOffset: " + Script(dOffset,4) + "\n";
//		sDescription += "Division (Angular): " + Script(iDivision_Angular) + " (offset: " + Script(_PI*dDiameter_Average/iDivision_Angular,4) + ")" + "\n";
//		sDescription += "Division (Radial): " + Script(iDivision_Radial) + " (offset: " + Script(dThickness_Ring/iDivision_Radial,4) + ")" + "\n";
//		sDescription += "Division (Longitudinal): " + Script(iDivision_Longitudinal) + " (offset: " + Script(dLength_Ring/iDivision_Longitudinal,4) + ")" + "\n";
//		sDescription += "Tube average_diameter: " + Script(dDiameter_Average,3) + "\n";
//		sDescription += "Tube thickness: " + Script(dThickness_Ring,3) + "\n";
//		sDescription += "Tube length: " + Script(dLength_Ring,3) + "\n";
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
