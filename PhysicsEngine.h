#ifndef PHYSICSENGINE_H
#define PHYSICSENGINE_H

#include <new>
#include <math.h>
#include <vector>
#include <algorithm>
#include <time.h>

#include <omp.h>

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp"//sina, glm is a column major implementation
// -Wl,--stack,10000000,--heap=10000000
#include "Canvas2D_CC.h"
#include "GridPoint.h"
#include "GridPoint_Factory.h"
#include "GridPoint_Mediator.h"
#include "Material_BC.h"
#include "MaterialPoint_BC.h"
#include "MaterialPoint_Factory_Classic_CC.h"
#include "MaterialPoint_Factory_CPDI_CC.h"
#include "ConstitutiveRelation.h"
#include "TimeLine.h"

#define _MAX_N_THREADS	8
#define _MAX_N_BODIES	1

class PhysicsEngine
{
	public:
		PhysicsEngine() {;}
		virtual ~PhysicsEngine();

		void initializeWorld_Classic_ESO(Canvas2D_CC *pCanvas, std::string sFileName_Log = "", std::string sFileName_Snapshot = "", std::string sDescription = "");

		void initializeWorld_Classic_Bar(void);
		void initializeWorld_Classic_HalfRing_Xiang_PlainStress(void);
		void initializeWorld_Classic_Ring_Fan(void);

		void initializeWorld_Classic_Cellular_Langrand(void);
		void initializeWorld_Classic_Cellular_Langrand_Hexagonal(void);

		void initializeWorld_Classic_Cellular_Shim_Square(void);

		void initializeWorld_Classic_Cellular_Graded(void);

		void initializeWorld_Classic_Foam(void);

		glm::dvec3 d3_Length_World = glm::dvec3(0.0, 0.0, 0.0);
		// MPM ----------------------------------------------------------------
		GridPoint_Mediator mpm_GP_Mediator_Thread[_MAX_N_THREADS];
		int	runSimulation_Classic_DoublePass_MPLocks(double dTimeIncrement_Total);
		int	runSimulation_Classic_DoublePass_MPLocks_Contact(double dTimeIncrement_Total);

		// simulation components
		int runSimulation_ResetGrid(void);
		int runSimulation_FindAGPs(int iThread);
		int runSimulation_M2G(int nThread);
		int runSimulation_IntegrateGrid(double dTimeIncrement);
		int runSimulation_DisplacementControl(void);
		int runSimulation_G2P_P2_SmallStrain(double dTimeIncrement);
		int runSimulation_G2P_P2_LargeStrain(double dTimeIncrement);
		// function to communicate with outside -------------------------------
		double getTime_Runtime(void) {return(d_Runtime_Total);}
		double getTime_Current(void) {return(d_Time);}
		double getTime_End(void) {return(d_TimeEnd);}
		double getTime_Increment(void) {return(d_TimeIncrement_Maximum);}
		double getTime_ConsoleInterval(void) {return(d_TimeConsole_Interval);}
		double getTime_SnapshotInterval(void) {return(d_TimeSnapshot_Interval);}
		// graphics interface -------------------------------------------------
		unsigned int	getCount_MaterialPoint(void) {return(allMaterialPoint.size());}
		unsigned int	getCount_MaterialPoint_CPDI(void) {return(allMaterialPoint_CPDI.size());}
		unsigned int 	getCount_GridPoint(void) {return(allGridPoint.size());}
		std::vector<MaterialPoint_BC *>			getMaterialPoints(void) {return(allMaterialPoint);}
		std::vector<MaterialPoint_CPDI_CC *>	getMaterialPoints_CPDI(void) {return(allMaterialPoint_CPDI);}
		std::vector<GridPoint *>				getGridPoints(void) {return(allGridPoint);}

		std::string str_Snapshot_FileName = _STR_SNAPFILE;
	protected:
		TimeLine m_TimeLine;
		double d_Mass_Minimum = 0.0;
		double d_DampingCoefficient = 0.0;

		double d_Runtime_Total = 0.0;
		std::array<double, 8> a_Runtime;

		std::vector<GridPoint *> allGridPoint;
		std::vector<Material_BC *> v_allMaterial;
		std::vector<MaterialPoint_BC *> allMaterialPoint;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Displacement_Monitor;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Displacement_Control;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Principal_Monitor;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Stress_Monitor;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Force_Monitor;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Momentum;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Monitor_Energy;

		std::vector<MaterialPoint_CPDI_CC *> allMaterialPoint_CPDI;
		std::vector<MaterialPoint_CPDI_CC *> v_MarkedMaterialPoints_CPDI_Displacement_Control;

		// multi-body related parameters
		//std::vector<GridPoint *> allGridPoint_Body[_MAX_N_BODIES];

		// parallelization related members
		std::vector<omp_lock_t *> v_GridPoint_Lock;

		// timing
		double d_Time = 0.0; // current simulation time
		double d_TimeIncrement_Maximum = 1.0e-3;
		double d_TimeEnd = 10.0;//5.0e-4;
		int i_TimeCycle = 0;

		double d_TimeConsole_Interval = 500.0*d_TimeIncrement_Maximum;
		double d_TimeSnapshot_Interval = 1.0e24;
		double d_TimeConsole_Last = -1.0e12; // before creation

		std::string str_Log_FileName = _STR_LOGFILE;
		void reportConsole(std::string sDescription = "");
	private:
};

#endif // PHYSICSENGINE_H
