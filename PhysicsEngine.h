#ifndef PHYSICSENGINE_H
#define PHYSICSENGINE_H

#include <math.h>
#include <vector>
#include<algorithm>
#include <time.h>

#include <omp.h>

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp"//sina, glm is a column major implementation

#include "GridPoint.h"
#include "GridPoint_Factory.h"
#include "GridPoint_Mediator.h"
#include "MaterialPoint_BC.h"
#include "MaterialPoint_Factory_Classic_CC.h"
#include "MaterialPoint_Factory_CPDI_CC.h"
#include "ConstitutiveRelation.h"
#include "TimeLine.h"

#define _MAX_N_THREADS 4

class PhysicsEngine
{
	public:
		PhysicsEngine() {};
		virtual ~PhysicsEngine();

		void	initializeWorld_Bar(void);
		void	initializeWorld_Bar_CPDI(void);
		void	initializeWorld_Ring(void);
		void	initializeWorld_CPDI_FullRing_Xiang_PlainStrain(void);
		void	initializeWorld_CPDI_HalfRing_Xiang_PlainStrain(void);
		void	initializeWorld_CPDI_HalfRing_Xiang_PlainStress_Modulus(void);
		void	initializeWorld_CPDI_HalfRing_Xiang_PlainStress_Runtime(void);
		void	initializeWorld_CPDI_HalfRing_Xiang_PlainStress(void);
		void	initializeWorld_CPDI_HalfRing_Xiang_FullLength(void);
		void 	initializeWorld_CPDI_HalfRing_Fan(void);
		void	initializeWorld_CPDI_HalfRing_Gupta_PlainStress(void);
		void	initializeWorld_CPDI_HalfRing_Xu_PlainStress(void);
		void	initializeWorld_CPDI_HalfRing_Shim_PlainStress_WaveSpeed(void);
		void	initializeWorld_Classic_HalfRing_Xiang_PlainStress(void);
		void	initializeWorld_Classic_Foam(void);
		void	initializeWorld_Classic_Foam_Bullet(void);
		void	initializeWorld_Classic_Foam_Ring(void);
		void	initializeWorld_Classic_Foam_HoneyComb(void);
		void	initializeWorld_QuarterRing_CPDI_Xiang(void);
		void	initializeWorld_Ring_CPDI_Xiang(void);
		void	initializeWorld_Ring_CPDI(void);
		void	initializeWorld_AuxeticSwisscheeseCell(void);
		void	initializeWorld_AuxeticPolygonCell(void);

		glm::dvec3 d3_Length_World = glm::dvec3(0.0, 0.0, 0.0);
		// MPM ----------------------------------------------------------------
		GridPoint_Mediator mpm_GP_Mediator_Thread[_MAX_N_THREADS];
		int		runSimulation_Classic_DoublePass_MP(double dTimeIncrement_Total);
		int		runSimulation_Classic_SinglePass_MP(double dTimeIncrement_Total);
		int		runSimulation_Classic_SinglePass_MP_Contact(double dTimeIncrement_Total);
		int		runSimulation_CPDI_SinglePass(double dTimeIncrement_Total);
		int		runSimulation_CPDI_SinglePass_MP(double dTimeIncrement_Total);
		int		runSimulation_CPDI_DoublePass_MP(double dTimeIncrement_Total);
		int		runSimulation_CPDI_SinglePass_MP_Locks(double dTimeIncrement_Total);

//		double d_Offset = 0.0;

		glm::dvec3 d3_Length_Grid_Kernel = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Length_Cell_Kernel = glm::dvec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Cells_Kernel = glm::ivec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Nodes_Kernel = glm::ivec3(0.0, 0.0, 0.0);

		// function to communicate with outside -------------------------------
		double getTime_Runtime(void) {return(d_Runtime_Total);}
		double getTime_Current(void) {return(d_Time);}
		double getTime_End(void) {return(d_TimeEnd);}
		double getTime_Increment(void) {return(d_TimeIncrement_Maximum);}
		double getTime_ConsoleInterval(void) {return(d_TimeConsole_Interval);}
		// graphics interface -------------------------------------------------
		unsigned int	getCount_MaterialPoint(void) {return(allMaterialPoint.size());}
		unsigned int	getCount_MaterialPoint_CPDI(void) {return(allMaterialPoint_CPDI.size());}
		unsigned int 	getCount_GridPoint(void) {return(allGridPoint.size());}
		std::vector<MaterialPoint_BC *>	getMaterialPoints(void) {return(allMaterialPoint);}
		std::vector<MaterialPoint_CPDI_CC *>	getMaterialPoints_CPDI(void) {return(allMaterialPoint_CPDI);}
		std::vector<GridPoint *>		getGridPoints(void) {return(allGridPoint);}
		std::vector<GridPoint *>		getGridPoints_Kernel(void) {return(v_GridPoint_Kernel);}
	protected:
		TimeLine m_TimeLine;
		double d_Mass_Minimum = 0.0;
		double d_DampingCoefficient = 0.0;

		double d_Runtime_Total = 0.0;
		std::array<double, 8> a_Runtime;

		std::vector<GridPoint *> allGridPoint;
		std::vector<GridPoint *> allGridPoint_Thread[_MAX_N_THREADS];
		std::vector<GridPoint *> v_GridPoint_Kernel;
		std::vector<MaterialPoint_BC *> allMaterialPoint;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Displacement_Monitor;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Displacement_Control;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Stress_Monitor;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Force_Monitor;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Momentum;
		std::vector<MaterialPoint_BC *> v_MarkedMaterialPoints_Monitor_Energy;

		std::vector<MaterialPoint_CPDI_CC *> allMaterialPoint_CPDI;
		std::vector<MaterialPoint_CPDI_CC *> v_MarkedMaterialPoints_CPDI_Displacement_Control;
		std::vector<MaterialPoint_CPDI_CC *> v_MarkedMaterialPoints_CPDI_Displacement_Monitor;

		std::vector<omp_lock_t *> v_GridPoint_Lock;

		double d_Time = 0.0; // current simulation time
		double d_TimeIncrement_Maximum = 1.0e-3;
		double d_TimeEnd = 10.0;//5.0e-4;
		int i_TimeCycle = 0;

		double d_TimeConsole_Interval = 500.0*d_TimeIncrement_Maximum;
		double d_TimeConsole_Last = -1.0e12; // before creation

		void reportConsole(std::string sDescription = "");
	private:
};

#endif // PHYSICSENGINE_H
