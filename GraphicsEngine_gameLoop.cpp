#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::gameLoop(bool bExit)
{
	int iSimulationStatus = 0;
	while(e_GameState != GameState::EXIT)
	{
		// removed for runtime analysis
		drawGame();

		// save snapshot
//		if(false)	//for runtime analysis
		if(f_Time - f_TimeSnapshot_LastSave > f_TimeSnapshot_Interval)
		{
			drawGame();
//			drawGame();

			f_TimeSnapshot_LastSave = f_Time;
//			std::string strFileName = _STR_SNAPFILE;
			std::string strFileName = mpm_PhysicsEngine->str_Snapshot_FileName;
			strFileName += Script((int)(f_Time*1000000),6);
//			strFileName += Script(i_TimeSnapshotCycle);
			strFileName += ".bmp";
			this->saveScreenshot(0, 0, i_ScreenWidth, i_ScreenHeight, strFileName.c_str());
			i_TimeSnapshotCycle += 1;
			std::cout << "GraphicsEngine: Screenshot saved at: " << Script(f_Time, 4) << std::endl;
		}

//		double dTimeIncrement_Request = 500*mpm_PhysicsEngine->getTime_Increment();
		double dTimeIncrement_Request = 0.5*mpm_PhysicsEngine->getTime_ConsoleInterval();
//		double dTimeIncrement_Request = mpm_PhysicsEngine->getTime_End();
		if(iSimulationStatus == 0)
		{
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_CPDI_SinglePass_MP(dTimeIncrement_Request);
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_CPDI_SinglePass_MPLocks(dTimeIncrement_Request);
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_CPDI_DoublePass_MPLocks(dTimeIncrement_Request);
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_CPDI_MultiBody_SinglePass_MPLocks(dTimeIncrement_Request);

//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_DoublePass_MPLocks_Contact(dTimeIncrement_Request);
			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_DoublePass_MPLocks(dTimeIncrement_Request);
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_DoublePass_MP(dTimeIncrement_Request);

//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_CPDI_SinglePass(dTimeIncrement_Request);

//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_SinglePass_MP(dTimeIncrement_Request);
//			iSimulationStatus = mpm_PhysicsEngine->runSimulation_Classic_SinglePass_MP_Contact(dTimeIncrement_Request);
		}
		if(bExit == true && iSimulationStatus == 1)
		{
			e_GameState = GameState::EXIT;
		}

		//processInput();
		this->f_Time = mpm_PhysicsEngine->getTime_Current();
	}
}
// ----------------------------------------------------------------------------
