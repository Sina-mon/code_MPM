#include <iostream>
#include <vector>

#include <omp.h>

#define STB_IMAGE_IMPLEMENTATION

#include "Canvas2D_CC.h"
#include "PhysicsEngine.h"
#include "GraphicsEngine.h"
#include "Definitions.h"

#include "ConstitutiveRelation.h"

bool	isLarger(Voxel_ST *V1, Voxel_ST *V2)
{
	return(V1->d_Objective > V2->d_Objective);
}
bool	isSmaller(Voxel_ST *V1, Voxel_ST *V2)
{
	return(V1->d_Objective < V2->d_Objective);
}

double	getObjective(std::vector<MaterialPoint_BC *> vMP)
{
	ConstitutiveRelation CR;

	double dObjective = 0.0;
	for(int index = 0; index < vMP.size(); index++)
	{
		MaterialPoint_BC *pMP = vMP[index];

		dObjective = pMP->d_Energy_Strain;// / pMP->d_Volume;

//		dObjective = CR.getState_J2(pMP->d6_Stress);

//		glm::dvec3 d3Principal = glm::abs(CR.getPrincipal(pMP->d6_Stress));
//		dObjective += glm::max(glm::max(d3Principal.x,d3Principal.y),d3Principal.z);
	}

	dObjective /= vMP.size();

	//float fObjective = thisMP->d_Energy_Strain;
	//float fJ2 = CR.getState_J2(thisMP->d6_Stress);
	//float fJ2 = glm::abs(thisMP->d6_Stress[0]);

	return dObjective;
}
/*
int main (int argc, char ** argv)
{// without graphics

	std::cout << "size: " << sizeof(MaterialPoint_BC) << std::endl;

	PhysicsEngine thePhysicsEngine;

//	thePhysicsEngine.initializeWorld_Classic_Cellular_Langrand();
	thePhysicsEngine.initializeWorld_Classic_Cellular_Shim_Square();

	double dTimeRequest = thePhysicsEngine.getTime_End();
	thePhysicsEngine.runSimulation_Classic_DoublePass_MPLocks(dTimeRequest);

	std::cout << "Execution finished, have a nice day!" << std::endl;

	return(0);
}
*/
/*
int main (int argc, char ** argv)
{ // with graphics
	//single simulation
	// physics engine initialization ------------------------------------------
	PhysicsEngine thePhysicsEngine;

//	thePhysicsEngine.initializeWorld_Classic_Bar();
//	thePhysicsEngine.initializeWorld_Classic_Cellular_Langrand();
	thePhysicsEngine.initializeWorld_Classic_Cellular_Shim_Square();

//	thePhysicsEngine.initializeWorld_Classic_Cellular_Graded();

	// graphics engine initialization -----------------------------------------
	GraphicsEngine theGraphicsEngine;
	// run simulation
	theGraphicsEngine.runVisualization(&thePhysicsEngine, false);

	std::cout << "Execution finished, have a nice day!" << std::endl;

	return(0);
}
*/

int main (int argc, char ** argv)
{
	// Create canvas, cantilever problem
//	double dOffset = 0.0010;
//	glm::dvec2 d2Size = glm::dvec2(0.100,0.100);
//	Canvas2D_CC Canvas2D(d2Size, dOffset);
//
//	Canvas2D.drawRectangle	(glm::dvec2(0.040,0.5*d2Size.y), glm::dvec2(0.080, 0.050), 0.0);
//	Canvas2D.setESORectangle(glm::dvec2(0.040,0.5*d2Size.y), glm::dvec2(0.080, 0.050), 0.0, true);
//	// loading location
//	Canvas2D.drawRectangle		(glm::dvec2(0.080+0.0*dOffset,0.5*d2Size.y), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0);
//	Canvas2D.setESORectangle	(glm::dvec2(0.080+0.0*dOffset,0.5*d2Size.y), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0, false);
//	Canvas2D.setLoadRectangle	(glm::dvec2(0.080+0.0*dOffset,0.5*d2Size.y), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0, true);

	// Create canvas, beam problem
	double dOffset = 0.00025;
	double dOpacity_Min = 1.0e-4;
	double dOpacity_Max = 1.0;
	glm::dvec2 d2Size = glm::dvec2(0.150,0.050);
	Canvas2D_CC Canvas2D(d2Size, dOffset);

//	Canvas2D.drawRectangle	(glm::dvec2(0.060,0.5*d2Size.y-0.5*0.040 + 0.5*0.012), glm::dvec2(0.120, 0.012), 0.0);
//	Canvas2D.setESORectangle(glm::dvec2(0.060,0.5*d2Size.y-0.5*0.040 + 0.5*0.012), glm::dvec2(0.120, 0.012), 0.0, true);
	Canvas2D.drawRectangle	(glm::dvec2(0.060,0.5*d2Size.y), glm::dvec2(0.120, 0.040), 0.0);
	Canvas2D.setESORectangle(glm::dvec2(0.060,0.5*d2Size.y), glm::dvec2(0.120, 0.040), 0.0, true);
	// loading location
	Canvas2D.drawRectangle		(glm::dvec2(0.0,0.5*d2Size.y - 0.020), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0);
	Canvas2D.setESORectangle	(glm::dvec2(0.0,0.5*d2Size.y - 0.020), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0, false);
	Canvas2D.setLoadRectangle	(glm::dvec2(0.0,0.5*d2Size.y - 0.020), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0, true);
	// Support location
	Canvas2D.drawRectangle		(glm::dvec2(0.120,0.5*d2Size.y - 0.020), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0);
	Canvas2D.setESORectangle	(glm::dvec2(0.120,0.5*d2Size.y - 0.020), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0, false);
	Canvas2D.setSupportRectangle(glm::dvec2(0.120,0.5*d2Size.y - 0.020), glm::dvec2(1.1*dOffset, 1.1*dOffset), 0.0, true);

	std::vector<Voxel_ST *> vVoxels = Canvas2D.getVoxels_Active(true);
	std::cout << "Canvas voxels: " << Script(Canvas2D.v_Voxels.size()) << " (" << Script(Canvas2D.u2_Size.x) << "," << Script(Canvas2D.u2_Size.y) << ")" << std::endl;
	std::cout << "Active voxels: " << Script(vVoxels.size()) << std::endl;

	for(unsigned int index = 0; index < vVoxels.size(); index++)
	{// start from lower portion
		double dTopEdge = 0.5*d2Size.y + 0.5*0.040 - 0.7*0.040;

		if(vVoxels[index]->d2_Position.y > dTopEdge)
		{
			vVoxels[index]->d_ESO_Opacity = dOpacity_Min;
		}
		else
			vVoxels[index]->d_ESO_Opacity = dOpacity_Max;
	}

	// ESO iterations
	double dRatio_Redundant = 0.1; // ratio in terms of active voxels
	double dRatio_Rejection = 0.7;
	double dRatio_Rejection_Increment = 0.0;
	double dRatio_Rejection_Volume_Upper = 0.1;
	double dRatio_Rejection_Volume_Lower = 0.0;
	double dRatio_Inclusion = 1.0 - dRatio_Rejection;
	double dInfluence_Radius = sqrt(2.0)*2.1*dOffset;
	int iActiveMPs = Canvas2D.getVoxels_Active(true).size();
	int iIteration_Local = 0;
	unsigned int iCounter_Removed = 0;
	unsigned int iCounter_Restored = 0;

	std::string sFilename_Log = "";
	std::string sFilename_Snapshot = "";
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time (&rawtime);
		timeinfo = localtime(&rawtime);

		// time format for filename
		strftime(buffer,80,"%Y%m%d%H%M%S",timeinfo);
		std::string strTime_Filename(buffer);

		sFilename_Log = strTime_Filename;
	}

	for(int iIteration_Global = 0; iIteration_Global < 10000; iIteration_Global++)
	{
		std::string sDescription = "";
		{
			sFilename_Snapshot = sFilename_Log + "_" + Script(iIteration_Global);

			sDescription += "######################################################################\n";
			//sDescription += "BESO process started on: \t" + strTime_Description + "\n";
			sDescription += "BESO iIteration_Global: \t" + Script(iIteration_Global) + "\n";
			sDescription += "BESO iIteration_Local: \t" + Script(iIteration_Local) + "\n";
			sDescription += "BESO dRatio_Rejection: \t" + Script(dRatio_Rejection, 3) + "\n";

			sDescription += "Removed in previous step: \t" + Script(iCounter_Removed) + "\n";
			sDescription += "Restored in previous step: \t" + Script(iCounter_Restored) + "\n";
			sDescription += "######################################################################\n";
		}
		// initialize physics engine
		PhysicsEngine thePhysicsEngine;
		thePhysicsEngine.initializeWorld_Classic_ESO(&Canvas2D, sFilename_Log, sFilename_Snapshot, sDescription);
		// initialize graphics engine
		GraphicsEngine theGraphicsEngine;
		// run simulation
		theGraphicsEngine.runVisualization(&thePhysicsEngine, true);

		// move objective values from simulation to canvas
		std::vector<MaterialPoint_BC *> vMaterialPoint = thePhysicsEngine.getMaterialPoints();
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

//			if(thisMP->b_Mark_ESO != true)
//				continue;

			std::vector<MaterialPoint_BC *> vMP;
			vMP.push_back(thisMP);
			double dObjective = getObjective(vMP);

			for(unsigned long int index_Voxel = 0; index_Voxel < Canvas2D.v_Voxels.size(); index_Voxel++)
			{
				if(thisMP->i_ID == Canvas2D.v_Voxels[index_Voxel]->u_ID)
				{
					Canvas2D.v_Voxels[index_Voxel]->d_Objective = 0.5*Canvas2D.v_Voxels[index_Voxel]->d_Objective + 0.5*dObjective;
				}
			}
		}
		// filter-smooth objective values
		Canvas2D.filterObjective_Smooth(dInfluence_Radius);// 2.1 seems to be better than 6.1

		std::vector<Voxel_ST *> vVoxels_Active = Canvas2D.getVoxels_Active(true);
		std::vector<Voxel_ST *> vVoxels_ESO = Canvas2D.getVoxels_ESO(true);
		// sort active voxels based on value, small to large
		{
			std::sort(vVoxels_ESO.begin(), vVoxels_ESO.end(), isSmaller);
			std::sort(vVoxels_Active.begin(), vVoxels_Active.end(), isSmaller);
		}
		// remove redundant voxels
		iCounter_Removed = 0;
		unsigned int iVoxels_Active = vVoxels_Active.size();
		unsigned int iVoxels_ESO = vVoxels_ESO.size();
		for(unsigned long int index_Voxel = 0; index_Voxel < dRatio_Rejection*iVoxels_ESO; index_Voxel++)
		{
//			if((double)iCounter_Removed > dRatio_Redundant*iVoxels_Active)
//				break;

			if(vVoxels_ESO[index_Voxel]->b_ESO == false)
				continue;

			if(vVoxels_ESO[index_Voxel]->d_ESO_Opacity > 0.8)
			{
				vVoxels_ESO[index_Voxel]->d_ESO_Opacity = dOpacity_Min;
				iCounter_Removed++;
			}
		}
//		std::cout << "iCounter_Removed: " << iCounter_Removed << std::endl;
		// restore essential voxels
		iCounter_Restored = 0;
		double dRatio_Inclusion = 1.0 - dRatio_Rejection;
		for(unsigned long int index_Voxel = dRatio_Rejection*iVoxels_ESO; index_Voxel < iVoxels_ESO; index_Voxel++)
		{
			if(vVoxels_ESO[index_Voxel]->b_ESO == false)
				continue;

			if(vVoxels_ESO[index_Voxel]->d_ESO_Opacity < 0.8)
			{
				vVoxels_ESO[index_Voxel]->d_ESO_Opacity = dOpacity_Max;
				iCounter_Restored++;
			}
		}
//		std::cout << "iCounter_Restored: " << iCounter_Restored << std::endl;
		if(iCounter_Removed < 0.001*iVoxels_Active || iIteration_Local == 200)
//		if(iCounter_Removed < 2)
		{
			dRatio_Rejection += dRatio_Rejection_Increment;
			iIteration_Local = 0;
		}
		else
			iIteration_Local++;

		//iActiveMPs = Canvas2D.getVoxels(true).size();
	}

	std::cout << "Execution finished, have a nice day!" << std::endl;

	return(0);
}
