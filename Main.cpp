#include <iostream>
#include <vector>

#include <omp.h>

#define STB_IMAGE_IMPLEMENTATION

#include "Canvas2D_CC.h"
#include "PhysicsEngine.h"
#include "GraphicsEngine.h"
#include "Definitions.h"

#include "ConstitutiveRelation.h"

double getObjective(std::vector<MaterialPoint_BC *> vMP)
{
	ConstitutiveRelation CR;

	double dObjective = 0.0;
	for(int index = 0; index < vMP.size(); index++)
	{
		MaterialPoint_BC *pMP = vMP[index];

		dObjective = CR.getState_J2(pMP->d6_Stress);

//		glm::dvec3 d3Principal = glm::abs(CR.getPrincipal(pMP->d6_Stress));
//		dObjective += glm::max(glm::max(d3Principal.x,d3Principal.y),d3Principal.z);
	}

	dObjective /= vMP.size();

	//float fObjective = thisMP->d_Energy_Strain;
	//float fJ2 = CR.getState_J2(thisMP->d6_Stress);
	//float fJ2 = glm::abs(thisMP->d6_Stress[0]);

	return dObjective;
}

int main (int argc, char ** argv)
{
/*	single simulation
	// physics engine initialization ------------------------------------------
	PhysicsEngine thePhysicsEngine;
	thePhysicsEngine.initializeWorld_Classic_Cellular_Langrand_Hexagonal();
	// graphics engine initialization -----------------------------------------
	GraphicsEngine theGraphicsEngine;
	// run simulation
	theGraphicsEngine.runVisualization(&thePhysicsEngine);
*/

	// Create canvas
	double dOffset = 0.0002;
	glm::dvec2 d2Size = glm::dvec2(0.020,0.020);
	Canvas2D_CC Canvas2D(d2Size, dOffset);

	Canvas2D.drawRectangle	(glm::dvec2(0.008,0.5*d2Size.y), glm::dvec2(0.016, 0.010), 0.0);
	Canvas2D.setESORectangle(glm::dvec2(0.008,0.5*d2Size.y), glm::dvec2(0.016, 0.010), 0.0, true);
	// loading location
	Canvas2D.setESORectangle	(glm::dvec2(0.016,0.5*d2Size.y), glm::dvec2(0.002, 0.002), 0.0, false);
	Canvas2D.setLoadRectangle	(glm::dvec2(0.016,0.5*d2Size.y), glm::dvec2(0.001, 0.001), 0.0, true);
	// boundary location
	Canvas2D.setESORectangle(glm::dvec2(0.0,0.25*d2Size.y), glm::dvec2(0.005, 0.002), 0.0, false);
	Canvas2D.setESORectangle(glm::dvec2(0.0,0.75*d2Size.y), glm::dvec2(0.005, 0.002), 0.0, false);

	std::vector<Voxel_ST> vVoxels = Canvas2D.getVoxels(true);
	std::cout << "Canvas voxels: " << Script(Canvas2D.v_Voxels.size()) << " (" << Script(Canvas2D.u2_Size.x) << "," << Script(Canvas2D.u2_Size.y) << ")" << std::endl;
	std::cout << "Active voxels: " << Script(vVoxels.size()) << std::endl;

	// ESO iterations
	double dRatio_Rejection = 0.5;
	double dRatio_Rejection_Increment = 0.02;
	double dRatio_Rejection_Volume_Upper = 10.1;
	double dRatio_Rejection_Volume_Lower = 0.0;
	double dInfluence_Radius = 0.0001;
	int iActiveMPs = Canvas2D.getVoxels(true).size();
	for(int iIteration = 0; iIteration < 500; iIteration++)
	{
		// initialize physics engine
		PhysicsEngine thePhysicsEngine;
		thePhysicsEngine.initializeWorld_Classic_ESO(&Canvas2D);
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
					Canvas2D.v_Voxels[index_Voxel]->d_Objective = dObjective;
				}
			}
		}
		// find maximum objective value
		double dObjective_Maximum = 0.0;
		unsigned int iVoxels_Active = Canvas2D.getCount_Active(true);
		unsigned int iVoxels_ESO = Canvas2D.getCount_ESO(true);
		for(unsigned long int index_Voxel = 0; index_Voxel < Canvas2D.v_Voxels.size(); index_Voxel++)
		{
			if(Canvas2D.v_Voxels[index_Voxel]->b_Active == false)
				continue;

//			if(Canvas2D.v_Voxels[index_Voxel]->b_ESO == false)
//				continue;

			dObjective_Maximum += Canvas2D.v_Voxels[index_Voxel]->d_Objective / iVoxels_Active;

//			if(Canvas2D.v_Voxels[index_Voxel]->d_Objective > dObjective_Maximum)
//				dObjective_Maximum = Canvas2D.v_Voxels[index_Voxel]->d_Objective;
		}
		// reset redundancies
		for(unsigned long int index_Voxel = 0; index_Voxel < Canvas2D.v_Voxels.size(); index_Voxel++)
		{
			if(Canvas2D.v_Voxels[index_Voxel]->b_Active == false)
				continue;

			Canvas2D.v_Voxels[index_Voxel]->b_Redundant == false;
		}
		// mark redundant voxels
		for(unsigned long int index_Voxel = 0; index_Voxel < Canvas2D.v_Voxels.size(); index_Voxel++)
		{
			if(Canvas2D.v_Voxels[index_Voxel]->b_Active == false)
				continue;

			if(Canvas2D.v_Voxels[index_Voxel]->b_ESO == false)
				continue;

			// find neighbouring voxels
			std::vector <unsigned long int> vNeighbours;// hold indices from the canvas vector
			for(unsigned long int index_Voxel2 = 0; index_Voxel2 < Canvas2D.v_Voxels.size(); index_Voxel2++)
			{
				// sina, this determines whether non-existing voxels are included in the averaging or not
				// sina, it seems that excluding them gives more reasonable results
				if(Canvas2D.v_Voxels[index_Voxel2]->b_Active == false)
					continue;

				double dDistance = glm::length(Canvas2D.v_Voxels[index_Voxel]->d2_Position - Canvas2D.v_Voxels[index_Voxel2]->d2_Position);

				if(dDistance < dInfluence_Radius)
					vNeighbours.push_back(index_Voxel2);
			}

			// calculate average objective for neighbouring voxels
			double dObjective = 0.0;
			for(unsigned long int index_Neighbour = 0; index_Neighbour < vNeighbours.size(); index_Neighbour++)
			{
				dObjective += Canvas2D.v_Voxels[vNeighbours[index_Neighbour]]->d_Objective / vNeighbours.size();
			}

			// mark only center if the average is redundant
			if(dObjective < dRatio_Rejection * dObjective_Maximum)
			{
				Canvas2D.v_Voxels[index_Voxel]->b_Redundant = true;
			}

			// mark all neighbours if the average is redundant
//			if(dObjective < dRatio_Rejection * dObjective_Maximum)
//			{
//				for(unsigned long int index_Neighbour = 0; index_Neighbour < vNeighbours.size(); index_Neighbour++)
//				{
//					Canvas2D.v_Voxels[vNeighbours[index_Neighbour]]->b_Redundant = true;
//				}
//			}
		}
		// remove redundant voxels
		unsigned int iCounter_Removed = 0;
//		unsigned int iVoxels_Active = Canvas2D.getCount_Active(true);
		for(unsigned long int index_Voxel = 0; index_Voxel < Canvas2D.v_Voxels.size(); index_Voxel++)
		{
//			if((double)iCounter_Removed/iVoxels_Active > dRatio_Rejection_Volume_Upper)
//				break;

			if(Canvas2D.v_Voxels[index_Voxel]->b_Active == false)
				continue;

			if(Canvas2D.v_Voxels[index_Voxel]->b_ESO == false)
				continue;

			if(Canvas2D.v_Voxels[index_Voxel]->b_Redundant == false)
				continue;

			Canvas2D.v_Voxels[index_Voxel]->d_Objective = 0.0;
			Canvas2D.v_Voxels[index_Voxel]->b_Active = false;
			Canvas2D.v_Voxels[index_Voxel]->b_ESO = false;
			Canvas2D.v_Voxels[index_Voxel]->b_Redundant = false;

			iCounter_Removed++;
		}
/*
		// remove redundant voxels
		unsigned long int iCounter_Removed = 0;
		unsigned long int iVoxels_Ative = Canvas2D.getVoxels(true).size();
		for(unsigned long int index_Voxel = 0; index_Voxel < Canvas2D.v_Voxels.size(); index_Voxel++)
		{
			if((double)iCounter_Removed/iVoxels_Ative > dRatio_Rejection_Volume)
				break;

			if(Canvas2D.v_Voxels[index_Voxel]->b_Active == false)
				continue;

			if(Canvas2D.v_Voxels[index_Voxel]->b_ESO == false)
				continue;

			// find neighbouring voxels
			std::vector <unsigned long int> vNeighbours;// hold indices from the canvas vector
			double dDistance_Influence = 2.0*dOffset;
			for(unsigned long int index_Voxel2 = 0; index_Voxel2 < Canvas2D.v_Voxels.size(); index_Voxel2++)
			{
				double dDistance = glm::length(Canvas2D.v_Voxels[index_Voxel]->d2_Position - Canvas2D.v_Voxels[index_Voxel2]->d2_Position);

				if(dDistance < dDistance_Influence)
					vNeighbours.push_back(index_Voxel2);
			}

			// calculate average objective for neighbouring voxels
			double dObjective = 0.0;
			for(unsigned long int index_Neighbour = 0; index_Neighbour < vNeighbours.size(); index_Neighbour++)
			{
				dObjective += Canvas2D.v_Voxels[vNeighbours[index_Neighbour]]->d_Objective / vNeighbours.size();
			}

			// remove all neighbours if the average is redundant
			if(dObjective < dRatio_Rejection * dObjective_Maximum)
			{
				for(unsigned long int index_Neighbour = 0; index_Neighbour < vNeighbours.size(); index_Neighbour++)
				{
					Canvas2D.v_Voxels[vNeighbours[index_Neighbour]]->d_Objective = 0.0;
					Canvas2D.v_Voxels[vNeighbours[index_Neighbour]]->b_Active = false;

					iCounter_Removed++;
				}
			}

//			if(Canvas2D.v_Voxels[index_Voxel]->d_Objective < dRatio_Rejection * dObjective_Maximum)
//			{
//				Canvas2D.v_Voxels[index_Voxel]->d_Objective = 0.0;
//				Canvas2D.v_Voxels[index_Voxel]->b_Active = false;
//			}
		}
*/
/*
		// find maximum objective value
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			if(thisMP->b_Mark_ESO != true)
				continue;

			std::vector<MaterialPoint_BC *> vMP;
			vMP.push_back(thisMP);

			double dObjective = getObjective(vMP);

			if(dObjective > dObjective_Maximum)
				dObjective_Maximum = dObjective;
		}
		// search for MP near the threshold
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// bypass MP that are not to be touched
			if(thisMP->b_Mark_ESO == false)
				continue;

			std::vector<MaterialPoint_BC *> vMP;

			// find influencing MPs
			for(int index_MP2 = 0; index_MP2 < vMaterialPoint.size(); index_MP2++)
			{
				MaterialPoint_BC *thisMP2 = vMaterialPoint[index_MP2];

				if(thisMP2->b_Mark_ESO == false)
					continue;

				double dDistance = glm::length(glm::dvec2(thisMP2->d3_Position) - glm::dvec2(thisMP->d3_Position));

				if(dDistance < dInfluence_Radius)
					vMP.push_back(thisMP2);
			}

			double dObjective = getObjective(vMP);

			if(dObjective < dObjective_Threshold_Lower * dObjective_Maximum)
			{
				for(int index_Voxel = 0; index_Voxel < Canvas2D.v_Voxels.size(); index_Voxel++)
				{
					if(thisMP->i_ID == Canvas2D.v_Voxels[index_Voxel]->u_ID && Canvas2D.v_Voxels[index_Voxel]->b_ESO == true)
					{
//						Canvas2D.v_Voxels[index_Voxel]->d_ESO_Opacity *= 0.5;
//
//						if(Canvas2D.v_Voxels[index_Voxel]->d_ESO_Opacity < glm::pow(0.5,3))
							Canvas2D.v_Voxels[index_Voxel]->b_Active = false;
					}
				}
			}
//			if(dObjective > dObjective_Threshold_Upper * dObjective_Maximum)
//			{
//				for(int index_Voxel = 0; index_Voxel < Canvas2D.v_Voxels.size(); index_Voxel++)
//				{
//					if(thisMP->i_ID == Canvas2D.v_Voxels[index_Voxel]->u_ID && Canvas2D.v_Voxels[index_Voxel]->b_ESO == true)
//					{
//						if(Canvas2D.v_Voxels[index_Voxel]->d_ESO_Opacity < 0.9)
//							Canvas2D.v_Voxels[index_Voxel]->d_ESO_Opacity *= 2.0;
//					}
//				}
//			}
		}
*/
		if((double)iCounter_Removed/iVoxels_Active <= dRatio_Rejection_Volume_Lower)
		{
			dRatio_Rejection += dRatio_Rejection_Increment;
		}
//		if(Canvas2D.getVoxels(true).size() == iActiveMPs)
//		{
//			dRatio_Rejection += dRatio_Rejection_Increment;
//		}
		std::cout << "iIteration: " << Script(iIteration) << std::endl;
		std::cout << "dRatio_Rejection: " << Script(dRatio_Rejection, 3) << std::endl;
		std::cout << "######################################################################" << std::endl;
		iActiveMPs = Canvas2D.getVoxels(true).size();
	}

	std::cout << "Execution finished, have a nice day!" << std::endl;

//		std::vector<Voxel_ST> vVoxels2 = Canvas.getVoxels(true);
//		std::cout << "Canvas voxels: " << Script(Canvas.v_Voxels.size()) << " (" << Script(Canvas.i3_Size.x) << "," << Script(Canvas.i3_Size.y) << "," << Script(Canvas.i3_Size.z) << ")" << std::endl;
//		std::cout << "Active voxels: " << Script(vVoxels2.size()) << std::endl;
//
//	std::cout << "here" << std::endl;
//	{
//		PhysicsEngine thePhysicsEngine;
//		thePhysicsEngine.initializeWorld_Classic_ESO(&Canvas);
//
//		GraphicsEngine theGraphicsEngine;
//
//		theGraphicsEngine.runVisualization(&thePhysicsEngine);
//	}
//	std::cout << "here" << std::endl;

	return(0);
}

