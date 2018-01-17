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

		glm::dvec3 d3Principal = glm::abs(CR.getPrincipal(pMP->d6_Stress));
		dObjective += glm::max(glm::max(d3Principal.x,d3Principal.y),d3Principal.z);
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
	double dOffset = 0.00025;
	glm::dvec2 d2Size = glm::dvec2(0.02,0.04);
	Canvas2D_CC Canvas2D(d2Size, dOffset);

	Canvas2D.drawRectangle(glm::dvec2(0.4*d2Size.x+0.001,0.5*d2Size.y), glm::dvec2(0.8*d2Size.x, 0.5*d2Size.y), 0.0);
	Canvas2D.setESORectangle(glm::dvec2(0.4*d2Size.x+0.001,0.5*d2Size.y), glm::dvec2(0.8*d2Size.x, 0.5*d2Size.y), 0.0, true);
	Canvas2D.setESORectangle(glm::dvec2(0.8*d2Size.x,0.5*d2Size.y), glm::dvec2(0.005, 0.005), 0.0, false);
	Canvas2D.setESORectangle(glm::dvec2(0.0*d2Size.x,0.5*d2Size.y), glm::dvec2(0.005, 0.5*d2Size.y), 0.0, false);

	std::vector<Voxel_ST> vVoxels = Canvas2D.getVoxels(true);
	std::cout << "Canvas voxels: " << Script(Canvas2D.v_Voxels.size()) << " (" << Script(Canvas2D.u2_Size.x) << "," << Script(Canvas2D.u2_Size.y) << ")" << std::endl;
	std::cout << "Active voxels: " << Script(vVoxels.size()) << std::endl;

	// ESO iterations
	double dObjective_Threshold = 0.02;
	double dObjective_Threshold_Increment = 0.02;
	double dInfluence_Radius = 10.0*dOffset;
	for(int iIteration = 0; iIteration < 100; iIteration++)
	{
		// initialize physics engine
		PhysicsEngine thePhysicsEngine;
		thePhysicsEngine.initializeWorld_Classic_ESO(&Canvas2D);
		// initialize graphics engine
		GraphicsEngine theGraphicsEngine;
		// run simulation
		theGraphicsEngine.runVisualization(&thePhysicsEngine, true);

		// find maximum objective variable
		std::vector<MaterialPoint_BC *> vMaterialPoint = thePhysicsEngine.getMaterialPoints();
		double dObjective_Maximum = 1.0e-12;
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
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			if(thisMP->b_Mark_ESO == false)
				continue;

			std::vector<MaterialPoint_BC *> vMP;

			// find influencing MPs
			for(int index_MP2 = 0; index_MP2 < vMaterialPoint.size(); index_MP2++)
			{
				MaterialPoint_BC *thisMP2 = vMaterialPoint[index_MP2];

				if(thisMP2->b_Mark_ESO == false)
					continue;

				double dDistance = glm::length(thisMP2->d3_Position - thisMP->d3_Position);

				if(dDistance < dInfluence_Radius)
					vMP.push_back(thisMP2);
			}

			double dObjective = getObjective(vMP);

			if(dObjective < dObjective_Threshold * dObjective_Maximum)
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
		}

		dObjective_Threshold += dObjective_Threshold_Increment;
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

