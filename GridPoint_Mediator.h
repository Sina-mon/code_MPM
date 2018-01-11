#ifndef GRIDPOINT_MEDIATOR_H
#define GRIDPOINT_MEDIATOR_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

#include "Definitions.h"
#include "MaterialPoint_CPDI_CC.h"
#include "GridPoint.h"
#include "GridPoint_Factory.h"

class GridPoint_Mediator
{
	public:
		// for multi-threading ----------------------------
		std::vector<unsigned int> v_adjacentGridPoints;
		double d_ShapeValue = 0.0;
		glm::dvec3 d3_ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);

		double getShapeValue(void) {return(d_ShapeValue);}
		glm::dvec3 getShapeGradient(void) {return(d3_ShapeGradient);}
		// ------------------------------------------------
		glm::dvec3 d3_Length_Grid = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 d3_Length_Cell = glm::dvec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Cells = glm::ivec3(0.0, 0.0, 0.0);
		glm::ivec3 i3_Node_Count = glm::ivec3(0.0, 0.0, 0.0);

		GridPoint_Mediator() {;}
		virtual ~GridPoint_Mediator() {;}

		void findAdjacentGridPoints(glm::dvec3 d3Position);
		void findAdjacentGridPoints_CPDI(MaterialPoint_CPDI_CC *pMP);
		void calculateBases_Classic(glm::dvec3 d3Position_MP, glm::dvec3 d3Position_GP);
		void calculateBases_CPDI(MaterialPoint_CPDI_CC *pMP, glm::dvec3 d3Position_GP);
		void findNeighborGridPoints(glm::dvec3 d3Position, glm::ivec3 i3Node_Count, glm::dvec3 d3Length_Cell, int iLayers);

	protected:
	private:
};

#endif // GRIDPOINT_MEDIATOR_H
