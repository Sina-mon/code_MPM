#ifndef GRIDPOINT_MEDIATOR_H
#define GRIDPOINT_MEDIATOR_H

#include <iostream>
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

		void findAdjacentGridPoints(glm::dvec3 d3Position)
		{
			v_adjacentGridPoints.resize(0);

			glm::ivec3 i3BottomLeftFar = d3Position / d3_Length_Cell;

			int index = 0;

			for(int ix = 0; ix < 2; ix++)
			{
				for(int iy = 0; iy < 2; iy++)
				{
					for(int iz = 0; iz < 2; iz++)
					{
						glm::ivec3 i3Node_Index = i3BottomLeftFar + glm::ivec3(ix, iy, iz);

						if(i3Node_Index.x < 0)	continue;
						if(i3Node_Index.y < 0)	continue;
						if(i3Node_Index.z < 0)	continue;

						if(i3Node_Index.x > i3_Node_Count.x-1)	continue;
						if(i3Node_Index.y > i3_Node_Count.y-1)	continue;
						if(i3Node_Index.z > i3_Node_Count.z-1)	continue;

						index = GridPoint_Factory::getIndex(i3Node_Index, i3_Node_Count);
						v_adjacentGridPoints.push_back(index);
					}
				}
			}
			// remove duplicates
//			std::sort( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() );
//			v_adjacentGridPoints.erase( std::unique( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() ), v_adjacentGridPoints.end() );
		}
		void findAdjacentGridPoints_CPDI(MaterialPoint_CPDI_CC *pMP)
		{// fina AGP for the corners of a 3D tetrahedra
			v_adjacentGridPoints.resize(0);

			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				glm::dvec3 d3Position = pMP->d3_Corner[index_Corner];

				glm::ivec3 i3BottomLeftFar = d3Position / d3_Length_Cell;

				for(int ix = 0; ix < 2; ix++)
				{
					for(int iy = 0; iy < 2; iy++)
					{
						for(int iz = 0; iz < 2; iz++)
						{
							glm::ivec3 i3Node_Index = i3BottomLeftFar + glm::ivec3(ix, iy, iz);

							if(i3Node_Index.x < 0)	continue;
							if(i3Node_Index.y < 0)	continue;
							if(i3Node_Index.z < 0)	continue;

							if(i3Node_Index.x > i3_Node_Count.x-1)	continue;
							if(i3Node_Index.y > i3_Node_Count.y-1)	continue;
							if(i3Node_Index.z > i3_Node_Count.z-1)	continue;

							int index = GridPoint_Factory::getIndex(i3Node_Index, i3_Node_Count);
							v_adjacentGridPoints.push_back(index);
						}
					}
				}
			}
			// remove duplicates
			std::sort( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() );
			v_adjacentGridPoints.erase( std::unique( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() ), v_adjacentGridPoints.end() );
		}

		void calculateBases_Classic(glm::dvec3 d3Position_MP, glm::dvec3 d3Position_GP)
		{
			glm::dvec3 d3Distance = d3Position_MP - d3Position_GP;

			glm::dvec3 d3ShapeValue = glm::dvec3(0.0, 0.0, 0.0);
			d3ShapeValue.x = 1.0 - fabs(d3Distance.x) / d3_Length_Cell.x;
			d3ShapeValue.y = 1.0 - fabs(d3Distance.y) / d3_Length_Cell.y;
			d3ShapeValue.z = 1.0 - fabs(d3Distance.z) / d3_Length_Cell.z;

			if(d3ShapeValue.x < 0.0 || d3ShapeValue.x > 1.0)
				d3ShapeValue.x = 0.0;

			if(d3ShapeValue.y < 0.0 || d3ShapeValue.y > 1.0)
				d3ShapeValue.y = 0.0;

			if(d3ShapeValue.z < 0.0 || d3ShapeValue.z > 1.0)
				d3ShapeValue.z = 0.0;

			// results to pointers
			this->d_ShapeValue = d3ShapeValue.x * d3ShapeValue.y * d3ShapeValue.z;

			this->d3_ShapeGradient.x = -d3ShapeValue.y * d3ShapeValue.z * glm::sign(d3Distance.x) / d3_Length_Cell.x;//sina, make sure these are formulated correctly
			this->d3_ShapeGradient.y = -d3ShapeValue.x * d3ShapeValue.z * glm::sign(d3Distance.y) / d3_Length_Cell.y;
			this->d3_ShapeGradient.z = -d3ShapeValue.x * d3ShapeValue.y * glm::sign(d3Distance.z) / d3_Length_Cell.z;
		}
		void calculateBases_CPDI(MaterialPoint_CPDI_CC *pMP, glm::dvec3 d3Position_GP)
		{
			// shape value ----------------------------------------------------
			double dWeight[4] = {0.25, 0.25, 0.25, 0.25};
			double dShapeValue[4] = {0.0, 0.0, 0.0, 0.0};

			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				glm::dvec3 d3Position_MP = pMP->d3_Corner[index_Corner];

				glm::dvec3 d3Distance = d3Position_MP - d3Position_GP;

				glm::dvec3 d3ShapeValue = glm::dvec3(0.0, 0.0, 0.0);
				d3ShapeValue.x = 1.0 - fabs(d3Distance.x) / d3_Length_Cell.x;
				d3ShapeValue.y = 1.0 - fabs(d3Distance.y) / d3_Length_Cell.y;
				d3ShapeValue.z = 1.0 - fabs(d3Distance.z) / d3_Length_Cell.z;

				if(d3ShapeValue.x < 0.0 || d3ShapeValue.x > 1.0)
					d3ShapeValue.x = 0.0;

				if(d3ShapeValue.y < 0.0 || d3ShapeValue.y > 1.0)
					d3ShapeValue.y = 0.0;

				if(d3ShapeValue.z < 0.0 || d3ShapeValue.z > 1.0)
					d3ShapeValue.z = 0.0;

				dShapeValue[index_Corner] = d3ShapeValue.x * d3ShapeValue.y * d3ShapeValue.z;
			}

			this->d_ShapeValue = 0.0;
			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
				this->d_ShapeValue += dWeight[index_Corner] * dShapeValue[index_Corner];

			// shape gradient -------------------------------------------------

			glm::dvec4 x = glm::dvec4(pMP->d3_Corner[0].x, pMP->d3_Corner[1].x, pMP->d3_Corner[2].x, pMP->d3_Corner[3].x);
			glm::dvec4 y = glm::dvec4(pMP->d3_Corner[0].y, pMP->d3_Corner[1].y, pMP->d3_Corner[2].y, pMP->d3_Corner[3].y);
			glm::dvec4 z = glm::dvec4(pMP->d3_Corner[0].z, pMP->d3_Corner[1].z, pMP->d3_Corner[2].z, pMP->d3_Corner[3].z);

			double x21 = (x[1] - x[0]); double x23 = (x[1] - x[2]); double x24 = (x[1] - x[3]);
			double x12 = (x[0] - x[1]); double x13 = (x[0] - x[2]); double x14 = (x[0] - x[3]);
			double x31 = (x[2] - x[0]); double x32 = (x[2] - x[1]); double x34 = (x[2] - x[3]);
			double x41 = (x[3] - x[0]); double x42 = (x[3] - x[1]); double x43 = (x[3] - x[2]);

			double y12 = (y[0] - y[1]); double y13 = (y[0] - y[2]); double y14 = (y[0] - y[3]);
			double y21 = (y[1] - y[0]); double y23 = (y[1] - y[2]); double y24 = (y[1] - y[3]);
			double y31 = (y[2] - y[0]); double y32 = (y[2] - y[1]); double y34 = (y[2] - y[3]);
			double y41 = (y[3] - y[0]); double y42 = (y[3] - y[1]); double y43 = (y[3] - y[2]);

			double z12 = (z[0] - z[1]); double z13 = (z[0] - z[2]); double z14 = (z[0] - z[3]);
			double z21 = (z[1] - z[0]); double z23 = (z[1] - z[2]); double z24 = (z[1] - z[3]);
			double z31 = (z[2] - z[0]); double z32 = (z[2] - z[1]); double z34 = (z[2] - z[3]);
			double z41 = (z[3] - z[0]); double z42 = (z[3] - z[1]); double z43 = (z[3] - z[2]);

			double dVolume = 0.0;
			dVolume += x21*(y23*z34 - y34*z23);
			dVolume += x32*(y34*z12 - y12*z34);
			dVolume += x43*(y12*z23 - y23*z12);
			dVolume /= 1.0;

			double da[4] = {0.0, 0.0, 0.0, 0.0};
			da[0] = y42*z32 - y32*z42;
			da[1] = y31*z43 - y34*z13;
			da[2] = y24*z14 - y14*z24;
			da[3] = y13*z21 - y12*z31;

			double db[4] = {0.0, 0.0, 0.0, 0.0};
			db[0] = x32*z42 - x42*z32;
			db[1] = x43*z31 - x13*z34;
			db[2] = x14*z24 - x24*z14;
			db[3] = x21*z13 - x31*z12;

			double dc[4] = {0.0, 0.0, 0.0, 0.0};
			dc[0] = x42*y32 - x32*y42;
			dc[1] = x31*y43 - x34*y13;
			dc[2] = x24*y14 - x14*y24;
			dc[3] = x13*y21 - x12*y31;

			glm::dvec3 d3Weight[4] = {glm::dvec3(0.0,0.0,0.0), glm::dvec3(0.0,0.0,0.0), glm::dvec3(0.0,0.0,0.0), glm::dvec3(0.0,0.0,0.0)};
			d3Weight[0].x = da[0]; d3Weight[1].x = da[1]; d3Weight[2].x = da[2]; d3Weight[3].x = da[3];
			d3Weight[0].y = db[0]; d3Weight[1].y = db[1]; d3Weight[2].y = db[2]; d3Weight[3].y = db[3];
			d3Weight[0].z = dc[0]; d3Weight[1].z = dc[1]; d3Weight[2].z = dc[2]; d3Weight[3].z = dc[3];

			this->d3_ShapeGradient = glm::dvec3(0.0, 0.0, 0.0);
			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				d3_ShapeGradient.x += dShapeValue[index_Corner] * d3Weight[index_Corner].x;
				d3_ShapeGradient.y += dShapeValue[index_Corner] * d3Weight[index_Corner].y;
				d3_ShapeGradient.z += dShapeValue[index_Corner] * d3Weight[index_Corner].z;
			}

			d3_ShapeGradient *= (1.0/dVolume);
		}
		void findNeighborGridPoints(glm::dvec3 d3Position, glm::ivec3 i3Node_Count, glm::dvec3 d3Length_Cell, int iLayers)
		{// similar to getAdjacentGridPoint (iLayer=0), but can return extra layers of neighboring grid points
			v_adjacentGridPoints.resize(0);

			glm::ivec3 i3BottomLeftFar = d3Position / d3Length_Cell;

			int index = 0;

			for(int ix = -iLayers; ix < 2+iLayers; ix++)
			{
				for(int iy = -iLayers; iy < 2+iLayers; iy++)
				{
					for(int iz = -iLayers; iz < 2+iLayers; iz++)
					{
						glm::ivec3 i3Node_Index = i3BottomLeftFar + glm::ivec3(ix, iy, iz);

						if(i3Node_Index.x < 0)	continue;
						if(i3Node_Index.y < 0)	continue;
						if(i3Node_Index.z < 0)	continue;

						if(i3Node_Index.x > i3Node_Count.x-1)	continue;
						if(i3Node_Index.y > i3Node_Count.y-1)	continue;
						if(i3Node_Index.z > i3Node_Count.z-1)	continue;

						index = GridPoint_Factory::getIndex(i3Node_Index, i3Node_Count);
						v_adjacentGridPoints.push_back(index);
					}
				}
			}
			// remove duplicates
//			std::sort( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() );
//			v_adjacentGridPoints.erase( std::unique( v_adjacentGridPoints.begin(), v_adjacentGridPoints.end() ), v_adjacentGridPoints.end() );
		}

	protected:
	private:
};

#endif // GRIDPOINT_MEDIATOR_H
