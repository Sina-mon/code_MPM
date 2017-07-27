#ifndef MATERIALPOINT_FACTORY_BC_H
#define MATERIALPOINT_FACTORY_BC_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "Definitions.h"

#include "Definitions.h"
#include "MaterialPoint_BC.h"
#include "MaterialPoint_Classic_CC.h"
#include "MaterialPoint_CPDI_CC.h"

class MaterialPoint_Factory_BC
{
	public:
		MaterialPoint_Factory_BC(){;}
		virtual ~MaterialPoint_Factory_BC(){;}

		bool isInside(glm::dvec3 d3Coordinate, std::vector<glm::dvec3> vVertex);
		bool isInside(glm::dvec3 d3Coordinate, glm::dvec3 d3Center, glm::dvec2 d2Radii); // is inside a ellipse

//		virtual MaterialPoint_BC *createMaterialPoint(glm::dvec3 d3Center, double dOffset) {std::cout << "MaterialPoint_BC::createMaterialPoint, method not implemented yet!" << std::endl;}

		virtual std::vector<MaterialPoint_BC *> createDomain_Cuboid(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_Cuboid, method not implemented yet!" << std::endl;}
		virtual std::vector<MaterialPoint_BC *> createDomain_Polygon(std::vector<glm::dvec3> vVertex, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_Polygon, method not implemented yet!" << std::endl;}
		virtual std::vector<MaterialPoint_BC *> createDomain_AuxeticCell_Polygon(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dDent, double dThickness, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_AuxeticCell_Polygon, method not implemented yet!" << std::endl;}
		virtual std::vector<MaterialPoint_BC *> createDomain_AuxeticCell_Swisscheese(glm::dvec3 d3Center, glm::dvec3 d3Dimension, glm::dvec2 d2Spacing, glm::dvec2 d2Radii, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_AuxeticCell_Swisscheese, method not implemented yet!" << std::endl;}
		virtual std::vector<MaterialPoint_BC *> createDomain_AuxeticCell(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_AuxeticCell, method not implemented yet!" << std::endl;}
		virtual std::vector<MaterialPoint_BC *> createDomain_AuxeticMesh(glm::dvec3 d3Origin, glm::ivec2 i2Array, glm::dvec3 d3Dimension, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_AuxeticMesh, method not implemented yet!" << std::endl;}
		virtual std::vector<MaterialPoint_BC *> createDomain_Sphere(glm::dvec3 d3Center, double dRadius_Outer, double dRadius_Inner, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_Sphere, method not implemented yet!" << std::endl;}
		virtual std::vector<MaterialPoint_BC *> createDomain_Tube(glm::dvec3 d3Center, glm::dvec3 d3Rotation, double dRadius_Outer, double dRadius_Inner, double dLength, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_Tube, method not implemented yet!" << std::endl;}
		virtual std::vector<MaterialPoint_BC *> createDomain_QuarterTube(glm::dvec3 d3Center, glm::dvec3 d3Rotation, double dRadius_Outer, double dRadius_Inner, double dLength, double dOffset) {std::cout << "MaterialPoint_BC::createDomain_QuarterTube, method not implemented yet!" << std::endl;}

		virtual std::string getScript(MaterialPoint_BC *thisMaterialPoint) {std::cout << "MaterialPoint_BC::getScript, method not implemented yet!" << std::endl;}
	protected:
	private:
};

#endif
