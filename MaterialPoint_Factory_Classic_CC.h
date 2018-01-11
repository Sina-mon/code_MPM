#ifndef MATERIALPOINT_FACTORY_CLASSIC_CC_H
#define MATERIALPOINT_FACTORY_CLASSIC_CC_H

#include "MaterialPoint_Factory_BC.h"

class MaterialPoint_Factory_Classic_CC : public MaterialPoint_Factory_BC
{
	public:
		MaterialPoint_Factory_Classic_CC() {;}
		virtual ~MaterialPoint_Factory_Classic_CC() {;}

		MaterialPoint_BC *createMaterialPoint(glm::dvec3 d3Center);

		virtual std::vector<MaterialPoint_BC *> createDomain_Cuboid(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset);
//		virtual std::vector<MaterialPoint_BC *> createDomain_Polygon(std::vector<glm::dvec3> vVertex, double dOffset);
//		virtual std::vector<MaterialPoint_BC *> createDomain_AuxeticCell_Polygon(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dDent, double dThickness, double dOffset);
//		virtual std::vector<MaterialPoint_BC *> createDomain_AuxeticCell_Swisscheese(glm::dvec3 d3Center, glm::dvec3 d3Dimension, glm::dvec2 d2Spacing, glm::dvec2 d2Radii, double dOffset);
//		virtual std::vector<MaterialPoint_BC *> createDomain_AuxeticCell(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset);
//		virtual std::vector<MaterialPoint_BC *> createDomain_AuxeticMesh(glm::dvec3 d3Origin, glm::ivec2 i2Array, glm::dvec3 d3Dimension, double dOffset);
//		virtual std::vector<MaterialPoint_BC *> createDomain_Tube(glm::dvec3 d3Center, glm::dvec3 d3Rotation, double dRadius_Outer, double dRadius_Inner, double dLength, double dOffset);

		virtual std::string getScript(MaterialPoint_BC *thisMaterialPoint);
	protected:

	private:
};

#endif // MATERIALPOINT_FACTORY_CLASSIC_CC_H
