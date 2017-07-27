#ifndef MATERIALPOINT_FACTORY_CPDI_CC_H
#define MATERIALPOINT_FACTORY_CPDI_CC_H

#include "MaterialPoint_Factory_BC.h"

class MaterialPoint_Factory_CPDI_CC : public MaterialPoint_Factory_BC
{
	public:
		MaterialPoint_Factory_CPDI_CC() {;}
		virtual ~MaterialPoint_Factory_CPDI_CC() {;}

		std::vector<MaterialPoint_BC *> createMaterialPoint(glm::dvec3 d3Center, double dOffset);

		std::vector<MaterialPoint_BC *> createCell_Hexahedron(glm::dvec3 d3Center, double dOffset);
		std::vector<MaterialPoint_BC *> createCell_Hexahedron(std::vector<glm::dvec3> vCorner);

		virtual std::vector<MaterialPoint_BC *> createDomain_Cuboid(glm::dvec3 d3Center, glm::dvec3 d3Dimension, double dOffset);
		virtual std::vector<MaterialPoint_BC *> createDomain_Tube(glm::dvec3 d3Center, glm::dvec3 d3Rotation, double dRadius_Outer, double dRadius_Inner, double dLength, double dOffset);
		virtual std::vector<MaterialPoint_BC *> createDomain_Tube_Smooth(glm::dvec3 d3Center, glm::dvec3 d3Rotation, double dRadius_Outer, double dRadius_Inner, double dLength, double dOffset);

		double getVolume(MaterialPoint_CPDI_CC *pMP)
		{
			double dVolume = 0.0;

			glm::dvec4 x = glm::dvec4(pMP->d3_Corner[0].x, pMP->d3_Corner[1].x, pMP->d3_Corner[2].x, pMP->d3_Corner[3].x);
			glm::dvec4 y = glm::dvec4(pMP->d3_Corner[0].y, pMP->d3_Corner[1].y, pMP->d3_Corner[2].y, pMP->d3_Corner[3].y);
			glm::dvec4 z = glm::dvec4(pMP->d3_Corner[0].z, pMP->d3_Corner[1].z, pMP->d3_Corner[2].z, pMP->d3_Corner[3].z);

			double x12 = (x[0] - x[1]); double x13 = (x[0] - x[2]); double x14 = (x[0] - x[3]);
			double x21 = (x[1] - x[0]); double x23 = (x[1] - x[2]); double x24 = (x[1] - x[3]);
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

			dVolume += x21*(y23*z34 - y34*z23);
			dVolume += x32*(y34*z12 - y12*z34);
			dVolume += x43*(y12*z23 - y23*z12);
			dVolume /= 6.0;

			return(dVolume);
		}

//		virtual std::string getScript(MaterialPoint_BC *thisMaterialPoint);
	protected:

	private:
};

#endif // MATERIALPOINT_FACTORY_CLASSIC_CC_H
