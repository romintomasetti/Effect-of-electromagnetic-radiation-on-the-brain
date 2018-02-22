/* This class contains all we need to build the mesh */
#ifndef GRIDCREATOR_H
#define GRIDCREATOR_H

#include <string>
#include <vector>
#include <map>

#include "Array_3D_Template.h"
#include "Node3DField.h"
#include "Node1DField.h"
#include "InputParser.h"
#include "Materials.h"


class GridCreator{
	public:
		// Dictionnary with the materials and the chosen unsigned char assigned to it:
		map<string,unsigned char> materialNameForMaterialID;
		// 3D array of nodes for electric field (those of the MPI process and of the neighboors):
		Array_3D_Template<Node3DField> nodesElec;
		// 3D array of nodes for magnetic field (those of the MPI process and of the neighboors):
		Array_3D_Template<Node3DField> nodesMagn;
		// 3D array of nodes for temperature    (those of the MPI process and of the neighboors):
		Array_3D_Template<Node1DField> nodesTemp;
		// DeltaX, DeltaY, DeltaZ:
		double deltaX = -1;
		double deltaY = -1;
		double deltaZ = -1;
		// Length of the domain along each direction:
		double lengthX = -1;
		double lengthY = -1;
		double lengthZ = -1;
		// Number of nodes along each direction (excluding neighboors' nodes):
		unsigned long nbrElts_X = 0;
		unsigned long nbrElts_Y = 0;
		unsigned long nbrElts_Z = 0;
		// Time increment:
		double deltaT = 0.0;
		// Number of nodes:
		size_t numberOfNodesInEachDir[3] = {0,0,0};
		size_t totalNumberOfNodes        = 0;

		InputParser &input_parser;
		Materials   &materials;
        ElectromagneticSource &elec_source;
//	public:
		// Constructor:
		GridCreator(InputParser &input_parser, Materials &materials):
		input_parser(input_parser), materials(materials){}
		// Destructor:
		~GridCreator(void);
	
		// Grid initialization:
		void meshInitialization();

		void test(double, double, double, double, double, double, double);


};

#endif
