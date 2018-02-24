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
#include "MPI_Initializer.h"
#include "vtl/vtl.h"
#include "vtl/vtlVec3.h"
#include "vtl/vtlSPoints.h"

class MPI_Initializer;

class GridCreator{
	public:
		// Origin of the whole grid, that is, of all the simulation (same for all subgrids):
		vtl::Vec3d originOfWholeSimulation = vtl::Vec3d(0.0,0.0,0.0);

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
		// Indices of the domain's first point
		std::vector<unsigned long> originIndices;
		// Time increment:
		double deltaT = 0.0;
		// Number of nodes:
		size_t numberOfNodesInEachDir[3] = {0,0,0};
		size_t totalNumberOfNodes        = 0; 

		size_t numberOfNodesInEachDirTemp[3] = {0,0,0};
		size_t totalNumberOfNodesTemp       = 0;

		InputParser 	&input_parser;
		Materials   	&materials;
		MPI_Initializer &MPI_communicator;
        
//	public:
		// Constructor:
		GridCreator(InputParser &input_parser,
					Materials &materials,
					MPI_Initializer &MPI_communicator)/*:
					input_parser(input_parser),
					materials(materials),
					MPI_communicator(MPI_communicator){}*/;
		// Destructor:
		~GridCreator(void);
	
		// Grid initialization:
		void meshInitialization();
		
		// Convert the local indices to global indices:
		void LocalToGlobal(unsigned long *, unsigned long *);

		// Assign to each node a material. For the first tests, it will 
		// just set all materials to be air. For later simulations, we will
		// have to read input files to link each node to a material.
		void assignToEachNodeAMaterial(void);

};

#endif
