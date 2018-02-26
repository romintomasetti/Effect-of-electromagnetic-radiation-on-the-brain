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
#include "vtl.h"
#include "vtlVec3.h"
#include "vtlSPoints.h"

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

		/*
		 * This function fills in the array of Node3DField nodes
		 * with the nodes of the face specified by the 'char' argument.
		 * Also, sizes of the 2d array are given. If they are wrong, 
		 * they are corrected (so we need a reference to them).
		 */
		void GetVecSend(char, Node3DField ***,size_t &,size_t&);
		//Fonction receive the value of the information of the face "char"
		/*
		 * This function fills in the array of the grid with data
		 * received from the neighboors. Direction specified by
		 * the 'char' argument. Sizes are given for error chacking.
		 */ 
		void SetVecRecv(char, Node3DField ***,size_t,size_t);

		/* Communication with other MPI processes */
		void communicateWithNeighboor(
					char /*direction*/,
					size_t /*size1_send*/,
					size_t /*size2_send*/,
					size_t /*size1_recv*/,
					size_t /*size2_recv*/,
                    Node3DField *** /*ElectricNodes_toSend*/,
                    Node3DField *** /*ElectricNodes_toRecv*/,
                    MPI_Request ** /*requests_MPI*/);
};

#endif
