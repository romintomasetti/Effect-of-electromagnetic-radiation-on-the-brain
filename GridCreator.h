/* This class contains all we need to build the mesh */
#ifndef GRIDCREATOR_H
#define GRIDCREATOR_H

#include <string>
#include <vector>
#include <map>

#include "Array_3D_Template.h"
#include "Node3DField.h"
#include "Node1DField.h"


class GridCreator{
	private:
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
		// Number of nodes along each direction (excluding neighboors' nodes):
		unsigned long nbrElts_X = 0;
		unsigned long nbrElts_Y = 0;
		unsigned long nbrElts_Z = 0;
		// Time increment:
		double deltaT = 0.0;
	public:
		// Constructor:
		GridCreator(map<string,unsigned char> materialNameForMaterialID){
			this->materialNameForMaterialID = materialNameForMaterialID;
		}
		// Destructor:
		~GridCreator(void);
	
		// Grid initialization:
		void meshInitialization(){};
};

#endif