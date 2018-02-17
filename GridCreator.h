/* This class contains all we need to build the mesh */
#ifndef GRIDCREATOR_H
#define GRIDCREATOR_H

#include <string>
#include "Array_3D_Template.h"
#include "Node.h"
#include <vector>
#include <map>

class GridCreator{
	private:
		// Dictionnary with the materials and the chosen unsigned char assigned to it:
		map<string,unsigned char> materialNameForMaterialID;
	public:
		// Constructor:
		GridCreator(map<string,unsigned char> materialNameForMaterialID){
			this->materialNameForMaterialID = materialNameForMaterialID;
		}
		// Destructor:
		~GridCreator(void);
		// Assign to each node there coordinate and corresponding material:
		void nodesInitialization(Array_3D_Template<Node> &nodes);
};

#endif