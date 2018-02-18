#include "GridCreator.h"
/* Initializes the coordinates and the material of all the nodes inside the mesh */
void GridCreator::nodesInitialization(Array_3D_Template<Node> &nodes){
	// Acquiring the size of 'nodes':
	vector<size_t> sizes = nodes.get_size_data();
	printf("Size of 'nodes' is (%ld,%ld,%ld).\n",sizes[0],sizes[1],sizes[2]);
	// What's next in this function ?
	// 1) Use this->materialNameForMaterialID to get the material ID from its name
	printf("ID of the material AIR   is %d.\n",materialNameForMaterialID["AIR"]);
	printf("ID of the material WATER is %d.\n",materialNameForMaterialID["WATER"]);
	// 2) Assign to each node their coordinates and the material, as well as 
	//    initial temperature or magn./elec. field if needed.
}
		
/* Destructor */
GridCreator::~GridCreator(void){
}