#include "GridCreator.h"

/* Grid initialization */
void GridCreator::meshInitialization(){
	// First, compute the number of nodes from delta(X,Y,Z) and length(X,Y,Z):
	this->numberOfNodesInEachDir[0] = (size_t) this->lengthX / this->deltaX;
	this->numberOfNodesInEachDir[1] = (size_t) this->lengthY / this->deltaY;
	this->numberOfNodesInEachDir[2] = (size_t) this->lengthZ / this->deltaZ;
	this->totalNumberOfNodes        = this->numberOfNodesInEachDir[0] *
										this->numberOfNodesInEachDir[1] * 
										this->numberOfNodesInEachDir[2];
	// Initialize nodesElec:
	this->nodesElec.set_size_data(this->numberOfNodesInEachDir[0],
								  this->numberOfNodesInEachDir[1],
								  this->numberOfNodesInEachDir[2]);
	// Initialize nodesMagn. It has two nodes more, in each direction:
	this->nodesMagn.set_size_data(this->numberOfNodesInEachDir[0]+2,
								  this->numberOfNodesInEachDir[1]+2,
								  this->numberOfNodesInEachDir[2]+2);
	// Now, go through each node and decide in which material they are: TODO
	// Missing also: initializtion of the heat mesh.
	cout << "GridCreator::meshInitialization::OUT\n";
}
		
/* Destructor */
GridCreator::~GridCreator(void){
	cout << "Destructor of Grid creator ok.\n";
}