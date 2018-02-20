#include "GridCreator.h"

/* Grid initialization */
void GridCreator::meshInitialization(){
	if(this->lengthX <= 0.0 || this->lengthY <= 0.0 || this->lengthZ <= 0.0){
		cout << "One of the lengths of the domain is still <= to 0 (Line ";
		printf("%d - FILE %s).Aborting.\n",__LINE__,__FILE__);
		abort();
	}
	// First, compute the number of nodes from delta(X,Y,Z) and length(X,Y,Z):
	this->numberOfNodesInEachDir[0] = (size_t) this->lengthX / this->deltaX +1;
	this->numberOfNodesInEachDir[1] = (size_t) this->lengthY / this->deltaY +1;
	this->numberOfNodesInEachDir[2] = (size_t) this->lengthZ / this->deltaZ +1;
	this->totalNumberOfNodes        = this->numberOfNodesInEachDir[0] *
										this->numberOfNodesInEachDir[1] * 
										this->numberOfNodesInEachDir[2];
	// Initialize nodesElec:
	cout << "Allocate nodesElec\n";
	this->nodesElec.set_size_data(this->numberOfNodesInEachDir[0],
								  this->numberOfNodesInEachDir[1],
								  this->numberOfNodesInEachDir[2]);
	// Initialize nodesMagn. It has two nodes more, in each direction:
	cout << "Allocate nodesMagn\n";
	this->nodesMagn.set_size_data(this->numberOfNodesInEachDir[0]+(size_t)1,
								  this->numberOfNodesInEachDir[1]+(size_t)1,
								  this->numberOfNodesInEachDir[2]+(size_t)1);
	cout << "Done." << endl;
	// Now, go through each node and decide in which material they are: TODO
	// Missing also: initializtion of the heat mesh.
	cout << "GridCreator::meshInitialization::OUT\n";
}
		
/* Destructor */
GridCreator::~GridCreator(void){
	cout << "Destructor of Grid creator ok.\n";
}

void GridCreator::test(double deltaX, double deltaY, double deltaZ, 
					   double lengthX, double lengthY, double lengthZ, 
					   double deltaT){
	this->deltaX = deltaX;
	this->deltaY = deltaY;
	this->deltaZ = deltaZ;

	this->lengthX = lengthX;
	this->lengthY = lengthY;
	this->lengthZ = lengthZ;
	
	this->deltaT = deltaT;
}
