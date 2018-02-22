#include "GridCreator.h"

/* Grid initialization */
void GridCreator::meshInitialization(){

	if(this->lengthX <= 0.0 || this->lengthY <= 0.0 || this->lengthZ <= 0.0){
		cout << "One of the lengths of the domain is still <= to 0 (Line ";
		printf("%d - FILE %s).Aborting.\n",__LINE__,__FILE__);
		std::abort();
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

	// Fonction à faire: obtenir les indices de début selon x,y,z de cette manière this->indices_x_first this->indices_y_first this->indices_z_first



	//Fonction obtenir mu et epsilon en fonction des indices données (globaux)  this->epsilon    this->mu

}

/*
 * Constructor of the class GridCreator.
 * Arguments are:
 * 		1)
 * 		2)
 * 		3)
 */
GridCreator::GridCreator(InputParser &input_parser,
						 Materials &materials,
						 MPI_Initializer &MPI_communicator):
						 input_parser(input_parser),
						 materials(materials),
						 MPI_communicator(MPI_communicator){
	// The DELTAS are common to any subgrid, so we can take that from the
	// object reference 'input_parser'. However, the lengths along each direction
	// cannot be taken from this object because this object gives the lengths
	// of the whole domain and we need the lengths of this particular subgrid !!
	this->deltaX = this->input_parser.deltaX;
	this->deltaY = this->input_parser.deltaY;
	this->deltaZ = this->input_parser.deltaZ;

	// To get the lengths of this subgrid, call MpiDivision:
	this->MPI_communicator.MpiDivision(*this);
}
		
/* Destructor */
GridCreator::~GridCreator(void){
	cout << "Destructor of Grid creator ok.\n";
}


 void GridCreator::LocalToGlobal(unsigned long *localIndices, unsigned long *globalIndices){
 	int i=0;
	//std::vector<int> globalIndices;

	for(i=0; i<3; i++){
		globalIndices[i] = localIndices[i] + this->originIndices[i];
		cout << globalIndices[i] << endl;
	}
 }

