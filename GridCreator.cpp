#include "GridCreator.h"
#include <omp.h>

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
	/* ASSIGN TO EACH NODE A MATERIAL */
	this->assignToEachNodeAMaterial();
	cout << "GridCreator::meshInitialization::ASSIGN_TO_EACH_NODE_A_MATERIAL::DONE\n";

	// Missing also: initializtion of the heat mesh.


	cout << "GridCreator::meshInitialization::OUT\n";





}

void GridCreator::assignToEachNodeAMaterial(void){
	// Print the type of simulation it is:
	cout << "GridCreator::assignToEachNodeAMaterial : ";
	cout << this->input_parser.get_SimulationType() << endl;
	// Now try to determine which type of simulation it is.
	// Then assign to each node its material and initial temperature.
	if(this->input_parser.get_SimulationType() == "USE_AIR_EVERYWHERE"){
		cout << "\n\nTYPE OF SIMULATION IS AIR EVERYWHERE\n\n";
		cout << "\n\nAssigning material and initial temperature for each node.\n";
		
		#pragma omp parallel for collapse(3)
		for(size_t K = 0 ; K < this->numberOfNodesInEachDir[2] ; K ++){
			for(size_t J = 0 ; J < this->numberOfNodesInEachDir[1] ; J ++ ){
				for(size_t I = 0 ; I < this->numberOfNodesInEachDir[0] ; I ++){
					/* Initialize material */
					this->nodesElec(I,J,K).material    = 
								this->materials.materialID_FromMaterialName["AIR"];
					this->nodesMagn(I,J,K).material    = 
								this->materials.materialID_FromMaterialName["AIR"];

					/* Initialize temperature */
					this->nodesElec(I,J,K).Temperature = 
								this->input_parser.GetInitTemp_FromMaterialName["AIR"];
					this->nodesMagn(I,J,K).Temperature = 
								this->input_parser.GetInitTemp_FromMaterialName["AIR"];
				}
			}
		}
		printf("At node(5,5,5) we have material %d.\n",this->nodesElec(0,0,0).material);
		cout << "This corresponds to " + 
			this->materials.materialName_FromMaterialID[this->nodesElec(0,0,0).material];
		cout << endl;
		printf("Initial temperature at this node is %f.\n",
			this->nodesElec(0,0,0).Temperature);
	}
}

/*
 * Constructor of the class GridCreator.
 * Arguments are:
 * 		1) Object input_parser contains information from input file.
 * 		2) Object materials contains information on all used materials.
 * 		3) Object MPI_communicator used for MPI communications.
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

