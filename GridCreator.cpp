#include "GridCreator.h"
#include "Node3DField.h"
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
	this->numberOfNodesInEachDirTemp[0] = this->numberOfNodesInEachDir[0];

	this->numberOfNodesInEachDir[1] = (size_t) this->lengthY / this->deltaY +1;
	this->numberOfNodesInEachDirTemp[1] = this->numberOfNodesInEachDir[1];

	this->numberOfNodesInEachDir[2] = (size_t) this->lengthZ / this->deltaZ +1;
	this->numberOfNodesInEachDirTemp[2] = this->numberOfNodesInEachDir[2];

	this->totalNumberOfNodes        = this->numberOfNodesInEachDir[0] *
										this->numberOfNodesInEachDir[1] * 
										this->numberOfNodesInEachDir[2];
	this->totalNumberOfNodesTemp = this->numberOfNodesInEachDirTemp[0] *
									this->numberOfNodesInEachDirTemp[1] *
									this->numberOfNodesInEachDirTemp[2];
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
	// Initialize temperature field:
	cout << "Allocate nodesTemp\n";
	this->nodesTemp.set_size_data(this->numberOfNodesInEachDirTemp[0],
									this->numberOfNodesInEachDirTemp[1],
									this->numberOfNodesInEachDirTemp[2]);

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


//  Function sending the values of the information of the face "char Direction" 
 Node3DField **GridCreator:GetVecSend(char Direction){
	 int i,j,k;


	 /* Determine the face needed */
	 
	 /* The x-direction is represented by 0 */
	 /* The y-direction is represented by 1 */
	 /* The y-direction is represented by 2 */

	 /* The table Table_send does not contain the empty places at the extremities
	 	==> we begin the index counting in GridCreator at 1 (in order to neglect the neighbours) */
	 /* GridCreator has numberOfNodesInEachDir[i]+1 nodes in the i-direction */

	 /* For each face, we try to send a table of class Node3D */
	 
	 
	 
	 /* Direction Y */
	 if (Direction == "W"  || Direction == "E"){
          Node3DField Table_send[this->numberOfNodesInEachDir[0]][this->numberOfNodesInEachDir[2]];
		  
		  /* Face W */
		  if (Direction == "W"){
             for(i=1;  i <= this->numberOfNodesInEachDir[0];  i++){
				 for(k=1;  k <= this->numberOfNodesInEachDir[2];  k++){
					 Table_send[i-1][k-1]=this->nodesElec(i,0,k);
				 }
			 }
		  }
		  
		  /* Face E */
		  else(Direction= "E"){
			 for(i=1;  i <= this->numberOfNodesInEachDir[0];  i++){
				 for(k=1;  k <= this->numberOfNodesInEachDir[2];  k++){
					 Table_send[i-1][k-1]=this->nodesElec(i,(this->numberOfNodesInEachDir[1])+1,k);
				 }
			 }
		  }
	 }
	 
	 
	 /* Direction X */
	 else if(Direction == "N" || Direction == "S"){
          Node3DField Table_send[this->numberOfNodesInEachDir[1]][this->numberOfNodesInEachDir[2]];
		  
		  /* Face N */
		  if(Direction == "N"){
			  for(j=1;  j <= this->numberOfNodesInEachDir[1];  j++){
				 for(k=1;  k <=this->numberOfNodesInEachDir[2];  k++){
					 Table_send[j-1][k-1]=this->nodesElec((this->numberOfNodesInEachDir[0])+1,j,k);
				 }
			 }
		  }
		  
		  /* Face S */
		  else(Direction == "S"){
			  for(j=1;  j <= this->numberOfNodesInEachDir[1];  j++){
				 for(k=1;  k <= this->numberOfNodesInEachDir[2];  k++){
					 Table_send[j-1][k-1]=this->nodesElec(0,j,k);
				 }
			 }
		  }
	 }
	 
	 
	 /* Direction Z */
	 else if(Direction == "U" || Direction == "D"){
		 Node3DField Table_send[this->numberOfNodesInEachDir[0]][this->numberOfNodesInEachDir[1]];
		 
		 /* Face U */
		 if(Direction == "U"){
          	for(i=1;  i <= this->numberOfNodesInEachDir[0];  i++){
				 for(j=1;j <= this->numberOfNodesInEachDir[1];j++){
					 Table_send[i-1][j-1]=this->nodesElec(i,j,(this->numberOfNodesInEachDir[2])+1);
				 }
			 }
		 }
		 
		 /* Face D */
		 else(Direction == "D"){
			 for(i=1;  i <= this->numberOfNodesInEachDir[0];  i++){
				 for(j=1;  j <= this->numberOfNodesInEachDir[1];  j++){
					 Table_send[i-1][j-1]=this->nodesElec(i,j,0);
				 }
			 }
		 }

	 }
	 return Table_send;
 }



// Fonction receive the value of the information of the face "char"
void SetVecReceive(char Direction, Node3DField **Table_receive){
	int i,j,k;


	/* Direction X */
	if(Direction == "N"){
		for(j=0;  j< this->numberOfNodesInEachDir[1];  j++){
			for(k=0;  k< this->numberOfNodesInEachDir[2];  k++){
				this->nodesElec((this->numberOfNodesInEachDir[1])+1,j+1,k+1) = Table_receive(j,k);
			}
		}
	}
	if(Direction == "S"){
		for(j=0;  j < this->numberOfNodesInEachDir[1];  j++){
			for(k=0;  k < this->numberOfNodesInEachDir[2];  k++){
				this->nodesElec(0,j+1,k+1) = Table_receive(j,k)
			}
		}
	}

	/* Direction Y */
	if(Direction == "W"){
		for(i=0;  i < this->numberOfNodesInEachDir[0];  i++){
			 for(k=0;  k < this->numberOfNodesInEachDir[2];  k++){
				 this->nodesElec(i+1,0,k+1) = Table_receive(i,k);
			 }
		}		
	}
	if(Direction == "E"){
		for(i=0;  i < this->numberOfNodesInEachDir[0];  i++){
			 for(k=0;  k < this->numberOfNodesInEachDir[2];  k++){
				 this->nodesElec(i+1,(this->numberOfNodesInEachDir[1])+1,k+1) = Table_receive(i,k);
			 }
		}
	}

	/* Direction Z */
	if(Direction == "U"){
		for(i=0;  i < this->numberOfNodesInEachDir[0];  i++){
			for(j=0;  j < this->numberOfNodesInEachDir[1];  j++){
				this->nodesElec(i+1,j+1,(this->numberOfNodesInEachDir[2])+1) = Table_receive(i,j)
			}
		}
	}
	if(Direction == "D"){
		for(i=0;  i<this->numberOfNodesInEachDir[0];  i++){
			for(j=0;  j< this->numberOfNodesInEachDir[1];  j++){
				this->nodesElec(i+1,j+1,0) = Table_receive(i,j)
			}
		}
	}
	if(Direction != "U" && Direction != "D" && Direction != "W" && Direction != "E" && Direction != "N" && Direction != "S"){
		printf("No face specify\n",__FILE__,__LINE__);
	}
}
