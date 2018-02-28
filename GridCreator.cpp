#include "GridCreator.h"
#include "Node3DField.h"
#include <omp.h>
#include <cmath>


/* Grid initialization */
void GridCreator::meshInitialization(){

	/*if(this->lengthX <= 0.0 || this->lengthY <= 0.0 || this->lengthZ <= 0.0){
		cout << "One of the lengths of the domain is still <= to 0 (Line ";
		printf("%d - FILE %s).Aborting.\n",__LINE__,__FILE__);
		std::abort();
	} */
	// First, compute the number of nodes from delta(X,Y,Z) and length(X,Y,Z):
	// The following is now computed directly inside MpiDivion (in mpi_communciator)
	//this->numberOfNodesInEachDir[0] = (size_t) (this->lengthX / this->deltaX) +1;
	this->numberOfNodesInEachDirTemp[0] = this->numberOfNodesInEachDir[0];

	//this->numberOfNodesInEachDir[1] = (size_t) (this->lengthY / this->deltaY) +1;
	this->numberOfNodesInEachDirTemp[1] = this->numberOfNodesInEachDir[1];

	//this->numberOfNodesInEachDir[2] = (size_t) (this->lengthZ / this->deltaZ) +1;
	this->numberOfNodesInEachDirTemp[2] = this->numberOfNodesInEachDir[2];

	this->totalNumberOfNodes        = this->numberOfNodesInEachDir[0] *
										this->numberOfNodesInEachDir[1] * 
										this->numberOfNodesInEachDir[2];
	this->totalNumberOfNodesTemp = this->numberOfNodesInEachDirTemp[0] *
									this->numberOfNodesInEachDirTemp[1] *
									this->numberOfNodesInEachDirTemp[2];

	printf("Size of electric mesh is (%ld,%ld,%ld).\n",this->numberOfNodesInEachDir[0],
					this->numberOfNodesInEachDir[1],this->numberOfNodesInEachDir[2]);
	
	// Initialize nodesElec:
	cout << "Allocate nodesElec\n";
	this->nodesElec.set_size_data(this->numberOfNodesInEachDir[0] + 2,
								  this->numberOfNodesInEachDir[1] + 2,
								  this->numberOfNodesInEachDir[2] + 2);
	// Initialize nodesMagn. It has one more node, in each direction:
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

	/* Initialize source */
	for(unsigned int I = 0 ; I < this->input_parser.source.get_number_of_sources() ; I ++){
		this->input_parser.source.computeNodesInsideSource(this->input_parser.lengthX,
													this->input_parser.lengthY,
													this->input_parser.lengthZ,
													this->input_parser.deltaX,
													this->input_parser.deltaY,
													this->input_parser.deltaZ,
													I);
	}

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
		int RANK = this->MPI_communicator.getRank();
		#pragma omp parallel for collapse(3)
		for(size_t K = 0 ; K < this->numberOfNodesInEachDir[2]+2 ; K ++){
			for(size_t J = 0 ; J < this->numberOfNodesInEachDir[1]+2 ; J ++ ){
				for(size_t I = 0 ; I < this->numberOfNodesInEachDir[0]+2 ; I ++){
					/* Initialize material */
					this->nodesElec(I,J,K).material    = 
								this->materials.materialID_FromMaterialName["AIR"];

					if(K < this->numberOfNodesInEachDir[2]+1 &&
						J < this->numberOfNodesInEachDir[1]+1 &&
						I < this->numberOfNodesInEachDir[0]+1){
						this->nodesMagn(I,J,K).material    = 
									this->materials.materialID_FromMaterialName["AIR"];
						this->nodesMagn(I,J,K).Temperature = 
								this->input_parser.GetInitTemp_FromMaterialName["AIR"];
						}

					/* Initialize temperature */
					this->nodesElec(I,J,K).Temperature = 
								this->input_parser.GetInitTemp_FromMaterialName["AIR"];
					
					if(K < this->numberOfNodesInEachDir[2] &&
						J < this->numberOfNodesInEachDir[1] &&
						I < this->numberOfNodesInEachDir[0])
						this->nodesTemp(I,J,K).field = RANK;
				}
			}
		}
		printf("At node(5,5,5) we have material %d.\n",this->nodesElec(0,0,0).material);
		cout << "This corresponds to " + 
			this->materials.materialName_FromMaterialID[this->nodesElec(0,0,0).material];
		cout << endl;
		printf("Initial temperature at this node is %f.\n",
			this->nodesElec(0,0,0).Temperature);
	}else if(this->input_parser.get_SimulationType() == "TEST_PARAVIEW"){
		cout << "Using paraview" << endl;
		
		for(size_t K = 0 ; K < this->numberOfNodesInEachDir[2]+1 ; K ++){
			for(size_t J = 0 ; J < this->numberOfNodesInEachDir[1]+1 ; J ++ ){
				for(size_t I = 0 ; I < this->numberOfNodesInEachDir[0]+1 ; I ++){
					unsigned long local[3];
					unsigned long global[3];
					local[0] = I;
					local[1] = J;
					local[2] = K;
					this->LocalToGlobal(local,global);
					this->nodesElec(I,J,K).field[0] = (double)global[0];
					this->nodesElec(I,J,K).field[1] = (double)global[1];
					this->nodesElec(I,J,K).field[2] = (double)global[2];
					this->nodesMagn(I,J,K).field[0] = (double)global[0];
					this->nodesMagn(I,J,K).field[1] = (double)global[1];
					this->nodesMagn(I,J,K).field[2] = (double)global[2];
					printf("%d::local(%ld,%ld,%ld) to global(%ld,%ld,%ld)\n",
						this->MPI_communicator.getRank(),
						local[0],local[1],local[2],global[0],global[1],global[2]);
					printf("nodesMagn(%ld,%ld,%ld) = %f\n",I,J,K,this->nodesMagn(I,J,K).field[0]);
				}
			}
		}
	}else if(this->input_parser.get_SimulationType() == "DEBUG_MPI_COMM"){
		int MPI_RANK = this->MPI_communicator.getRank()+1;
		for(size_t K = 0 ; K < this->numberOfNodesInEachDir[2]+1 ; K ++){
			for(size_t J = 0 ; J < this->numberOfNodesInEachDir[1]+1 ; J ++ ){
				for(size_t I = 0 ; I < this->numberOfNodesInEachDir[0]+1 ; I ++){
					this->nodesElec(I,J,K).field[0] = MPI_RANK;
					this->nodesElec(I,J,K).field[1] = MPI_RANK;
					this->nodesElec(I,J,K).field[2] = MPI_RANK;
				}
			}
		}
	}else{
		printf("GridCreator::assignToEachNodeAMaterial::ERROR\n");
		printf("\tWrong simulation type ! Has %s.\n",
			this->input_parser.get_SimulationType().c_str());
		std::abort();
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


 void GridCreator::LocalToGlobal(unsigned long *localIndices,
 						 unsigned long *globalIndices){
 	int i=0;

	for(i=0; i<3; i++){
		globalIndices[i] = localIndices[i] + this->originIndices[i];
	}
 }


//  Function sending the values of the information of the face "char Direction" 
void GridCreator::GetVecSend(char Direction,
				double **array_send,
				size_t &size1,size_t &size2){
	 size_t i,j,k;
	 size_t counter=0;
	 printf("GetVecSend IN : MPI %d OMP %d  direction %c %ld %ld \n", this->MPI_communicator.getRank(),
	 			omp_get_thread_num(), Direction, size1, size2 );

	 /* Determine the face needed */
	 
	 /* The x-direction is represented by 0 */
	 /* The y-direction is represented by 1 */
	 /* The y-direction is represented by 2 */

	 /* The table Table_send does not contain the empty places at the extremities
	 	==> we begin the index counting in GridCreator at 1 (in order to neglect the neighbours) */
	 /* GridCreator has numberOfNodesInEachDir[i]+1 nodes in the i-direction */

	 /* For each face, we try to send a table of class Node3D */
	 
	 
	 
	 /* Direction Y */
	 if (Direction == 'W'  || Direction == 'E'){
		// Check the sizes:
		if(size1 != this->numberOfNodesInEachDir[0]){
			printf("Aborting at line %d (direction %c, MPI %d), file %s\n",
				__LINE__,Direction,this->MPI_communicator.getRank(),__FILE__);
			abort();
		}
		if(size2 != this->numberOfNodesInEachDir[2]){
			printf("Aborting at line %d, file %s\n",__LINE__,__FILE__);
			abort();
		}
		  
		/* Face W */
		if (Direction == 'W'){
			int rank = this->MPI_communicator.getRank();
			for(i=1;  i <= this->numberOfNodesInEachDir[0];  i++){
				for(k=1;  k <= this->numberOfNodesInEachDir[2];  k++){
					(*array_send)[counter++]=this->nodesElec(i,1,k).field[0];
					(*array_send)[counter++]=this->nodesElec(i,1,k).field[1];
					(*array_send)[counter++]=this->nodesElec(i,1,k).field[2];
					/*printf("MPI %d :: send(%ld,%ld,%ld)=(%f,%f,%f)\n",
						rank,i,(size_t)1,k,
						this->nodesElec(i,1,k).field[0],
						this->nodesElec(i,1,k).field[1],
						this->nodesElec(i,1,k).field[2]);*/
				}
			}
		}
		  
		/* Face E */
		else/*(Direction == 'E')*/{
			int rank = this->MPI_communicator.getRank();
			for(i=1;  i <= this->numberOfNodesInEachDir[0];  i++){
				for(k=1;  k <= this->numberOfNodesInEachDir[2];  k++){
					(*array_send)[counter++]=this->nodesElec(i,(this->numberOfNodesInEachDir[1]),k).field[0];
					(*array_send)[counter++]=this->nodesElec(i,(this->numberOfNodesInEachDir[1]),k).field[1];
					(*array_send)[counter++]=this->nodesElec(i,(this->numberOfNodesInEachDir[1]),k).field[2];
					/*printf("MPI %d :: send(%ld,%ld,%ld)=(%f,%f,%f)\n",
						rank,i,(this->numberOfNodesInEachDir[1]),k,
						this->nodesElec(i,(this->numberOfNodesInEachDir[1]),k).field[0],
						this->nodesElec(i,(this->numberOfNodesInEachDir[1]),k).field[1],
						this->nodesElec(i,(this->numberOfNodesInEachDir[1]),k).field[2]);*/
				}
			}
		}

		if(counter != size1*size2*3){
			printf("Size1 = %ld, size2 = %ld, counter = %ld\n",size1,size2,counter);
			std::abort();
		}
	 }
	 
	 
	 /* Direction X */
	 else if(Direction == 'N' || Direction == 'S'){
		// Check the sizes:
		if(size1 != this->numberOfNodesInEachDir[1]){
			printf("Aborting at line %d (Direction %c, MPI %d), file %s\n",__LINE__,
				Direction,this->MPI_communicator.getRank(),__FILE__);
			printf("Expected %ld but has %ld\n",this->numberOfNodesInEachDir[1],
				size1);
			abort();
		}
		if(size2 != this->numberOfNodesInEachDir[2]){
			printf("Aborting at line %d, file %s\n",__LINE__,__FILE__);
			abort();
		} 
		/* Face N */
		if(Direction == 'N'){
			for(j=1;  j <= this->numberOfNodesInEachDir[1];  j++){
				for(k=1;  k <=this->numberOfNodesInEachDir[2];  k++){
					(*array_send)[counter++]=this->nodesElec((this->numberOfNodesInEachDir[0])+1,j,k).field[0];
					(*array_send)[counter++]=this->nodesElec((this->numberOfNodesInEachDir[0])+1,j,k).field[1];
					(*array_send)[counter++]=this->nodesElec((this->numberOfNodesInEachDir[0])+1,j,k).field[2];
				}
			}
		}
		
		/* Face S */
		else/*(Direction == 'S')*/{
			for(j=1;  j <= this->numberOfNodesInEachDir[1];  j++){
				for(k=1;  k <= this->numberOfNodesInEachDir[2];  k++){
					(*array_send)[counter++]=this->nodesElec(0,j,k).field[0];
					(*array_send)[counter++]=this->nodesElec(0,j,k).field[1];
					(*array_send)[counter++]=this->nodesElec(0,j,k).field[2];
				}
			}
		}

		if(counter != size1*size2*3){
			printf("Size1 = %ld, size2 = %ld, counter = %ld\n",size1,size2,counter);
			std::abort();
		}
	 }
	 
	 
	 /* Direction Z */
	 else if(Direction == 'U' || Direction == 'D'){
		// Check the sizes:
		if(size1 != this->numberOfNodesInEachDir[0]){
			printf("Aborting at line %d (Direction %c, MPI %d), file %s\n",__LINE__,Direction,
				this->MPI_communicator.getRank(),__FILE__);
			abort();
		}
		if(size2 != this->numberOfNodesInEachDir[1]){
			printf("Aborting at line %d, file %s\n",__LINE__,__FILE__);
			abort();
		}
		
		//Node3DField Table_send[this->numberOfNodesInEachDir[0]][this->numberOfNodesInEachDir[1]];
		 
		/* Face U */
		if(Direction == 'U'){
		for(i=1;  i <= this->numberOfNodesInEachDir[0];  i++){
				for(j=1;j <= this->numberOfNodesInEachDir[1];j++){
					(*array_send)[counter++]=this->nodesElec(i,j,(this->numberOfNodesInEachDir[2])+1).field[0];
					(*array_send)[counter++]=this->nodesElec(i,j,(this->numberOfNodesInEachDir[2])+1).field[1];
					(*array_send)[counter++]=this->nodesElec(i,j,(this->numberOfNodesInEachDir[2])+1).field[2];
				}
			}
		}
		
		
		/* Face D */
		else/*(Direction == 'D')*/{
			for(i=1;  i <= this->numberOfNodesInEachDir[0];  i++){
				for(j=1;  j <= this->numberOfNodesInEachDir[1];  j++){
					(*array_send)[counter++]=this->nodesElec(i,j,0).field[0];
					(*array_send)[counter++]=this->nodesElec(i,j,0).field[1];
					(*array_send)[counter++]=this->nodesElec(i,j,0).field[2];
				}
			}
		}

		if(counter != size1*size2*3){
			printf("Size1 = %ld, size2 = %ld, counter = %ld\n",size1,size2,counter);
			std::abort();
		}

	 }
	printf("GetVecSend OUT : MPI %d OMP %d  direction %c %ld %ld %ld %ld\n", this->MPI_communicator.getRank(),
	 			omp_get_thread_num(), Direction, size1, size2, size1*size2*3, counter);
 }



// Fonction receive the value of the information of the face "char"
void GridCreator::SetVecRecv(char Direction, double **Table_receive,size_t size1,size_t size2){
	printf("SetVecRecv IN : MPI %d OMP %d  direction %c %ld %ld \n", this->MPI_communicator.getRank(),
	 			omp_get_thread_num(), Direction, size1, size2 );
	size_t i,j,k;
	size_t count=0;;

	/* Direction X */
	if(Direction == 'N'){
		for(j=0;  j< this->numberOfNodesInEachDir[1];  j++){
			for(k=0;  k< this->numberOfNodesInEachDir[2];  k++){
				this->nodesElec((this->numberOfNodesInEachDir[1])+1,j+1,k+1).field[0] = (*Table_receive)[count++];
				this->nodesElec((this->numberOfNodesInEachDir[1])+1,j+1,k+1).field[1] = (*Table_receive)[count++];
				this->nodesElec((this->numberOfNodesInEachDir[1])+1,j+1,k+1).field[2] = (*Table_receive)[count++];
			}
		}
	}
	if(Direction == 'S'){
		for(j=0;  j < this->numberOfNodesInEachDir[1];  j++){
			for(k=0;  k < this->numberOfNodesInEachDir[2];  k++){
				this->nodesElec(0,j+1,k+1).field[0] = (*Table_receive)[count++];
				this->nodesElec(0,j+1,k+1).field[1] = (*Table_receive)[count++];
				this->nodesElec(0,j+1,k+1).field[2] = (*Table_receive)[count++];
			}
		}
	}

	

	/* Direction Y */
	if(Direction == 'W'){
		int rank = this->MPI_communicator.getRank();
		for(i=0;  i < this->numberOfNodesInEachDir[0];  i++){
			 for(k=0;  k < this->numberOfNodesInEachDir[2];  k++){
				 this->nodesElec(i+1,0,k+1).field[0] = (*Table_receive)[count++];
				 this->nodesElec(i+1,0,k+1).field[1] = (*Table_receive)[count++];
				 this->nodesElec(i+1,0,k+1).field[2] = (*Table_receive)[count++];
				/*printf("MPI %d :: recv(%ld,%ld,%ld)=(%f,%f,%f)\n",
					rank,i+1,(size_t)0,k+1,
					this->nodesElec(i+1,0,k+1).field[0],
					this->nodesElec(i+1,0,k+1).field[1],
					this->nodesElec(i+1,0,k+1).field[2]);*/
			 }
		}		
	}
	if(Direction == 'E'){
		int rank = this->MPI_communicator.getRank();
		for(i=0;  i < this->numberOfNodesInEachDir[0];  i++){
			 for(k=0;  k < this->numberOfNodesInEachDir[2];  k++){
				this->nodesElec(i+1,(this->numberOfNodesInEachDir[1])+1,k+1).field[0] = (*Table_receive)[count++];
				this->nodesElec(i+1,(this->numberOfNodesInEachDir[1])+1,k+1).field[1] = (*Table_receive)[count++];
				this->nodesElec(i+1,(this->numberOfNodesInEachDir[1])+1,k+1).field[2] = (*Table_receive)[count++];
				/*printf("MPI %d :: recv(%ld,%ld,%ld)=(%f,%f,%f)\n",
					rank,i+1,this->numberOfNodesInEachDir[1]+1,k+1,
					this->nodesElec(i+1,(this->numberOfNodesInEachDir[1])+1,k+1).field[0],
					this->nodesElec(i+1,(this->numberOfNodesInEachDir[1])+1,k+1).field[1],
					this->nodesElec(i+1,(this->numberOfNodesInEachDir[1])+1,k+1).field[2]);*/
			 }
		}
	}

	/* Direction Z */
	if(Direction == 'U'){
		for(i=0;  i < this->numberOfNodesInEachDir[0];  i++){
			for(j=0;  j < this->numberOfNodesInEachDir[1];  j++){

				this->nodesElec(i+1,j+1,(this->numberOfNodesInEachDir[2])+1).field[0] = (*Table_receive)[count++];
				this->nodesElec(i+1,j+1,(this->numberOfNodesInEachDir[2])+1).field[1] = (*Table_receive)[count++];
				this->nodesElec(i+1,j+1,(this->numberOfNodesInEachDir[2])+1).field[2] = (*Table_receive)[count++];
			}
		}
	}
	if(Direction == 'D'){
		for(i=0;  i<this->numberOfNodesInEachDir[0];  i++){
			for(j=0;  j< this->numberOfNodesInEachDir[1];  j++){
				this->nodesElec(i+1,j+1,0).field[0] = (*Table_receive)[count++];
				this->nodesElec(i+1,j+1,0).field[1] = (*Table_receive)[count++];
				this->nodesElec(i+1,j+1,0).field[2] = (*Table_receive)[count++];
			}
		}
	}
	if(Direction != 'U' && Direction != 'D' && Direction != 'W' && Direction != 'E' 
		&& Direction != 'N' && Direction != 'S'){
		printf("GridCreator::SetVecReceive::ERROR\n");
		printf("\tNo face specified, in %s at line %d.\n",__FILE__,__LINE__);
	}
	if(count != size1*size2*3){
		printf("Size1=%ld, Size2=%ld, count=%ld\n",size1,size2,count);
		abort();
	}
	printf("SetVecRecv OUT : MPI %d OMP %d  direction %c %ld %ld %ld %ld\n", this->MPI_communicator.getRank(),
	 			omp_get_thread_num(), Direction, size1, size2, size1*size2*3, count );
}


/////////////////////////
/// MPI COMMUNICATION ///
/////////////////////////
void GridCreator::communicateWithNeighboor(
					char direction,
					char type,
					size_t size1_send,
					size_t size2_send,
					size_t size1_recv,
					size_t size2_recv,
                    double **ElectricNodes_toSend,
                    double **ElectricNodes_toRecv,
                    MPI_Request **requests_MPI){
	/* Get OMP thread num */
	int omp_thread = omp_get_thread_num();
	if(omp_thread < 0 || omp_thread > 5){
		printf("GridCreator::communicateWithNeighboor::ERROR (MPI %d, OMP %d)\n",
			this->MPI_communicator.getRank(),omp_thread);
		std::abort();
	}
	int omp_Neighboor = -1;

	printf("> communicateWithNeighboor :: MPI %d :: OMP %d SIZE %ld \n",
		this->MPI_communicator.getRank(),
		omp_thread, size1_send*size2_send*3);

	int SIZE_VEC = 3;
	if(type == 'V'){
		// Vector, has 3 components:
		SIZE_VEC = 3;
	}else if(type == 'S'){
		// Scalar, has 1 component:
		SIZE_VEC = 1;
	}else{
		printf("GridCreator::communicateWithNeighboor::ERROR (line %d)\n",
			__LINE__);
		printf("\t'type' is either 'S' (scalar field) or 'V' (vector field).\n");
		printf("\tAborting.\n");
		std::abort();
	}
	    /////////////////////////////////////////////////
        /// CONVENTION FOR COMMUNICATION              ///
        /// 1) OMP thread(0) communicates with SOUTH. ///
        /// 2) OMP thread(1) communicates with NORTH. ///
        /// 3) OMP thread(2) communicates with WEST.  ///
        /// 3) OMP thread(3) communicates with EAST.  ///
        /// 4) OMP thread(4) communicates with DOWN.  ///
        /// 5) OMP thread(5) communicates with UP.    ///
        /////////////////////////////////////////////////
	int DECR = -1;
	if(omp_thread == 0){
		omp_Neighboor=1;
		DECR = 0;
		printf("Coucou from omp thread %d, from MPI %d.\n",omp_thread,
			this->MPI_communicator.getRank());
	}else if(omp_thread == 1){
		omp_Neighboor = 0;
		DECR = 2;
		printf("Coucou from omp thread %d, from MPI %d.\n",omp_thread,
			this->MPI_communicator.getRank());
	}else if(omp_thread == 2){
		omp_Neighboor =3;
		DECR = 4;
		printf("Coucou from omp thread %d, from MPI %d.\n",omp_thread,
			this->MPI_communicator.getRank());
	}else if(omp_thread == 3){
		omp_Neighboor = 2;
		DECR = 6;
		printf("Coucou from omp thread %d, from MPI %d.\n",omp_thread,
			this->MPI_communicator.getRank());
	}else if(omp_thread == 4){
		omp_Neighboor = 5;
		DECR = 8;
		printf("Coucou from omp thread %d, from MPI %d.\n",omp_thread,
			this->MPI_communicator.getRank());
	}else if(omp_thread == 5){
		omp_Neighboor = 4;
		DECR = 10;
		printf("Coucou from omp thread %d, from MPI %d.\n",omp_thread,
			this->MPI_communicator.getRank());
	}else{
		printf("GridCreator::communicateWithNeighboor::ERROR\n");
		printf("\tIn MPI %d, at line %d, in file %s\n",this->MPI_communicator.getRank(),
			__LINE__,__FILE__);
		std::abort();
	}


	/* NEIGHBOOR IS EVEN, I AM ODD */
	if((this->MPI_communicator.RankNeighbour[omp_thread]%2) == 0 &&
	   (this->MPI_communicator.getRank()%2) != 0){
		   	printf(">>>>>>>>>>>>>>>>>>>NeigEven :: OMP %d , RANK(%d)=%d, MYRANK=%d (DECR %d)\n",omp_thread,omp_thread,
				this->MPI_communicator.RankNeighbour[omp_thread],
				this->MPI_communicator.getRank(),DECR);

			/* We first send and then receive */
			/*MPI_Isend((*ElectricNodes_toSend),size1_send*size2_send*SIZE_VEC,
				MPI_DOUBLE,
				this->MPI_communicator.RankNeighbour[omp_thread],omp_thread,
				MPI_COMM_WORLD,&(*requests_MPI)[DECR+0]);*/
			MPI_Send((*ElectricNodes_toSend),size1_send*size2_send*SIZE_VEC,
				MPI_DOUBLE,
				this->MPI_communicator.RankNeighbour[omp_thread],omp_thread,
				MPI_COMM_WORLD);
			

			printf("MPI %d::MPI_Isend::done (request at %d).\n",
				this->MPI_communicator.getRank(),DECR);
			/*MPI_Irecv((*ElectricNodes_toRecv),size1_recv*size2_recv*SIZE_VEC,
				MPI_DOUBLE,this->MPI_communicator.RankNeighbour[omp_thread],
				omp_Neighboor,MPI_COMM_WORLD,&(*requests_MPI)[DECR+1]);*/
			MPI_Recv((*ElectricNodes_toRecv),size1_recv*size2_recv*SIZE_VEC,
				MPI_DOUBLE,this->MPI_communicator.RankNeighbour[omp_thread],
				omp_Neighboor,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			printf("MPI %d::MPI_Irecv::done (request at %d).\n",
				this->MPI_communicator.getRank(),DECR+1);

	/* NEIGHBOOR IS ODD, I AM EVEN */
	}else if( this->MPI_communicator.RankNeighbour[omp_thread]%2 != 0 &&
		this->MPI_communicator.getRank()%2 == 0){
			printf(">>>>>>>>>>>>>>>>>>>NeighOdd :: OMP %d , RANK(%d)=%d, MYRANK=%d (DECR %d)\n",omp_thread,omp_thread,
				this->MPI_communicator.RankNeighbour[omp_thread],
				this->MPI_communicator.getRank(),DECR);

			/* We first receive and then send */
			/*MPI_Irecv((*ElectricNodes_toRecv),size1_recv*size2_recv*SIZE_VEC,
				MPI_DOUBLE,this->MPI_communicator.RankNeighbour[omp_thread],
				omp_Neighboor,MPI_COMM_WORLD,&(*requests_MPI)[DECR+1]);*/
			MPI_Recv((*ElectricNodes_toRecv),size1_recv*size2_recv*SIZE_VEC,
				MPI_DOUBLE,this->MPI_communicator.RankNeighbour[omp_thread],
				omp_Neighboor,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			printf("MPI %d::MPI_Irecv::done (request at %d).\n",
				this->MPI_communicator.getRank(),DECR+1);

			/*MPI_Isend((*ElectricNodes_toSend),size1_send*size2_send*SIZE_VEC,MPI_DOUBLE,
				this->MPI_communicator.RankNeighbour[omp_thread],omp_thread,
				MPI_COMM_WORLD,&(*requests_MPI)[DECR+0]);*/
			MPI_Send((*ElectricNodes_toSend),size1_send*size2_send*SIZE_VEC,MPI_DOUBLE,
				this->MPI_communicator.RankNeighbour[omp_thread],omp_thread,
				MPI_COMM_WORLD);

			printf("MPI %d::MPI_Isend::done (request at %d).\n",
				this->MPI_communicator.getRank(),DECR);

	/* NEIGHBOOR LARGER THAN ME */
	}else if( this->MPI_communicator.RankNeighbour[omp_thread] > 
		this->MPI_communicator.getRank()){
			printf("Larger :: OMP %d , RANK(%d)=%d, MYRANK=%d\n",omp_thread,omp_thread,
				this->MPI_communicator.RankNeighbour[omp_thread],
				this->MPI_communicator.getRank());

			/* First receive and then send */
			/*MPI_Irecv((*ElectricNodes_toRecv),size1_recv*size2_recv*SIZE_VEC,
				MPI_DOUBLE,this->MPI_communicator.RankNeighbour[omp_thread],
				omp_Neighboor,MPI_COMM_WORLD,&(*requests_MPI)[DECR+1]);*/
			MPI_Recv((*ElectricNodes_toRecv),size1_recv*size2_recv*SIZE_VEC,
				MPI_DOUBLE,this->MPI_communicator.RankNeighbour[omp_thread],
				omp_Neighboor,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			printf("MPI %d::MPI_Irecv::done (request at %d).\n",
				this->MPI_communicator.getRank(),DECR+1);

			/*MPI_Isend((*ElectricNodes_toSend),size1_send*size2_send*SIZE_VEC,
				MPI_DOUBLE,
				this->MPI_communicator.RankNeighbour[omp_thread],omp_thread,
				MPI_COMM_WORLD,&(*requests_MPI)[DECR+0]);*/
			MPI_Send((*ElectricNodes_toSend),size1_send*size2_send*SIZE_VEC,
				MPI_DOUBLE,
				this->MPI_communicator.RankNeighbour[omp_thread],omp_thread,
				MPI_COMM_WORLD);
			printf("MPI %d::MPI_Isend::done (request at %d).\n",
				this->MPI_communicator.getRank(),DECR);

	/* NEIGHBOOR SMALLER THAN ME */
	}else if( this->MPI_communicator.RankNeighbour[omp_thread] < 
		this->MPI_communicator.getRank()){
			printf("Smaller :: OMP %d , RANK(%d)=%d, MYRANK=%d\n",omp_thread,omp_thread,
				this->MPI_communicator.RankNeighbour[omp_thread],
				this->MPI_communicator.getRank());
			/* First send and then receive */
			/*MPI_Isend((*ElectricNodes_toSend),size1_send*size2_send*SIZE_VEC,MPI_DOUBLE,
				this->MPI_communicator.RankNeighbour[omp_thread],omp_thread,
				MPI_COMM_WORLD,&(*requests_MPI)[DECR+0]);*/
			MPI_Send((*ElectricNodes_toSend),size1_send*size2_send*SIZE_VEC,MPI_DOUBLE,
				this->MPI_communicator.RankNeighbour[omp_thread],omp_thread,
				MPI_COMM_WORLD);
			printf("MPI %d::MPI_Isend::done (request at %d).\n",
				this->MPI_communicator.getRank(),DECR);

			/*MPI_Irecv((*ElectricNodes_toRecv),size1_recv*size2_recv*SIZE_VEC,
				MPI_DOUBLE,this->MPI_communicator.RankNeighbour[omp_thread],
				omp_Neighboor,MPI_COMM_WORLD,&(*requests_MPI)[DECR+1]);*/
			MPI_Recv((*ElectricNodes_toRecv),size1_recv*size2_recv*SIZE_VEC,
				MPI_DOUBLE,this->MPI_communicator.RankNeighbour[omp_thread],
				omp_Neighboor,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			printf("MPI %d::MPI_Irecv::done (request at %d).\n",
				this->MPI_communicator.getRank(),DECR+1);

	}else{
		printf("comminucateWithNeighboor not yet implemented, abort!\n");
		std::abort();
	}

	
}