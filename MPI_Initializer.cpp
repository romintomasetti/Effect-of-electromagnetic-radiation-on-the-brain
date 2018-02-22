#include "MPI_Initializer.h"
#include <iostream>

// Constructor:
MPI_Initializer::MPI_Initializer(int argc, char *argv[],int required){
	this->required = required;
	#if DEBUG > 1
	cout << "MPI_Initializer::constructor::IN" << endl;
	#endif
	
	/* MPI initialization: */
	int provided = INT_MIN;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &(provided) );
	this->provided = provided;
	
	// MPI Query: check that the required level of thread support is effectively the 
	// provided level of thread support:
	if(this->provided.get() != this->required.get()){
		cout << "MPI_Initializer::constructor::\n\tRequired level of thread support is ";
		cout << this->required.get() << " and the provided level of thread support is ";
		cout << this->provided.get() << ".\n";
		// If the required level of thread support is higher than the provided one, abort:
		if(this->required.get() > this->provided.get()){
			// The current machine cannot provided the required level of
			// thread support. Ask what to do.
			cout << "The required level of thread support is higher than the provided one.\n";
			cout << "For the communication between the MPI processes to work flowlessly";
			cout << ", MPI_THREAD_MULTIPLE (or at least MPI_THREAD_SERIALIZED).\n";
			cout << "You have thread support ";
			
			abort();
		}
	}
	
	// Get the ID of the MPI process:
	int ID_MPI_Process = INT_MIN;
	MPI_Comm_rank( MPI_COMM_WORLD, &ID_MPI_Process);
	this->ID_MPI_Process = ID_MPI_Process;
	
	// Get the number of MPI processes:
	int number_of_MPI_Processes = INT_MIN;
	MPI_Comm_size( MPI_COMM_WORLD, &number_of_MPI_Processes);
	this->number_of_MPI_Processes = number_of_MPI_Processes;
	cout << "number of MPI processes is " << this->number_of_MPI_Processes.get() << endl;
	#if DEBUG > 1
	cout << "MPI_Initializer::constructor::OUT" << endl;
	#endif
}

// Destructor:
MPI_Initializer::~MPI_Initializer(void){
	#if DEBUG > 1
	cout << "MPI_Initializer::destructor::IN" << endl;
	#endif
	MPI_Finalize();
	#if DEBUG > 1
	cout << "MPI_Initializer::destructor::OUT" << endl;
	#endif
}

// Is this MPI process the root one?
// If the process is not the root one, than send INT_MIN.
int MPI_Initializer::isRootProcess(void){
	if(this->ID_MPI_Process.get() == ROOT_PROCESSOR){
		return this->ID_MPI_Process.get();
	}else{
		return INT_MIN;
	}
}

// Get rank/ID of the MPI process:
int MPI_Initializer::getRank(void){
	return this->ID_MPI_Process.get();
}


void MPI_Initializer::MpiDivision(GridCreator &subGrid, int nbProc, int myRank){  

	double Lx = subGrid.input_parser.lengthX;
	double Ly = subGrid.input_parser.lengthY;
	double Lz = subGrid.input_parser.lengthZ;

	double deltaX = subGrid.deltaX;
	double deltaY = subGrid.deltaY;
	double deltaZ = subGrid.deltaZ;


	int N = (int) pow(nbProc, 1.0/3.0);
	std::vector<double> mpiExtremity;


	// Cubic case
	if(N*N*N == nbProc){
		double LxLocal = Lx/N;
		double LyLocal = Ly/N;
		double LzLocal = Lz/N;
		mpiExtremity.push_back((myRank%N)*LxLocal);
		mpiExtremity.push_back(((myRank%N)+1)*LxLocal);
		mpiExtremity.push_back( (((int)(myRank/N) - (int) (myRank/(N*N))*N )) *LyLocal);
		mpiExtremity.push_back((((int)(myRank/N) - ((int) (myRank/(N*N) ))*N )+1)*LyLocal);
		mpiExtremity.push_back( myRank/ (N*N) * LzLocal);
		mpiExtremity.push_back(((myRank/(N*N))+1)*LzLocal);
	}

	//Impair case
	else if(nbProc %2 != 0){
		double LxLocal = Lx/nbProc;
		double LyLocal = Ly;
		double LzLocal = Lz;

		mpiExtremity.push_back(myRank*LxLocal);
		mpiExtremity.push_back((myRank+1)*LxLocal);
		mpiExtremity.push_back( 0);
		mpiExtremity.push_back(LyLocal);
		mpiExtremity.push_back(0);
		mpiExtremity.push_back(LzLocal);
		 // Coordinates of all subdivisions in the order -> Lx, Ly, Lz
	}

	//Pair case
	else{
		double LxLocal = 2*Lx/nbProc;
		double LyLocal = Ly/2;
		double LzLocal = Lz;
		mpiExtremity.push_back(myRank%(nbProc/2)*LxLocal);
		mpiExtremity.push_back(((myRank%(nbProc/2))+1)*LxLocal);
		mpiExtremity.push_back(((int)(2*myRank/nbProc))*LyLocal);
		mpiExtremity.push_back(((int)((2*myRank/nbProc))+1)*LyLocal);
		mpiExtremity.push_back(0);
		mpiExtremity.push_back(LzLocal); // Coordinates of all subdivisions in the order -> Lx, Ly, Lz
		cout << ((int)(2*myRank/nbProc)) << endl;

	}
    
    // transformation of global lengths to indices and save them in subdomain
	subGrid.originIndices.push_back(mpiExtremity[0]/deltaX);
	subGrid.originIndices.push_back(mpiExtremity[2]/deltaY);
	subGrid.originIndices.push_back(mpiExtremity[4]/deltaZ);
	// Length of the subdomains
	subGrid.lengthX = mpiExtremity[1]-mpiExtremity[0];
	subGrid.lengthY = mpiExtremity[3]-mpiExtremity[2];
	subGrid.lengthZ = mpiExtremity[5]-mpiExtremity[4];

}

