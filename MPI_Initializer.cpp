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
			if(this->provided.get() == MPI_THREAD_SINGLE ){
				cout << "MPI_THREAD_SINGLE. Aborting because it is too low.\n";
				abort();
			}else if(this->provided.get() == MPI_THREAD_FUNNELED ){
				cout << "MPI_THREAD_FUNNELED. Aborting because it is too low.\n";
				abort();
			}else if(this->provided.get() == MPI_THREAD_SERIALIZED ){
				cout << "MPI_THREAD_SERIALIZED. The communications will be severely affected";
				cout << " but it is OK.\n";
			}else if(this->required.get() == MPI_THREAD_MULTIPLE ){
				cout << "MPI_THREAD_MULTIPLE. This is the best thread support you can have.";
				cout << " We will intensively rely on that to communicate.\n";
			}
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


void MPI_Initializer::MpiDivision(GridCreator &subGrid){  

	// Retrieve the number of MPI proceses and the ID of the current MPI process:
	int nbProc = this->number_of_MPI_Processes.get();
	int myRank = this->ID_MPI_Process.get();

	// Retrieve the length of the whole domain along each direction:
	double Lx = subGrid.input_parser.lengthX_WholeDomain_Electro;
	double Ly = subGrid.input_parser.lengthY_WholeDomain_Electro;
	double Lz = subGrid.input_parser.lengthZ_WholeDomain_Electro;

	// Retrieve the spatial step of the electromagnetic grid:
	double deltaX = subGrid.input_parser.deltaX_Electro;
	double deltaY = subGrid.input_parser.deltaY_Electro;
	double deltaZ = subGrid.input_parser.deltaZ_Electro;

	size_t nbr_nodes_X = (size_t) Lx / deltaX + 1;
	size_t nbr_nodes_Y = (size_t) Ly / deltaY + 1;
	size_t nbr_nodes_Z = (size_t) Lz / deltaZ + 1;


	int N = (int) pow(nbProc, 1.0/3.0);
	std::vector<double> mpiExtremity;

	/* 
	RankNeighbour[0] = SOUTH (along the opposite direction of the x-axis)
	RankNeighbour[1] = NORTH (along the direction of the x-axis)
	RankNeighbour[2] = WEST (along the opposite direction of the y-axis)
	RankNeighbour[3] = EAST (along the direction of the y-axis)
	RankNeighbour[4] = DOWN (along the opposite direction of the z-axis)
	RankNeighbour[5] = UP (along the direction of the z-axis)
	*/


	// Cubic case
	if(N*N*N == nbProc){
		/*
		double LxLocal = Lx/N;
		double LyLocal = Ly/N;
		double LzLocal = Lz/N;
		mpiExtremity.push_back((myRank%N)*LxLocal);
		mpiExtremity.push_back(((myRank%N)+1)*LxLocal);
		mpiExtremity.push_back( (((int)(myRank/N) - (int) (myRank/(N*N))*N )) *LyLocal);
		mpiExtremity.push_back((((int)(myRank/N) - ((int) (myRank/(N*N) ))*N )+1)*LyLocal);
		mpiExtremity.push_back( myRank/ (N*N) * LzLocal);
		mpiExtremity.push_back(((myRank/(N*N))+1)*LzLocal);
		*/

		/* myRank%N gives the current position on the x-axis */
		/* (((int)(myRank/N) - (int) (myRank/(N*N))*N )) gives the current position on the y-axis */
		/* myRank/ (N*N) gives the current position on the z-axis */

		int PositionOnX = myRank%N;
		int PositionOnY = (((int)(myRank/N) - (int) (myRank/(N*N))*N ));
		int PositionOnZ = myRank/ (N*N);

		size_t nbr_nodes_local_X = (size_t) nbr_nodes_X / N;
		if(nbr_nodes_local_X * N != nbr_nodes_X){
			if(PositionOnX == N-1){
				subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X + (-nbr_nodes_local_X * N + nbr_nodes_X);
			}else{
				subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X;
			}
		}else{
			subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X;
		}

		size_t nbr_nodes_local_Y = (size_t) nbr_nodes_Y / N;
		if(nbr_nodes_local_Y * N != nbr_nodes_Y){
			if(PositionOnY == N-1){
				subGrid.numberOfNodesInEachDir[1] = nbr_nodes_local_Y + (-nbr_nodes_local_Y * N + nbr_nodes_Y);
			}else{
				subGrid.numberOfNodesInEachDir[1] = nbr_nodes_local_Y;
			}
		}else{
			subGrid.numberOfNodesInEachDir[1] = nbr_nodes_local_Y;
		}

		size_t nbr_nodes_local_Z = (size_t) nbr_nodes_Z / N;
		if(nbr_nodes_local_Z * N != nbr_nodes_Z){
			if(PositionOnZ == N-1){
				subGrid.numberOfNodesInEachDir[2] = nbr_nodes_local_Z + (-nbr_nodes_local_Z * N + nbr_nodes_Z);
			}else{
				subGrid.numberOfNodesInEachDir[2] = nbr_nodes_local_Z;
			}
		}else{
			subGrid.numberOfNodesInEachDir[2] = nbr_nodes_local_Z;
		}
		
		subGrid.originIndices.push_back( PositionOnX * nbr_nodes_local_X );
		subGrid.originIndices.push_back( PositionOnY * nbr_nodes_local_Y );
		subGrid.originIndices.push_back( PositionOnZ * nbr_nodes_local_Z );

		if(N == 1)
		{
			this->RankNeighbour[0] = -1;
			this->RankNeighbour[1] = -1;
			this->RankNeighbour[2] = -1;
			this->RankNeighbour[3] = -1;
			this->RankNeighbour[4] = -1;
			this->RankNeighbour[5] = -1;
		}
		else
		{
		/* --------------------------------------------------------------------------------------------------- */
		/* We do the x component */
			if(PositionOnX == 0)
			{
				this->RankNeighbour[0] = -1;
				this->RankNeighbour[1] = myRank+1;
			}
			else if(PositionOnX == N-1)
			{
				this->RankNeighbour[0] = myRank-1;
				this->RankNeighbour[1] = -1;
			}
			else
			{
				this->RankNeighbour[0] = myRank-1;
				this->RankNeighbour[1] = myRank+1;
			}

			/* We do the y component */
			if(PositionOnY == 0)
			{
				this->RankNeighbour[2] = -1;
				this->RankNeighbour[3] = myRank+N; 
			}
			else if(PositionOnY == N-1)
			{
				this->RankNeighbour[2] = myRank-N;
				this->RankNeighbour[3] = -1;
			}
			else
			{
				this->RankNeighbour[2] = myRank-N;
				this->RankNeighbour[3] = myRank+N;
			}

			/* We do the z component */
			if(PositionOnZ == 0)
			{
				this->RankNeighbour[4] = -1;
				this->RankNeighbour[5] = myRank+N*N;
			}
			else if(PositionOnZ == N-1)
			{
				this->RankNeighbour[4] = myRank-N*N;
				this->RankNeighbour[5] = -1;
			}
			else
			{
				this->RankNeighbour[4] = myRank-N*N;
				this->RankNeighbour[5] = myRank+N*N;
			}
			/* --------------------------------------------------------------------------------------------------------------- */
		}
	}

	// ODD case
	else if(nbProc %2 != 0){

		// Coordinates of all subdivisions in the order -> Lx, Ly, Lz
		/*
		mpiExtremity.push_back(myRank*LxLocal);
		mpiExtremity.push_back((myRank+1)*LxLocal);
		mpiExtremity.push_back( 0);
		mpiExtremity.push_back(LyLocal);
		mpiExtremity.push_back(0);
		mpiExtremity.push_back(LzLocal);
		*/

		int PositionOnX = myRank;
		int PositionOnY = 0;
		int PositionOnZ = 0;

		size_t nbr_nodes_local_X = (size_t) nbr_nodes_X / nbProc;
		if(nbr_nodes_local_X * nbProc != nbr_nodes_X){
			if(PositionOnX == nbProc-1){
				subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X + nbr_nodes_X - nbr_nodes_local_X*nbProc;
			}else{
				subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X;
			}
		}else{
			subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X;
		}

		subGrid.numberOfNodesInEachDir[1] = nbr_nodes_Y;
		subGrid.numberOfNodesInEachDir[2] = nbr_nodes_Z;

		subGrid.originIndices.push_back( PositionOnX * nbr_nodes_local_X );
		subGrid.originIndices.push_back( PositionOnY * nbr_nodes_Y );
		subGrid.originIndices.push_back( PositionOnZ * nbr_nodes_Z );

		/* We do the x component */
		if(PositionOnX == 0)
		{
			this->RankNeighbour[0] = -1;
			this->RankNeighbour [1] = myRank+1;
		}
		else if(PositionOnX == nbProc-1)
		{
			this->RankNeighbour[0] = myRank-1;
			this->RankNeighbour[1] = -1;
		}
		else
		{
			this->RankNeighbour[0] = myRank-1;
			this->RankNeighbour[1] = myRank+1;
		}

		/* We put -1 everywhere because we have only a separation along the x-axis */
		this->RankNeighbour[2] = -1;
		this->RankNeighbour[3] = -1;
		this->RankNeighbour[4] = -1;
		this->RankNeighbour[5] = -1;

	}

	// EVEN case
	else{
		/*
		mpiExtremity.push_back(myRank%(nbProc/2)*LxLocal);
		mpiExtremity.push_back(((myRank%(nbProc/2))+1)*LxLocal);
		mpiExtremity.push_back(((int)(2*myRank/nbProc))*LyLocal);
		mpiExtremity.push_back(((int)((2*myRank/nbProc))+1)*LyLocal);
		mpiExtremity.push_back(0);
		mpiExtremity.push_back(LzLocal); // Coordinates of all subdivisions in the order -> Lx, Ly, Lz
		*/

		int PositionOnX = myRank%(nbProc/2);
		int PositionOnY = ((int)(2*myRank/nbProc));
		int PositionOnZ = 0;

		size_t nbr_nodes_local_X = (size_t) 2*nbr_nodes_X / (nbProc);
		size_t nbr_nodes_local_Y = (size_t) nbr_nodes_Y / 2;
		size_t nbr_nodes_local_Z = nbr_nodes_Z;

		if(nbr_nodes_local_X * (nbProc/2) != nbr_nodes_X){
			if(PositionOnX == nbProc/2 - 1){
				subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X + nbr_nodes_X - nbr_nodes_local_X * (nbProc/2);
			}else{
				subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X;
			}
		}else{
			subGrid.numberOfNodesInEachDir[0] = nbr_nodes_local_X;
		}

		if(nbr_nodes_local_Y * 2 != nbr_nodes_X){
			if(PositionOnY == 1){
				subGrid.numberOfNodesInEachDir[1] = nbr_nodes_local_Y + nbr_nodes_Y - nbr_nodes_local_Y * 2;
			}else{
				subGrid.numberOfNodesInEachDir[1] = nbr_nodes_local_Y;
			}
		}else{
			subGrid.numberOfNodesInEachDir[1] = nbr_nodes_local_Y;
		}

		subGrid.numberOfNodesInEachDir[2] = nbr_nodes_local_Z;

		subGrid.originIndices.push_back( PositionOnX * nbr_nodes_local_X );
		subGrid.originIndices.push_back( PositionOnY * nbr_nodes_local_Y );
		subGrid.originIndices.push_back( PositionOnZ * nbr_nodes_local_Z );

		/* We do the x component */
		if(nbProc == 2){
			this->RankNeighbour[0] = -1;
			this->RankNeighbour[1] = -1;
			cout << myRank <<"!!!!!!!!!!!!!!!!!!!!!! NTP !!! RankNeighbour x: " <<this->RankNeighbour[0] << " ; " <<this->RankNeighbour[1] << endl;
		}
		else{
			if(PositionOnX == 0)
			{
				this->RankNeighbour[0] = -1;
				this->RankNeighbour[1] = myRank + 1;
				
				cout << myRank <<"!!!!!!!!!!!!!!!!!!!!!! NTM !!! RankNeighbour x: " <<this->RankNeighbour[0] << " ; " <<this->RankNeighbour[1] << endl;

			}
			else if(PositionOnX == nbProc/2 -1)
			{
				this->RankNeighbour[0] = myRank-1;
				this->RankNeighbour[1] = -1;
				cout << myRank <<"!!!!!!!!!!!!!!!!!!!!!!RankNeighbour x: " <<this->RankNeighbour[0] << " ; "<<this->RankNeighbour[1] << endl;	
			}
			else
			{
				this->RankNeighbour[0] = myRank-1;
				this->RankNeighbour[1] = myRank+1;
				cout << myRank<< "!!!!!!!!!!!!!!!!!!!!!!RankNeighbour x: " <<this->RankNeighbour[0] << " ; "<<this->RankNeighbour[1] << endl;
			}
		}
		/* We do the y component */
		if(PositionOnY == 0)
		{
			this->RankNeighbour[2] = -1;
			this->RankNeighbour[3] = myRank+(nbProc/2);
			cout << myRank <<"!!!!!!!!!!!!!!!!!!!!!!RankNeighbour y: " <<this->RankNeighbour[2] << " ; "<< this->RankNeighbour[3] << endl;
		}
		else if(PositionOnY == 1)
		{
			this->RankNeighbour[2] = myRank-(nbProc/2);
			this->RankNeighbour[3] = -1;
			cout << myRank <<"!!!!!!!!!!!!!!!!!!!!!!RankNeighbour y: " <<this->RankNeighbour[2] << " ; "<<this->RankNeighbour[3] << endl;
		}
		else
		{
			printf("MPI_Initializer::MpiDivision\n");
			printf("This case is not envisageable, at line %d in file %s.\n",__LINE__,__FILE__);
			abort();

		}

		/* We do the z component */
		this->RankNeighbour[4] = -1;
		this->RankNeighbour[5] = -1;
		cout << myRank<< "!!!!!!!!!!!!!!!!!!!!!!RankNeighbour z: " <<this->RankNeighbour[4] << " ; "<<this->RankNeighbour[5] << endl;

	}

    // transformation of global lengths to indices and save them in subdomain
	/*
	subGrid.originIndices.push_back(round(mpiExtremity[0]/deltaX));
	subGrid.originIndices.push_back(round(mpiExtremity[2]/deltaY));
	subGrid.originIndices.push_back(round(mpiExtremity[4]/deltaZ));
	*/
	cout << myRank<< "---------------------->  x: " <<subGrid.originIndices[0] << endl;
	cout << myRank<< "---------------------->  y: " <<subGrid.originIndices[1] << endl;
	cout << myRank<< "---------------------->  z: " <<subGrid.originIndices[2] << endl;

	// Length of the subdomains
	/*subGrid.lengthX = mpiExtremity[1]-mpiExtremity[0];
	subGrid.lengthY = mpiExtremity[3]-mpiExtremity[2];
	subGrid.lengthZ = mpiExtremity[5]-mpiExtremity[4];*/
	

}

bool MPI_Initializer::SendDataToNeighboor(double *vectorToSend,
										  size_t lengthToSend,
										  unsigned char direction){
	// We can communicate in 6 directions
	// By default, it is assumed that everything worked fine.
	// Thus, return true.
	return true;
}

void MPI_Initializer::MPI_DIVISION(GridCreator_NEW & subGrid){
	
	// Retrieve the number of MPI proceses and the ID of the current MPI process:
	int nbProc = subGrid.MPI_communicator.number_of_MPI_Processes.get();
	int myRank = subGrid.MPI_communicator.ID_MPI_Process.get();

	// Retrieve the length of the whole domain along each direction, EM grid:
	double Lx = subGrid.input_parser.lengthX_WholeDomain_Electro;
	double Ly = subGrid.input_parser.lengthY_WholeDomain_Electro;
	double Lz = subGrid.input_parser.lengthZ_WholeDomain_Electro;

	// Retrieve the length of the whole domain along each direction, THERMAL grid:
	double Lx_thermal = subGrid.input_parser.lengthX_WholeDomain_Thermal;
	double Ly_thermal = subGrid.input_parser.lengthY_WholeDomain_Thermal;
	double Lz_thermal = subGrid.input_parser.lengthZ_WholeDomain_Thermal;

	// Retrieve the spatial step of the electromagnetic grid:
	double deltaX = subGrid.input_parser.deltaX_Electro;
	double deltaY = subGrid.input_parser.deltaY_Electro;
	double deltaZ = subGrid.input_parser.deltaZ_Electro;

	// Retrieve the spatial step of the thermal grid:
	double delta_thermal = subGrid.input_parser.delta_Thermal;

	size_t nbr_nodes_X =  Lx / deltaX + 1;
	size_t nbr_nodes_Y =  Ly / deltaY + 1;
	size_t nbr_nodes_Z =  Lz / deltaZ + 1;

	size_t nbr_nodes_thermal_X = Lx_thermal / delta_thermal + 1;
	size_t nbr_nodes_thermal_Y = Ly_thermal / delta_thermal + 1;
	size_t nbr_nodes_thermal_Z = Lz_thermal / delta_thermal + 1;


	int N = (int) pow(nbProc, 1.0/3.0);

	/* 
	 * RankNeighbour[0] = SOUTH (along the opposite direction of the x-axis)
	 * RankNeighbour[1] = NORTH (along the direction of the x-axis)
	 * RankNeighbour[2] = WEST (along the opposite direction of the y-axis)
	 * RankNeighbour[3] = EAST (along the direction of the y-axis)
	 * RankNeighbour[4] = DOWN (along the opposite direction of the z-axis)
	 * RankNeighbour[5] = UP (along the direction of the z-axis) 
	 */


	// Cubic case
	if(N*N*N == nbProc){

		/**
		 * 1) For the X component of the electric field, if I am not the last MPI process along the 
		 * 		X direction, add one node to E_x in the X direction.
		 * 2) For the Y component of the electric field, add one node in the Y direction if I am
		 * 		*NOT* the last MPI process in the Y direction.
		 * 3) For the Z component of the electric field, if I am not the last MPI process along the
		 * 		Z direction, add one node to E_z along Z.
		 * 4) For the X component of the magnetic field, add one node along Y and Z if the MPI process
		 * 		is not the last one along Y or Z, respectively.
		 * 5) For the Y component of the magnetic field, if I am not the last MPI process in the 
		 * 		X or Z direction, add one node in the X or Z direction, respectively.
		 * 6) For the Z component of the magnetic field, if I am not the last MPI process in the
		 * 		X or Y direction, add one node in the X or Y direction, respectively.
		 */

		/* myRank%N gives the current position on the x-axis */
		/* (((int)(myRank/N) - (int) (myRank/(N*N))*N )) gives the current position on the y-axis */
		/* myRank/ (N*N) gives the current position on the z-axis */

		int PositionOnX = myRank%N;
		int PositionOnY = (((int)(myRank/N) - (int) (myRank/(N*N))*N ));
		int PositionOnZ = myRank/ (N*N);

		subGrid.MPI_communicator.MPI_POSITION[0] = (char) PositionOnX;
		subGrid.MPI_communicator.MPI_POSITION[1] = (char) PositionOnY;
		subGrid.MPI_communicator.MPI_POSITION[2] = (char) PositionOnZ;

		/// In the following, you must do (nbProc-1) because MPI rank starts at zero !
		subGrid.MPI_communicator.MPI_MAX_POSI[0] = (char) ( (nbProc-1) % N );
		subGrid.MPI_communicator.MPI_MAX_POSI[1] = (char) ( (int) ( (nbProc-1) / N) 
																- (int) ( (nbProc-1) / (N*N) ) * N );
		subGrid.MPI_communicator.MPI_MAX_POSI[2] = (char) ( (nbProc-1) / (N*N));

		size_t nbr_nodes_local_X         = nbr_nodes_X / N;
		size_t nbr_nodes_local_thermal_X = nbr_nodes_thermal_X / N;

		if(nbr_nodes_local_X * N != nbr_nodes_X){
			/// If I am the last MPI, I take the remainder of the 'not perfect' division.
			if(PositionOnX == N-1){
				subGrid.sizes_EH[0] = nbr_nodes_local_X + 
					(-nbr_nodes_local_X * N + nbr_nodes_X);
			}else{
				subGrid.sizes_EH[0] = nbr_nodes_local_X;
			}
		}else{
			subGrid.sizes_EH[0] = nbr_nodes_local_X;
		}

		if(nbr_nodes_local_thermal_X * N != nbr_nodes_thermal_X){
			if(PositionOnX == N-1){
				subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X +
						(-nbr_nodes_local_thermal_X * N + nbr_nodes_thermal_X);
			}else{
				subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X;
			}
		}else{
			subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X;
		}

		size_t nbr_nodes_local_Y         = nbr_nodes_Y / N;
		size_t nbr_nodes_local_thermal_Y = nbr_nodes_thermal_Y / N;

		if(nbr_nodes_local_Y * N != nbr_nodes_Y){
			if(PositionOnY == N-1){
				subGrid.sizes_EH[1] = nbr_nodes_local_Y + 
					(-nbr_nodes_local_Y * N + nbr_nodes_Y);
			}else{
				subGrid.sizes_EH[1] = nbr_nodes_local_Y;
			}
		}else{
			subGrid.sizes_EH[1] = nbr_nodes_local_Y;
		}

		if(nbr_nodes_local_thermal_Y * N != nbr_nodes_thermal_Y){
			if(PositionOnY == N-1){
				subGrid.size_Thermal[1] = nbr_nodes_local_thermal_Y +
						(-nbr_nodes_local_thermal_Y * N + nbr_nodes_thermal_Y);
			}else{
				subGrid.size_Thermal[1] = nbr_nodes_local_thermal_Y;
			}
		}else{
			subGrid.size_Thermal[1] = nbr_nodes_local_thermal_Y;
		}

		size_t nbr_nodes_local_Z         = nbr_nodes_Z / N;
		size_t nbr_nodes_local_thermal_Z = nbr_nodes_thermal_Z / N;

		if(nbr_nodes_local_Z * N != nbr_nodes_Z){
			if(PositionOnZ == N-1){
				subGrid.sizes_EH[2] = nbr_nodes_local_Z + 
					(-nbr_nodes_local_Z * N + nbr_nodes_Z);
			}else{
				subGrid.sizes_EH[2] = nbr_nodes_local_Z;
			}
		}else{
			subGrid.sizes_EH[2] = nbr_nodes_local_Z;
		}

		if(nbr_nodes_local_thermal_Z * N != nbr_nodes_thermal_Z){
			if(PositionOnZ == N-1){
				subGrid.size_Thermal[2] = nbr_nodes_local_thermal_Z +
						(-nbr_nodes_local_thermal_Z * N + nbr_nodes_thermal_Z);
			}else{
				subGrid.size_Thermal[2] = nbr_nodes_local_thermal_Z;
			}
		}else{
			subGrid.size_Thermal[2] = nbr_nodes_local_thermal_Z;
		}
		
		// Origin indices for the electromagnetic grid:
		subGrid.originIndices_Electro[0] = PositionOnX * nbr_nodes_local_X;
		subGrid.originIndices_Electro[1] = PositionOnY * nbr_nodes_local_Y;
		subGrid.originIndices_Electro[2] = PositionOnZ * nbr_nodes_local_Z;
		// Origin indices for the thermal grid:
		subGrid.originIndices_Thermal[0] = PositionOnX * nbr_nodes_local_thermal_X;
		subGrid.originIndices_Thermal[1] = PositionOnY * nbr_nodes_local_thermal_Y;
		subGrid.originIndices_Thermal[2] = PositionOnZ * nbr_nodes_local_thermal_Z;

		if(N == 1)
		{
			this->RankNeighbour[0] = -1;
			this->RankNeighbour[1] = -1;
			this->RankNeighbour[2] = -1;
			this->RankNeighbour[3] = -1;
			this->RankNeighbour[4] = -1;
			this->RankNeighbour[5] = -1;
		}
		else
		{
		/* --------------------------------------------------------------------------------------------------- */
		/* We do the x component */
			if(PositionOnX == 0)
			{
				this->RankNeighbour[0] = -1;
				this->RankNeighbour[1] = myRank+1;
			}
			else if(PositionOnX == N-1)
			{
				this->RankNeighbour[0] = myRank-1;
				this->RankNeighbour[1] = -1;
			}
			else
			{
				this->RankNeighbour[0] = myRank-1;
				this->RankNeighbour[1] = myRank+1;
			}

			/* We do the y component */
			if(PositionOnY == 0)
			{
				this->RankNeighbour[2] = -1;
				this->RankNeighbour[3] = myRank+N; 
			}
			else if(PositionOnY == N-1)
			{
				this->RankNeighbour[2] = myRank-N;
				this->RankNeighbour[3] = -1;
			}
			else
			{
				this->RankNeighbour[2] = myRank-N;
				this->RankNeighbour[3] = myRank+N;
			}

			/* We do the z component */
			if(PositionOnZ == 0)
			{
				this->RankNeighbour[4] = -1;
				this->RankNeighbour[5] = myRank+N*N;
			}
			else if(PositionOnZ == N-1)
			{
				this->RankNeighbour[4] = myRank-N*N;
				this->RankNeighbour[5] = -1;
			}
			else
			{
				this->RankNeighbour[4] = myRank-N*N;
				this->RankNeighbour[5] = myRank+N*N;
			}
			/* --------------------------------------------------------------------------------------------------------------- */
		}
	}

	// ODD case
	else if(nbProc %2 != 0){

		/**
		 *  If the number of MPI process is odd, 'cut' along the X direction.
		 *  All MPI processes *MUST* add one node to the electric field along X.
		 *  *ONLY* the last MPI process writes the line of zeros.
		 *  For all MPI *BUT* the last one: 
		 * 		subGrid.MPI_communicator.must_add_one_E_X_along_XYZ[0] = true;
		 *  For the last MPI process:
		 * 		subGrid.MPI_communicator.must_add_one_E_X_along_XYZ[0] = false; (default value)
		 *  Regarding the magnetic field:
		 * 			1) For the X magnetic field, nothing to do.
		 * 			2) For the Y magnetic field, must add one in the X direction 
		 * 				for all *BUT* the last MPI process.
		 * 			3) For the Z magnetic field, must add one in the X direction
		 * 				for all *BUT* the last MPI process.
		 */

		int PositionOnX = myRank;
		int PositionOnY = 0;
		int PositionOnZ = 0;

		subGrid.MPI_communicator.MPI_POSITION[0] = (char) PositionOnX;
		subGrid.MPI_communicator.MPI_POSITION[1] = (char) PositionOnY;
		subGrid.MPI_communicator.MPI_POSITION[2] = (char) PositionOnZ;

		subGrid.MPI_communicator.MPI_MAX_POSI[0] = nbProc-1;
		subGrid.MPI_communicator.MPI_MAX_POSI[1] = 0;
		subGrid.MPI_communicator.MPI_MAX_POSI[2] = 0;


		size_t nbr_nodes_local_X         = nbr_nodes_X / nbProc;
		size_t nbr_nodes_local_thermal_X = nbr_nodes_thermal_X / nbProc;

		if(nbr_nodes_local_X * nbProc != nbr_nodes_X){
			if(PositionOnX == nbProc-1){
				subGrid.sizes_EH[0] = nbr_nodes_local_X + 
					nbr_nodes_X - nbr_nodes_local_X*nbProc;
			}else{
				subGrid.sizes_EH[0] = nbr_nodes_local_X;
			}
		}else{
			subGrid.sizes_EH[0] = nbr_nodes_local_X;
		}

		if(nbr_nodes_local_thermal_X * nbProc != nbr_nodes_thermal_X){
			if(PositionOnX == nbProc-1){
				subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X
						+ nbr_nodes_thermal_X - nbr_nodes_local_thermal_X * nbProc;
			}else{
				subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X;
			}
		}else{
			subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X;
		}

		subGrid.sizes_EH[1] = nbr_nodes_Y;
		subGrid.sizes_EH[2] = nbr_nodes_Z;

		subGrid.size_Thermal[1] = nbr_nodes_thermal_Y;
		subGrid.size_Thermal[2] = nbr_nodes_thermal_Z;

		// Origin indices for the electromagnetc grid:
		subGrid.originIndices_Electro[0] = PositionOnX * nbr_nodes_local_X;
		subGrid.originIndices_Electro[1] = PositionOnY * nbr_nodes_Y;
		subGrid.originIndices_Electro[2] = PositionOnZ * nbr_nodes_Z;
		// Origin indices for the thermal grid:
		subGrid.originIndices_Thermal[0] = PositionOnX * nbr_nodes_local_thermal_X;
		subGrid.originIndices_Thermal[1] = PositionOnY * nbr_nodes_thermal_Y;
		subGrid.originIndices_Thermal[2] = PositionOnZ * nbr_nodes_thermal_Z;


		/* We do the x component */
		if(PositionOnX == 0)
		{
			/// The first is for South neighboor, the second is for North neighboor.
			this->RankNeighbour[0] = myRank+1;
			this->RankNeighbour[1] = -1;
		}
		else if(PositionOnX == nbProc-1)
		{
			/// The first is for South neighboor, the second is for North neighboor.
			this->RankNeighbour[0] = -1;
			this->RankNeighbour[1] = myRank-1;
		}
		else
		{
			/// The first is for South neighboor, the second is for North neighboor.
			this->RankNeighbour[0] = myRank+1;
			this->RankNeighbour[1] = myRank-1;
		}

		/* We put -1 everywhere because we have only a separation along the x-axis */
		this->RankNeighbour[2] = -1;
		this->RankNeighbour[3] = -1;
		this->RankNeighbour[4] = -1;
		this->RankNeighbour[5] = -1;

	}

	// EVEN case
	else{

		int PositionOnX = myRank%(nbProc/2);
		int PositionOnY = ((int)(2*myRank/nbProc));
		int PositionOnZ = 0;

		subGrid.MPI_communicator.MPI_POSITION[0] = (char) PositionOnX;
		subGrid.MPI_communicator.MPI_POSITION[1] = (char) PositionOnY;
		subGrid.MPI_communicator.MPI_POSITION[2] = (char) PositionOnZ;

		subGrid.MPI_communicator.MPI_MAX_POSI[0] = (char) (nbProc-1)%(nbProc/2);
		subGrid.MPI_communicator.MPI_MAX_POSI[1] = (char) ((int)(2*(nbProc-1)/nbProc));
		subGrid.MPI_communicator.MPI_MAX_POSI[2] = 0;


		size_t nbr_nodes_local_X = 2*nbr_nodes_X / (nbProc);
		size_t nbr_nodes_local_Y = nbr_nodes_Y / 2;
		size_t nbr_nodes_local_Z = nbr_nodes_Z;

		size_t nbr_nodes_local_thermal_X = 2*nbr_nodes_thermal_X / nbProc;
		size_t nbr_nodes_local_thermal_Y = nbr_nodes_thermal_Y / 2;
		size_t nbr_nodes_local_thermal_Z = nbr_nodes_thermal_Z;

		if(nbr_nodes_local_X * (nbProc/2) != nbr_nodes_X){
			if(PositionOnX == nbProc/2 - 1){
				subGrid.sizes_EH[0] = nbr_nodes_local_X + 
					nbr_nodes_X - nbr_nodes_local_X * (nbProc/2);
			}else{
				subGrid.sizes_EH[0] = nbr_nodes_local_X;
			}
		}else{
			subGrid.sizes_EH[0] = nbr_nodes_local_X;
		}

		if(nbr_nodes_local_thermal_X * (nbProc/2) != nbr_nodes_thermal_X){
			if(PositionOnX == nbProc/2 - 1){
				subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X +
						nbr_nodes_thermal_X - nbr_nodes_local_thermal_X * (nbProc/2);
			}else{
				subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X;
			}
		}else{
			subGrid.size_Thermal[0] = nbr_nodes_local_thermal_X;
		}

		if(nbr_nodes_local_Y * 2 != nbr_nodes_Y){
			if(PositionOnY == 1){
				subGrid.sizes_EH[1] = nbr_nodes_local_Y + nbr_nodes_Y 
					- nbr_nodes_local_Y * 2;
			}else{
				subGrid.sizes_EH[1] = nbr_nodes_local_Y;
			}
		}else{
			subGrid.sizes_EH[1] = nbr_nodes_local_Y;
		}


		if(nbr_nodes_local_thermal_Y * 2 != nbr_nodes_thermal_Y){
			if(PositionOnY == 1){
				subGrid.size_Thermal[1] = nbr_nodes_local_thermal_Y
						+ nbr_nodes_thermal_Y - nbr_nodes_local_thermal_Y * 2;
			}else{
				subGrid.size_Thermal[1] = nbr_nodes_local_thermal_Y;
			}
		}else{
			subGrid.size_Thermal[1] = nbr_nodes_local_thermal_Y;
		}

		subGrid.sizes_EH[2]     = nbr_nodes_local_Z;
		subGrid.size_Thermal[2] = nbr_nodes_local_thermal_Z; 

		// Origin indices for the electromagnetic grid:
		subGrid.originIndices_Electro[0] = PositionOnX * nbr_nodes_local_X;
		subGrid.originIndices_Electro[1] = PositionOnY * nbr_nodes_local_Y;
		subGrid.originIndices_Electro[2] = PositionOnZ * nbr_nodes_local_Z;
		// Origin indices for the thermal grid:
		subGrid.originIndices_Thermal[0] = PositionOnX * nbr_nodes_local_thermal_X;
		subGrid.originIndices_Thermal[1] = PositionOnY * nbr_nodes_local_thermal_Y;
		subGrid.originIndices_Thermal[2] = PositionOnZ * nbr_nodes_local_thermal_Z;

		/* We do the x component */
		if(nbProc == 2){
			this->RankNeighbour[0] = -1;
			this->RankNeighbour[1] = -1;
		}
		else{
			/// First MPI process along X direction:
			if(PositionOnX == 0)
			{
				/// South neighboor:
				this->RankNeighbour[0] = myRank + 1;
				this->RankNeighbour[1] = -1;
				

			}
			/// Last MPI process along X direction:
			else if(PositionOnX == nbProc/2 -1)
			{
				this->RankNeighbour[0] = -1;
				this->RankNeighbour[1] =  myRank-1;
				cout << myRank <<"!!!!!!!!!!!!!!!!!!!!!!RankNeighbour x: " <<this->RankNeighbour[0] << " ; "<<this->RankNeighbour[1] << endl;	
			}
			else
			{
				this->RankNeighbour[0] = myRank+1;
				this->RankNeighbour[1] = myRank-1;
				cout << myRank<< "!!!!!!!!!!!!!!!!!!!!!!RankNeighbour x: " <<this->RankNeighbour[0] << " ; "<<this->RankNeighbour[1] << endl;
			}
		}
		/* We do the y component */
		if(PositionOnY == 0)
		{
			this->RankNeighbour[2] = -1;
			this->RankNeighbour[3] = myRank+(nbProc/2);
			cout << myRank <<"!!!!!!!!!!!!!!!!!!!!!!RankNeighbour y: " <<this->RankNeighbour[2] << " ; "<< this->RankNeighbour[3] << endl;
		}
		else if(PositionOnY == 1)
		{
			this->RankNeighbour[2] = myRank-(nbProc/2);
			this->RankNeighbour[3] = -1;
			cout << myRank <<"!!!!!!!!!!!!!!!!!!!!!!RankNeighbour y: " <<this->RankNeighbour[2] << " ; "<<this->RankNeighbour[3] << endl;
		}
		else
		{
			printf("MPI_Initializer::MpiDivision\n");
			printf("This case is not envisageable, at line %d in file %s.\n",__LINE__,__FILE__);
			abort();

		}

		/* We do the z component */
		this->RankNeighbour[4] = -1;
		this->RankNeighbour[5] = -1;
		cout << myRank<< "!!!!!!!!!!!!!!!!!!!!!!RankNeighbour z: " <<this->RankNeighbour[4] << " ; "<<this->RankNeighbour[5] << endl;

	}


	/**
	 * Check if last along X.
	 */
	if(subGrid.MPI_communicator.MPI_POSITION[0] == subGrid.MPI_communicator.MPI_MAX_POSI[0]){
		/// I am the last MPI along X.

		/// For E_x:
		subGrid.MPI_communicator.must_add_one_to_E_X_along_XYZ[0] = true;

		/// For E_y:
		subGrid.MPI_communicator.must_add_one_to_E_Y_along_XYZ[0] = false;

		/// For E_z:
		subGrid.MPI_communicator.must_add_one_to_E_Z_along_XYZ[0] = false;

		/// For H_x:
		subGrid.MPI_communicator.must_add_one_to_H_X_along_XYZ[0] = false;

		/// For H_y:
		subGrid.MPI_communicator.must_add_one_to_H_Y_along_XYZ[0] = true;

		/// For H_z:
		subGrid.MPI_communicator.must_add_one_to_H_Z_along_XYZ[0] = true;

	}

	if(subGrid.MPI_communicator.MPI_POSITION[1] == subGrid.MPI_communicator.MPI_MAX_POSI[1]){
		/// I am the last along Y.

		/// For E_x:
		subGrid.MPI_communicator.must_add_one_to_E_X_along_XYZ[1] = false;
		
		/// For E_y:
		subGrid.MPI_communicator.must_add_one_to_E_Y_along_XYZ[1] = true;

		/// For E_z:
		subGrid.MPI_communicator.must_add_one_to_E_Z_along_XYZ[1] = false;

		/// For H_x:
		subGrid.MPI_communicator.must_add_one_to_H_X_along_XYZ[1] = true;

		/// For H_y:
		subGrid.MPI_communicator.must_add_one_to_H_Y_along_XYZ[1] = false;

		/// For H_z:
		subGrid.MPI_communicator.must_add_one_to_H_Z_along_XYZ[1] = true;
	}

	if(subGrid.MPI_communicator.MPI_POSITION[2] == subGrid.MPI_communicator.MPI_MAX_POSI[2]){
		/// I am the last along Y.
		
		/// For E_x:
		subGrid.MPI_communicator.must_add_one_to_E_X_along_XYZ[2] = false;
		
		/// For E_y:
		subGrid.MPI_communicator.must_add_one_to_E_Y_along_XYZ[2] = false;

		/// For E_z:
		subGrid.MPI_communicator.must_add_one_to_E_Z_along_XYZ[2] = true;

		/// For H_x:
		subGrid.MPI_communicator.must_add_one_to_H_X_along_XYZ[2] = true;

		/// For H_y:
		subGrid.MPI_communicator.must_add_one_to_H_Y_along_XYZ[2] = true;

		/// For H_z:
		subGrid.MPI_communicator.must_add_one_to_H_Z_along_XYZ[2] = false;
	}

	printf("MPI %d ->  originIndices_Electro (%zu,%zu,%zu)\n",myRank,
		subGrid.originIndices_Electro[0],
		subGrid.originIndices_Electro[1],
		subGrid.originIndices_Electro[2]);
	printf("MPI %d ->  originInices_thermal  (%zu,%zu,%zu)\n",myRank,
		subGrid.originIndices_Thermal[0],
		subGrid.originIndices_Thermal[1],
		subGrid.originIndices_Thermal[2]);

}
