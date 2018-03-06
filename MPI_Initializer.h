/* This class provides all necessary functions to initialize MPI */
#ifndef MPI_INITIALIZER_H
#define MPI_INITIALIZER_H

#include <iostream>
#include <mpi.h>
#include <limits.h>
#include "SetOnceVariable_Template.h"
#include "GridCreator.h"
#include "GridCreator_NEW.h"

class GridCreator_NEW;

using namespace std;

// ID of the root processor:
const int ROOT_PROCESSOR = 0;

class MPI_Initializer{
	private:
		// Level of provided thread support:
		SetOnceVariable_Template<int> provided;
		// Level of required thread support:
		SetOnceVariable_Template<int> required;
		// ID of the MPI process:
		SetOnceVariable_Template<int> ID_MPI_Process;
		// Number of MPI processes:
		SetOnceVariable_Template<int> number_of_MPI_Processes;
	public:
		const int rootProcess = ROOT_PROCESSOR;
		// Rank of the MPI neighboors:
		int RankNeighbour[6];

		// Constructor:
		MPI_Initializer(int argc, char *argv[], int required);
		// Destructor:
		~MPI_Initializer(void);
		// Is this MPI process the root one?
		int isRootProcess(void);
		// Get rank/ID of the MPI process:
		int getRank(void);

		int getNumberOfMPIProcesses(void){
			if(this->number_of_MPI_Processes.get_alreadySet() == true)
				return this->number_of_MPI_Processes.get();
			else{
				printf("MPI_Initializer::getNumberOfMPIProcesses::ERROR\n");
				printf("The number of MPI processes has not been set yet.\n");
				printf("Aborting.\n");
				std::abort();
			}
		}

		void MpiDivision(GridCreator &);

		void MPI_DIVISION(GridCreator_NEW & /*subgrid*/);

		bool SendDataToNeighboor(double *,size_t,unsigned char);

};

#endif