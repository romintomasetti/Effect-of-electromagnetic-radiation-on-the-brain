/* This class provides all necessary functions to initialize MPI */
#ifndef MPI_INITIALIZER_H
#define MPI_INITIALIZER_H

#include <iostream>
#include <mpi.h>
#include <limits.h>
#include "SetOnceVariable_Template.h"

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
		SetOnceVariable_Template<int> ID_MPI_Process = INT_MIN;
		// Number of MPI processes:
		SetOnceVariable_Template<int> number_of_MPI_Processes = INT_MIN;
	public:
		// Constructor:
		MPI_Initializer(int argc, char *argv[], int required);
		// Destructor:
		~MPI_Initializer(void);
		// Is this MPI process the root one?
		int isRootProcess(void);
		// Get rank/ID of the MPI process:
		int getRank(void);
		
};

#endif