#include <new>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <climits>
#include <omp.h>
#include <mpi.h>
#include <unistd.h>

#include "Materials.h"
#include "Array_3D_Template.h"
#include "Node3DField.h"
#include "GridCreator.h"
#include "MPI_Initializer.h"
#include "SetOnceVariable_Template.h"
#include "InputParser.h"
#include "ElectromagneticSource.h"

#define PARALLELISM_OMP_ENABLED 1

using namespace std;

/*
 * Questions pour Boman:
 * 	1) Utiliser double *values ou vector<double>
 * 	2) Utiliser get_value ou bien passer values en public pour faire values[...] -> Rapidité ??
 */





int main(int argc, char *argv[]){
	/* First of all, initialize MPI because if it fails, the program must immediately be stopped. */
	MPI_Initializer MPI_communicator(argc,argv,MPI_THREAD_MULTIPLE);
	printf("\n---------\nMPI rank is %d and isRoot %d.\n--------\n",MPI_communicator.getRank(),
		  MPI_communicator.isRootProcess());
	
	/* The variable allMat will store the materials' properties */
	Materials allMat;
	allMat.getPropertiesFromFile("MaterialProperties.csv");
	allMat.printAllProperties();
	allMat.printNumberOfTempLinePerMat();
	cout << "For material 0 at 25K, property 1 is ";
	cout << allMat.getProperty(25.,0,1) << endl;
	
	/* Small OPENMP example */
	if(PARALLELISM_OMP_ENABLED)
		cout << "OMP enabled.\n";
	
	#pragma omp parallel if(PARALLELISM_OMP_ENABLED)
	{
		#pragma omp for
		for(int i = 0 ; i < 10 ; i ++)
			printf("%d, ",i);
		#pragma omp barrier
	}
	cout << endl;
	
	cout << "Calling input file parser...\n";
	string filenameInput = "TESTS/testSourceCenteredInCube.input";
	InputParser input_parser;
	input_parser.defaultParsingFromFile(filenameInput);
	
	map<std::string,std::string> test222 = input_parser.get_outputNames();
	cout << test222["output"] << test222["error"] << test222["profile"] << endl;

	cout << "Calling GridCreator constructor" << endl;
	GridCreator mesher(input_parser,allMat,MPI_communicator);
	mesher.meshInitialization();
	
	
	
	
	
	
	
	cout << "Calling all the destructors.\n";
	
	/* Phylosophy : */
	// First create the MPI communicator
	// Read the input file
	// Then create the mesh on the mpi process
	
	return 0;
}





