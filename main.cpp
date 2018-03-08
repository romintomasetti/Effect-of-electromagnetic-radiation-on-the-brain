#include <new>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstdlib>
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
#include "AlgoElectro.h"

#include "vtl.h"
#include "vtlSPoints.h"

#include "InterfaceToParaviewer.h"

#include "GridCreator_NEW.h"
#include "AlgoElectro_NEW.hpp"

#define PARALLELISM_OMP_ENABLED 1

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KNRM  "\x1B[0m"
#define KBLU  "\x1B[34m"
#define KYEL  "\x1B[33m"

#include "ProfilingClass.h"

using namespace std;

/*
 * Questions pour Boman:
 * 	1) Utiliser double *values ou vector<double>
 * 	2) Utiliser get_value ou bien passer values en public pour faire values[...] -> Rapidité ??
 */





int main(int argc, char *argv[]){

	/* SET OMP_DYNAMIC */
	if(const char *omp_dynamic_env = std::getenv("OMP_DYNAMIC")){
		// Already declared. Check it is false.
		if(std::strcmp(omp_dynamic_env,"false") == 0){
			printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
		}else{
			std::string set_env = "OMP_DYNAMIC=false";
			putenv(&set_env[0]);
			printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
		}
	}else{
		// OMP_DYNAMIC was not declared. Declare it.
		std::string set_env = "OMP_DYNAMIC=false";
		putenv(&set_env[0]);
		printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
	}

	ProfilingClass profiler;

	/* First of all, initialize MPI because if it fails, the program must immediately be stopped. */
	MPI_Initializer MPI_communicator(argc,argv,MPI_THREAD_MULTIPLE);
	printf("\n---------\nMPI rank is %d and isRoot %d.\n--------\n",MPI_communicator.getRank(),
		  MPI_communicator.isRootProcess());
	
	/* The variable allMat will store the materials' properties */
	Materials allMat;
	allMat.getPropertiesFromFile("data_air.csv");//MaterialProperties.csv
	allMat.printAllProperties();
	cout << "Print number of temp per mat::IN" << endl;
	allMat.printNumberOfTempLinePerMat();
	cout << "Print number of temp per mat::OUT" << endl;
	cout << "For material 0 at 25K, property 1 is ";
	cout << allMat.getProperty(25.,0,1) << endl;

	
	/* Small OPENMP example */
	if(PARALLELISM_OMP_ENABLED)
		cout << "OMP enabled.\n";
	
	omp_set_num_threads(6);
	
	cout << "Calling input file parser...\n";
	string filenameInput = "TESTS/testSourceCenteredInCube.input";
	InputParser input_parser;
	input_parser.defaultParsingFromFile(filenameInput);

	
	std::string profilingName = string();
	profilingName = input_parser.get_outputNames()["profile"];
	profilingName.append("_RANK");
	profilingName.append(std::to_string(MPI_communicator.getRank()));
	profiler.setOutputFileName(profilingName);
	

	
	
	printf("Initial temperature for AIR is %f.\n",
				input_parser.GetInitTemp_FromMaterialName["AIR"]);

	map<std::string,std::string> test222 = input_parser.get_outputNames();
	cout << test222["output"] << test222["error"] << test222["profile"] << endl;

	//cout << "Calling GridCreator constructor" << endl;
	GridCreator mesher(input_parser,allMat,MPI_communicator);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	

	GridCreator_NEW gridTest(input_parser,allMat,MPI_communicator,profiler);
	gridTest.meshInitialization();

	InterfaceToParaviewer interfaceToWriteOutput(mesher,MPI_communicator,gridTest);
	interfaceToWriteOutput.convertAndWriteData(0);


	//AlgoElectro algoElectromagn;
	
	//algoElectromagn.update(mesher,interfaceToWriteOutput);
	
	AlgoElectro_NEW algoElectro_newTst;
	algoElectro_newTst.update(gridTest);
	
	cout << "Calling all the destructors.\n";
	
	/* Phylosophy : */
	// First create the MPI communicator
	// Read the input file
	// Then create the mesh on the mpi process
	
	return 0;
}





