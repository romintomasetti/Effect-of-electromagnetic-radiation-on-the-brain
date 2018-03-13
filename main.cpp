/*! \mainpage My Personal Index Page
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *  
 * etc...
 */



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

	//////////////////
	/**
	 * Super important for speed of output writing !
	 */
	omp_set_nested(1);
	//////////////////

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
	cout << "INPUT PARSER HAS FINISHED HIS JOBS." << endl;

	
	std::string profilingName;
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
	cout << "Mesh init\n";
	gridTest.meshInitialization();

	InterfaceToParaviewer interfaceToWriteOutput(
			mesher,
			MPI_communicator,
			gridTest,true /*is_grid_creator_new*/);
	interfaceToWriteOutput.convertAndWriteData(0,"THERMAL");
	interfaceToWriteOutput.convertAndWriteData(0,"ELECTRO");

	printf("ABORTING IN MAIN LINE 157\n");
	MPI_Barrier(MPI_COMM_WORLD);

	printf("\n\nMPI %d : Ex(%zu,%zu,%zu) | Ey(%zu,%zu,%zu) | Ez(%zu,%zu,%zu)"
					"Hx(%zu,%zu,%zu) | Hy(%zu,%zu,%zu) | Hz(%zu,%zu,%zu)\n\n",
			MPI_communicator.getRank(),
			gridTest.size_Ex[0],
			gridTest.size_Ex[1],
			gridTest.size_Ex[2],
			gridTest.size_Ey[0],
			gridTest.size_Ey[1],
			gridTest.size_Ey[2],
			gridTest.size_Ez[0],
			gridTest.size_Ez[1],
			gridTest.size_Ez[2],
			gridTest.size_Hx[0],
			gridTest.size_Hx[1],
			gridTest.size_Hx[2],
			gridTest.size_Hy[0],
			gridTest.size_Hy[1],
			gridTest.size_Hy[2],
			gridTest.size_Hz[0],
			gridTest.size_Hz[1],
			gridTest.size_Hz[2]
			);

	MPI_Barrier(MPI_COMM_WORLD);
	//MPI_Abort(MPI_COMM_WORLD,-1);


	//AlgoElectro algoElectromagn; 
	
	//algoElectromagn.update(mesher,interfaceToWriteOutput);
	
	AlgoElectro_NEW algoElectro_newTst;
	algoElectro_newTst.update(gridTest,interfaceToWriteOutput);
	
	cout << "Calling all the destructors.\n";
	
	/* Phylosophy : */
	// First create the MPI communicator
	// Read the input file
	// Then create the mesh on the mpi process
	
	return 0;
}





