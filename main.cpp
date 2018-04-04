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

#include <vector>

#include "Materials.h"
#include "Array_3D_Template.h"
#include "Node3DField.h"
#include "MPI_Initializer.h"
#include "SetOnceVariable_Template.h"
#include "InputParser.h"
#include "ElectromagneticSource.h"

#include "vtl.h"
#include "vtlSPoints.h"

#include "InterfaceToParaviewer.h"

#include "GridCreator_NEW.h"
#include "AlgoElectro_NEW.hpp"

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KNRM  "\x1B[0m"
#define KBLU  "\x1B[34m"
#define KYEL  "\x1B[33m"

#include "ProfilingClass.h"

#include "MemoryUsage.hpp"



using namespace std;


void check_input_file_name_given(int argc, char *argv[],map<std::string,std::string> &inputs);	

int main(int argc, char *argv[]){

	omp_set_nested(1);
	omp_set_dynamic(0);

	map<std::string,std::string> inputs;

	check_input_file_name_given(argc, argv,inputs);	

	ProfilingClass profiler;

	/* First of all, initialize MPI because if it fails, the program must immediately be stopped. */
	MPI_Initializer MPI_communicator(argc,argv,MPI_THREAD_MULTIPLE);
	#ifndef NDEBUG
		printf("\n---------\nMPI rank is %d and isRoot %d.\n--------\n",MPI_communicator.getRank(),
			MPI_communicator.isRootProcess());
	#endif

	/* Call the input file parser, input file name given as an argument: */
	#ifndef NDEBUG
		cout << "Calling input file parser...\n";
	#endif
	string filenameInput = inputs["-inputfile"];
	InputParser input_parser;
	int MPI_RANK = MPI_communicator.getRank();
	input_parser.defaultParsingFromFile(filenameInput,MPI_RANK);
	#ifndef NDEBUG
		cout << "INPUT PARSER HAS FINISHED HIS JOBS." << endl;
	#endif
	
	/* The material object stores all the material properties */
	Materials allMat;
	allMat.getPropertiesFromFile(
		input_parser.material_data_files[0],
		input_parser.material_data_files[1]
	);
	//allMat.printAllProperties();
	/*cout << "Print number of temp per mat::IN" << endl;
	allMat.printNumberOfTempLinePerMat();
	cout << "Print number of temp per mat::OUT" << endl;
	cout << "For material 0 at 25K, property 1 is ";
	cout << allMat.getProperty(25.,0,1) << endl;*/

	
	

	
	std::string profilingName;
	profilingName = input_parser.get_outputNames()["profile"];
	profilingName.append("_RANK");
	profilingName.append(std::to_string(MPI_communicator.getRank()));
	profiler.setOutputFileName(profilingName);
	

	
	
	/*printf("Initial temperature for AIR is %f.\n",
				input_parser.GetInitTemp_FromMaterialName["AIR"]);*/	

	GridCreator_NEW gridTest(input_parser,allMat,MPI_communicator,profiler);

	//cout << "Mesh init\n";
	gridTest.meshInitialization();

	//gridTest.Display_size_fields();


	InterfaceToParaviewer interfaceToWriteOutput(
			MPI_communicator,
			gridTest
		);
		
	interfaceToWriteOutput.convertAndWriteData(0,"THERMAL");
	interfaceToWriteOutput.convertAndWriteData(0,"ELECTRO");

	
	AlgoElectro_NEW algoElectro_newTst;
	algoElectro_newTst.update(gridTest,interfaceToWriteOutput);

	profiler.probeMaxRSS();
	profiler.writeToOutputFile();
	
	#ifndef NDEBUG
		cout << "Calling all the destructors.\n";
	#endif
	
	return 0;
}

/**
 * Check that a valid input file name is given as an argument:
 */
void check_input_file_name_given(int argc, char *argv[],map<std::string,std::string> &inputs){

	for(int I = 0 ; I < argc ; I ++){
		//printf("Arg[%d] = %s\n",I,argv[I]);

		/// If argv[I] is "-inputfile", use it:
		if(strcmp(argv[I],"-inputfile") == 0 && I < argc -1 ){

			inputs.insert(std::pair<std::string,std::string>("-inputfile",argv[++I]));
			
		}else if(strcmp(argv[I],"-inputfile") == 0 && I < argc ){

			fprintf(stderr,"In %s :: ERROR :: you give '-inputfile' but nothing after!\n",
					__FUNCTION__);
			fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
			#ifdef MPI_COMM_WORLD
				MPI_Abort(MPI_COMM_WORLD,-1);
			#else
				abort();
			#endif
		}

	}

}



