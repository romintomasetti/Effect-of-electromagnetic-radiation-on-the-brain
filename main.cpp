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

#include <boost/lexical_cast.hpp>

#include <vector>

#include "Materials.h"
#include "MPI_Initializer.h"
#include "SetOnceVariable_Template.h"
#include "InputParser.h"
#include "ElectromagneticSource.h"

#include "UTILS/directory_searching.hpp"

#include "vtl_romin.h"
#include "vtlSPoints_romin.h"

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

#include "DISC_INTEGR/discrete_integration_util.hpp"

#include "boost/format.hpp"

#include "algo_thermo.hpp"

using namespace std;

void check_input_file_name_given(int argc, char *argv[],map<std::string,std::string> &inputs);	

int main(int argc, char *argv[]){
	
	unsigned int VERBOSITY = 5;

	char cwd[1024];
	if (getcwd(cwd, sizeof(cwd)) != NULL)
		fprintf(stdout, "Current working dir: %s\n", cwd);
	else
		perror("getcwd() error");
	
	omp_set_nested(1);
	omp_set_dynamic(0);

	std::map<std::string,std::string> inputs;
	std::map<std::string,std::string>::iterator it;

	check_input_file_name_given(argc, argv,inputs);	
	
	it = inputs.find("-v");
	if(it != inputs.end()){
		/// User provided a verbosity argument.
		VERBOSITY = std::stoul(inputs["-v"]);
	}

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
	InputParser input_parser(MPI_communicator.getRank());
	int MPI_RANK = MPI_communicator.getRank();
	input_parser.defaultParsingFromFile(filenameInput,MPI_RANK);
	#ifndef NDEBUG
		cout << "INPUT PARSER HAS FINISHED HIS JOBS." << endl;
	#endif
	
	/* The material object stores all the material properties */
	Materials allMat(
		VERBOSITY,
		MPI_communicator.getRank(),
		MPI_communicator.rootProcess,
		input_parser.material_data_directory.string(),
		input_parser.GetInitTemp_FromMaterialName
	);

	//allMat.printf_unified_material_list();
	//printf("Number of materials is %zu.\n",allMat.unified_material_list.size());
	//abort();

	
	std::string profilingName;
	profilingName = input_parser.get_outputNames()["profile"];
	profilingName.append("_RANK");
	profilingName.append(std::to_string(MPI_communicator.getRank()));
	profiler.setOutputFileName(profilingName);
	


	GridCreator_NEW gridTest(
		VERBOSITY,
		input_parser,
		allMat,
		MPI_communicator,
		profiler);

	//cout << "Mesh init\n";
	gridTest.meshInitialization();

	//gridTest.Display_size_fields();


	InterfaceToParaviewer interfaceToWriteOutput(
			MPI_communicator,
			gridTest
		);
		
	interfaceToWriteOutput.convertAndWriteData(0,"THERMAL");
	interfaceToWriteOutput.convertAndWriteData(0,"ELECTRO");

	
	AlgoElectro_NEW algoElectro_newTst(VERBOSITY);

	// printf("SUCCESS Coucou1\n");
	// std::vector<double> Ex (8,1); // The vector will contain 8 times the value 1
	// std::vector<double> Ey (8,1); // The vector will contain 8 times the value 1
	// std::vector<double> Ez (8,1); // The vector will contain 8 times the value 1

	// printf("SUCCESS Coucou2\n");
	// std::vector<double> Norm = algoElectro_newTst.ComputeNormE2(Ex, Ey, Ez);
	// // std::vector<double> Norm = algoElectro_newTst.ComputeNormE(gridTest);
	// printf("SUCCESS Coucou3\n");
	// // abort();

	printf("SUCCESS\n");
	printf("You successfully reach the point before the update and the norm computations\n");
	algoElectro_newTst.update(gridTest,interfaceToWriteOutput);
	printf("SUCCESS\n");
	printf("You successfully reach the point after the update and before the norm computations\n");
	std::vector<double> Norm = algoElectro_newTst.ComputeNormEsquareBIS(gridTest);
	printf("SUCCESS\n");
	printf("You successfully reach the point after the update and the norm computations\n");
	
	// Call thermal solver:
	//algo_thermo(argc, argv, gridTest);
	

	profiler.probeMaxRSS();
	profiler.writeToOutputFile();
	
	#ifndef NDEBUG
		cout << "Calling all the destructors.\n";
	#endif
	
	return EXIT_SUCCESS;
}

/**
 * Check that a valid input file name is given as an argument:
 */
void check_input_file_name_given(int argc, char *argv[],map<std::string,std::string> &inputs){

	for(int I = 0 ; I < argc ; I ++){
		if(strcmp(argv[I],"./main") == 0 || strcmp(argv[I],"main") == 0)
			continue;
		if(strstr(argv[I],"main") != NULL ){
			continue;
		}
		//printf("Arg[%d] = %s\n",I,argv[I]);

		/// If argv[I] is "-inputfile", use it:
		if(strcmp(argv[I],"-inputfile") == 0 && I < argc -1 ){

			inputs.insert(std::pair<std::string,std::string>("-inputfile",argv[++I]));
		}else if(strcmp(argv[I],"-v") == 0 && I < argc - 1){
			/// Verify that the argument to -v is a number:
			try{
				boost::lexical_cast<double>(argv[I+1]);
			}catch(boost::bad_lexical_cast &){
				DISPLAY_ERROR_ABORT("The argument to '-v' should be an unsigned int but has %s.",
									argv[I+1]);
			}
			inputs.insert(std::pair<std::string,std::string>("-v",argv[++I]));
		}else{
			DISPLAY_ERROR_ABORT(
				"The input argument %s is not known.",
				argv[I]
			);
		}

	}

}



