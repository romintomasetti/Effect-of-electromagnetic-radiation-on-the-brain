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
	// This variable can be set only once:
	SetOnceVariable_Template<int> setOnce;
	//SetOnceVariable_Template<int> setOnceBis(5);
	cout << "setOnce is " << setOnce.get() << endl;
	setOnce = 1;
	setOnce = 2;
	setOnce = 3;
	cout << "setOnce is " << setOnce.get() << endl;
	setOnce.~SetOnceVariable_Template<int>();
	/* First of all, initialize MPI because if it fails, the program must immediately be stopped. */
	MPI_Initializer MPI_communicator(argc,argv,MPI_THREAD_MULTIPLE);
	printf("\n---------\nMPI rank is %d and isRoot %d.\n--------\n",MPI_communicator.getRank(),
		  MPI_communicator.isRootProcess());
	
	// The variable allMat will store the materials' properties.
	Materials allMat;
	allMat.getPropertiesFromFile("MaterialProperties.csv");
	allMat.printAllProperties();
	allMat.printNumberOfTempLinePerMat();
	cout << "For material 0 at 25K, property 1 is ";
	cout << allMat.getProperty(25.,0,1) << endl;

	/* Each point (also called node hereinafter) in the grid has:
		1) 3 coordinates                             (3 doubles)
		2) 2*3 fields components (Ex,Ey,Ez,Hx,Hy,Hz) (6 doubles)
		3) 1 temperature                             (1 double)
		4) 1 number to determine the material        (1 unsigned char provided there are not more than 255 materials)
		
		We thus have 10*4 + 1 bytes per node.
		*/
	
	
	
	
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
	
	cout << "Calling parser...\n";
	string filenameInput = "TESTS/testSourceCenteredInCube.input";
	InputParser input_parser;
	input_parser.defaultParsingFromFile(filenameInput);
	
	cout << "Calling GridCreator constructor" << endl;
	GridCreator mesher(input_parser,allMat);
	mesher.test(1,1,1,10,10,10,1);
	mesher.meshInitialization();
	
	
	
	
	
	ElectromagneticSource source;
	source.setLengths(2,2,2);
	source.setCenter (5,5,5);
	source.computeNodesInsideSource(10,10,10,1,1,1);
	if(source.isInsideSource(5,5,5)){
		cout << "Node(5,5,5) is inside the source.\n";
	}
	if(source.isInsideSource(4,4,4)){
		cout << "Node(4,4,4) is inside the source.\n";
	}
	if(source.isInsideSource(3,3,3)){
		cout << "Node(3,3,3) is inside the source.\n";
	}else{
		cout << "Node(3,3,3) is NOT inside the source.\n";
	}
	if(source.isInsideSource(6,6,6)){
		cout << "Node(6,6,6) is inside the source.\n";
	}
	
	cout << "Calling all the destructors.\n";
	
	/* Phylosophy : */
	// First create the MPI communicator
	// Read the input file
	// Then create the mesh on the mpi process
	
	return 0;
}





