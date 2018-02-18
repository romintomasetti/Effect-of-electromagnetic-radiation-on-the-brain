#include <new>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <climits>
#include <omp.h>
#include <mpi.h>

#include "Array_3D.h"
#include "Materials.h"
#include "Array_3D_Template.h"
#include "Node.h"
#include "GridCreator.h"
#include "MPI_Initializer.h"
#include "SetOnceVariable_Template.h"

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
	cout << "setOnce is " << setOnce.get() << endl;
	setOnce = 1;
	setOnce = 2;
	setOnce = 3;
	cout << "setOnce is " << setOnce.get() << endl;
	
	/* First of all, initialize MPI because if it fails, the program must immediately be stopped. */
	MPI_Initializer MPI_communicator(argc,argv,MPI_THREAD_MULTIPLE);
	printf("\n---------\nMPI rank is %d and isRoot %d.\n--------\n",MPI_communicator.getRank(),
		  MPI_communicator.isRootProcess());
	
	
	// The variable allMat will store the materials' properties.
	Materials allMat;
	allMat.getPropertiesFromFile("MaterialProperties.csv");
	allMat.printAllProperties();

	/* Each point (also called node hereinafter) in the grid has:
		1) 3 coordinates                             (3 doubles)
		2) 2*3 fields components (Ex,Ey,Ez,Hx,Hy,Hz) (6 doubles)
		3) 1 temperature                             (1 double)
		4) 1 number to determine the material        (1 unsigned char provided there are not more than 255 materials)
		
		We thus have 10*4 + 1 bytes per node.
		*/
	// nodes is a vector of the class Array_3D_Template (3D array emulated by
	// a 1D array) containing all the nodes. Each node is an intance of the 
	// class Node and has its own properties.
	Array_3D_Template<Node> nodes(5,5,5);
	// Changing the temperature of the node (5,5,5):
	nodes(5,5,5).temperature = 5;
	cout << "Node(5,5,5) : temperature : " << nodes(5,5,5).temperature << endl;
	// Changing the field Ex of the node(5,5,5):
	nodes(5,5,5).electricField[1] = 1.2;
	cout << "Node(5,5,5) : Ex          : " << nodes(5,5,5).electricField[1] << endl;
	
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
	
	
	GridCreator mesher(allMat.get_dictionnary_MaterialToID());
	mesher.nodesInitialization(nodes);
	
	return 0;
}





