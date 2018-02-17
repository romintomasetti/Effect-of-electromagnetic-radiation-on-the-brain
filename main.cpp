#include <new>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <climits>
#include <omp.h>

#include "Array_3D.h"
#include "Materials.h"
#include "Array_3D_Template.h"
#include "Node.h"

#define PARALLELISM_OMP_ENABLED 1

using namespace std;

/*
 * Questions pour Boman:
 * 	1) Utiliser double *values ou vector<double>
 * 	2) Utiliser get_value ou bien passer values en public pour faire values[...] -> Rapidité ??
 */



int main(){
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
	}
}





