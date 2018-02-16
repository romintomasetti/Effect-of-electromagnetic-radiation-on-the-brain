#include <new>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Array_3D.h"
#include "Materials.h"

using namespace std;


/*
 * Questions pour Boman:
 * 	1) Utiliser double *values ou vector<double>
 * 	2) Utiliser get_value ou bien passer values en public pour faire values[...] -> Rapidité ??
 */



int main(){
	Materials allMat;
	allMat.getPropertiesFromFile("MaterialProperties.csv");
	allMat.printAllProperties();

	Array_3D testArray;
	
	testArray.set_sizes(2,2,2);
	testArray.init(2.0);
	cout << "Values(2,2,2) is " << testArray.get_value(1,1,1) << "\n";
	testArray.set_value(1,1,1,5.2);
	cout << "Values(2,2,2) is " << testArray.get_value(1,1,1) << "\n";
	
	vector<unsigned long> A = {0,0,0};
	
	cout << " Values(2,2,2) is " << testArray.get_value(1,1,1) << "\n";

		
	
	testArray.get_sizes(A);
	cout << "A : " << A[0] << A[1] << A[2] << "\n";
	
	cout << "Values(2,2,2) is " << testArray.get_value(1,1,1) << "\n";

	return EXIT_SUCCESS;
}





