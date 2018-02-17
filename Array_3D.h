#include <new>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

class Array_3D{
	// Size along x,y,z 
	vector<unsigned long> SIZES = {0,0,0};
	// Vector containing all the values.
	double *values = NULL;
	public:
		// Initialize the sizes:
		string set_sizes(unsigned long, unsigned long, unsigned long);
		string get_sizes(vector<unsigned long> &);
		// Allocate memory and assign a default value.
		string init(double);
		// Get/set methods:
		double get_value(unsigned long, unsigned long, unsigned long);
		bool   set_value(unsigned long, unsigned long, unsigned long, double);
		// Free memory:
		string free(void);
};
