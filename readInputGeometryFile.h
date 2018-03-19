#ifndef READINPUTGEOMETRYFILE_H
#define READINPUTGEOMETRYFILE_H

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include <cctype>

#include <sstream>

using namespace std;

/**
 * @brief This function reads the geometry file provided as an argument.
 */
unsigned int* read_input_geometry_file(std::string filename){

	ifstream geoFile(filename.c_str());

	vector<unsigned int> vector;

	if(!geoFile.is_open()){
		printf("%s :: ERROR :: Cannot open the geometry file %s. Aborting.\n",
			__FUNCTION__,filename.c_str());
		abort();
	}

	string line;
	
	while ( getline (geoFile,line) )
	{
		/// Check that the line *STARTS* with a number. If not, continue.
		if(!isdigit(line[0])){
			std::cout << line << '\n';
			continue;
		}else{
			/// Read the line and add the numbers:
			stringstream iss( line );
			unsigned int number;
			while ( iss >> number )
  				vector.push_back( number );
		}		
	}

	
	geoFile.close();	

	return &vector[0];
}

#endif
