#include "readInputGeometryFile.h"

CREATE_GEOMETRY_API unsigned int* read_input_geometry_file(std::string filename, size_t *size_read){

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
			while ( iss >> number ){
  				vector.push_back( number );
			}
		}		
	}

	
	geoFile.close();


	unsigned int* vector_to_return = (unsigned int*) calloc(vector.size(),sizeof(unsigned int));
	for(size_t K = 0 ; K < vector.size() ; K++)
		vector_to_return[K] = vector[K];

	*size_read = vector.size();

	return vector_to_return;
}

