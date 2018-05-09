#include "vtl.h"
#include "vtlSPoints.h"
#include "vtl_spoints.h"

#include "ThreeDim_linear_interp.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>

#include <stdio.h>

#include <fstream>
#include <iterator>
#include <algorithm>

#include <cstdint>

#define TO_SIZE_T(N) N < 0 ? 0 : N

using namespace std;  

size_t GetFileSize(std::string filename, size_t sizeInBytes)
{
    FILE *f = NULL;
    f = fopen(filename.c_str() , "r");
    if( f == NULL ){
        printf("File %s cannot be opened.\n",filename.c_str());
        abort();
    }
    fseek(f, 0, SEEK_END);
    size_t len = (size_t)ftell(f)/sizeInBytes;
    fclose(f);
    return len;
}

template <typename T >
std::vector<T> readBinaryFile(std::string filename, size_t sizeInBytes)
{

	std::ifstream in(filename, std::ios::binary | std::ios::in );

	size_t file_size = GetFileSize(filename, sizeInBytes);

	std::vector<T> ret(file_size+1);

	char current;
	size_t counter = 0;
	while( in.good() ) {
		in.read((char*)&current,1);
		ret[counter] = (T)current;
		counter++;
	}	

	ret.pop_back();

    return ret;
}

int main(int argc, char *argv[]){

    // Acquire data:
	std::string filename = "../../BrainScans/DATA/subject20_crisp_v.rawb";
	std::vector<double> data = readBinaryFile<double>(filename,sizeof(uint8_t));

    vector<double> data_step    = {0.5E-3,0.5E-3  ,0.5E-3};
    vector<size_t> data_size    = {362   , 434    , 362  };

    if( data_step[0] != data_step[1] || data_step[0] != data_step[2] || data_step[1] != data_step[2]){
        printf("Spatial steps must be equal.\n");
        abort();
    }
    if( data.size() != data_size[0] * data_size[1] * data_size[2] ){
        printf("Data.size() == %zu but I want %zu.\n",data.size(),data_size[0] * data_size[1] * data_size[2]);
        abort();
    }

    // Grid object used to save the binary file to a .vti file:
    SPoints grid;
 
    // Initialize global grid:
    grid.o   = Vec3d(0.,0.,0.);

    grid.np1 = Vec3i(0,0,0);

    grid.np2 = Vec3i(
        data_size[0],data_size[1],data_size[2]
    );

    grid.dx = Vec3d(
        data_step[0],data_step[1],data_step[2]
    );

    grid.cscalars["material"] = &data;

    printf("> Number of voxels in the binary file is %zu, number of cells in grid is %zu.\n",
        data.size(),grid.nbc());

    // Save material:
    export_spoints_XML("materials", 0, grid, grid, Zip::ZIPPED);

    return EXIT_SUCCESS;
}
