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

template <typename T>
size_t GetFileSize(std::string filename)
{
    FILE *f = NULL;
    f = fopen(filename.c_str() , "r");
    if( f == NULL ){
        printf("File %s cannot be opened.\n",filename.c_str());
        abort();
    }
    fseek(f, 0, SEEK_END);
    size_t len = (size_t)(ftell(f)/sizeof(T));
    fclose(f);
    return len;
}

template <typename T >
std::vector<T> readBinaryFile(std::string filename)
{

	std::ifstream in(filename, std::ios::binary | std::ios::in );

    size_t file_size = GetFileSize<T>(filename);
    
    std::vector<T> ret(file_size+1);

    in.read(reinterpret_cast<char*>(&ret[0]), ret.size()*sizeof(T));

    in.close();

    ret.pop_back();
    
    return ret;
}

int main(int argc, char *argv[]){

    // Acquire data:
	//std::string filename = "../../BrainScans/DATA/subject20_crisp_v.rawb";
    std::string filename = "../../BrainScans/DATA/materials_100by119by100.bin";
	std::vector<float> data_float = readBinaryFile<float>(filename);
    
    printf(">>> Reading binary file %s... Success !\n",filename.c_str());
    
    std::vector<double> data(data_float.size(),0);
    
    for(size_t I = 0 ; I < data.size() ; I ++){
        data[I] = (double)data_float[I];
    }

    vector<double> L = {
        0.181,
        0.217,
        0.181
    };
    
    //vector<size_t> data_size    = {362   , 434    , 362  };
    vector<size_t> data_size    = {100,119,100  };
    vector<double> data_step    = {
        1/((double)(data_size[0]-1)/L[0]),
        1/((double)(data_size[1]-1)/L[1]),
        1/((double)(data_size[2]-1)/L[2])
    };

    if( data_step[0] != data_step[1] || data_step[0] != data_step[2] || data_step[1] != data_step[2]){
        printf(">>> Spatial steps are not equal.\n");
        printf(">>> deltas(%.10lf,%.10lf,%.10lf)\n",data_step[0],data_step[1],data_step[2]);
        printf(">>> I take the mean, which is %.10lf.\n",
            accumulate( data_step.begin(), data_step.end(), 0.0)/data_step.size()); 
        double mean = accumulate( data_step.begin(), data_step.end(), 0.0)/data_step.size();
        data_step[0] = mean;
        data_step[1] = mean;
        data_step[2] = mean;
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
    export_spoints_XML("materials_100_from_outBin", 0, grid, grid, Zip::ZIPPED);

    return EXIT_SUCCESS;
}
