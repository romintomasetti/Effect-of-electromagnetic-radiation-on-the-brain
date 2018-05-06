#include "vtl.h"
#include "vtlSPoints.h"
#include "vtl_spoints.h"

#include "ThreeDim_linear_interp.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>

#include "sys/types.h"
#include "sys/sysinfo.h"

#include <stdio.h>

#include <fstream>
#include <iterator>
#include <algorithm>

#include <cstdint>

#define TO_SIZE_T(N) N < 0 ? 0 : N

using namespace std;

bool double_is_int_and_not_zero(double val) {
   double absolute = abs( val );
   if( floor( absolute) == 0. ){
       printf("floor(%.9g) is zero.\n",floor( absolute));
       abort();
   }
   return absolute == floor(absolute);
}

size_t GetFileSize(std::string filename)
{
    FILE *f;
    f = fopen(filename.c_str() , "r");
    fseek(f, 0, SEEK_END);
    size_t len = (size_t)ftell(f);
    fclose(f);
    return len;
}

template <typename T >
std::vector<T> readBinaryFile(string filename)
{

	ifstream in(filename, std::ios::binary | std::ios::in );

	size_t file_size = GetFileSize(filename);

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

    // Get physical memory on computer:
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    size_t totalPhysMem = memInfo.totalram;
    //Multiply in next statement to avoid int overflow on right hand side...
    totalPhysMem *= memInfo.mem_unit;

    printf("Total RAM on system is %lf MBytes.\n",((double)totalPhysMem)/1E6);

    // Global grid parameters:
    vector<double> L  = {0.2  ,0.3   ,0.2   };
    vector<double> dx = {0.5E-3,0.5E-3,0.5E-3};

    // Grid object:
    SPoints grid;

    // Initialize global grid:
    grid.o   = Vec3d(0.,0.,0.);
    grid.np1 = Vec3i(0,0,0);
    // Plus one because we will use cells in Paraview:
    grid.np2 = Vec3i(
        (size_t) (L[0]/dx[0])+1,
        (size_t) (L[1]/dx[1])+1,
        (size_t) (L[2]/dx[2])+1
    );
    printf("Number of nodes in grid is %zu.\n",grid.nbc());
    if( (double)grid.nbc() * 8 > 0.9 * (double)totalPhysMem ){
        printf("Would require too much RAM. Abort.\n");
        abort();
    }
    grid.dx = Vec3d(
        dx[0],
        dx[1],
        dx[2]
    );

    // Acquire data:
	std::string filename = "../../BrainScans/subject20_crisp_v.rawb";
	std::vector<double> data = readBinaryFile<double>(filename);

    vector<double> data_step    = {0.5E-3,0.5E-3  ,0.5E-3};
    vector<size_t> data_size    = {362   , 434    , 362  };
    vector<double> data_center  = {L[0]/2,L[1]/2  ,L[2]/2};
    vector<size_t> data_centerN = {
        (size_t) (data_center[0]/grid.dx[0]),
        (size_t) (data_center[1]/grid.dx[1]),
        (size_t) (data_center[2]/grid.dx[2])
    };

    if( data_step[0] != data_step[1] || data_step[0] != data_step[2] || data_step[1] != data_step[2]){
        printf("Spatial steps must be equal.\n");
        abort();
    }
    if( grid.dx[0] != grid.dx[1] || grid.dx[0] != grid.dx[2] || grid.dx[1] != grid.dx[2]){
        printf("Spatial steps must be equal.\n");
        abort();
    }

    double ratio = data_step[0] / grid.dx[0];

    if( ! double_is_int_and_not_zero(ratio) ){
        printf("Ratio is %.25g, which is not integer.\n",ratio);
        abort();
    }else{
        printf("Ratio between grids is %zu.\n",(size_t)ratio);
    }

    vector<size_t> data_new_size = {
        data_size[0] * (size_t)ratio - 1,
        data_size[1] * (size_t)ratio - 1,
        data_size[2] * (size_t)ratio - 1
    };

    ThreeDim_linear_interp trilin(
        data,
        data_size,
        data_step,
        data_new_size,
        dx,
        (size_t)ratio
    );

    std::vector<double> temp;
    std::vector<double> &new_data = temp;
    
    if(ratio > 1){
        
        temp = trilin.trilinearInterp();
    }else{
        new_data = data;
        data_new_size = {
            data_size[0],
            data_size[1],
            data_size[2]
        };
    }

    

    

    // Creation of the material grid on each MPI process:
    size_t mynbp = grid.nbp();
    size_t mynbc = grid.nbc();
    Vec3i np     = grid.np();
    Vec3i nc     = grid.nc();

    printf("My number of nodes is (%zu,%zu,%zu).\n",
        np[0],np[1],np[2]);
    printf("My number of cells is (%zu,%zu,%zu).\n",
        nc[0],nc[1],nc[2]); 


    std::vector<double> material(mynbc);

    grid.cscalars["material"] = &material;

    size_t Nx_inf = (size_t) TO_SIZE_T( (double)data_centerN[0]-(double)(data_new_size[0])/2. );
    size_t Nx_sup = (size_t) TO_SIZE_T( (double)data_centerN[0]+(double)(data_new_size[0])/2. );
    size_t Ny_inf = (size_t) TO_SIZE_T( (double)data_centerN[1]-(double)(data_new_size[1])/2. );
    size_t Ny_sup = (size_t) TO_SIZE_T( (double)data_centerN[1]+(double)(data_new_size[1])/2. );
    size_t Nz_inf = (size_t) TO_SIZE_T( (double)data_centerN[2]-(double)(data_new_size[2])/2. );
    size_t Nz_sup = (size_t) TO_SIZE_T( (double)data_centerN[2]+(double)(data_new_size[2])/2. );

    printf("Limits of the form: (%zu,%zu,%zu,%zu,%zu,%zu).\n",
        Nx_inf,
        Nx_sup,
        Ny_inf,
        Ny_sup,
        Nz_inf,
        Nz_sup
    );

    // Fill in local material vector:
    for(size_t K = 0 ; K < nc[2] ; K ++){
        for(size_t J = 0 ; J < nc[1] ; J ++){
            for(size_t I = 0 ; I < nc[0] ; I ++){

                size_t index = I + nc[0] * ( J + nc[1] * K );

                if( index >= material.size() ){
                    printf("Out of bound for material.\n");
                    abort();
                }

                if(    I >= Nx_inf
                    && I < Nx_sup){
                        if(    J >= Ny_inf
                            && J < Ny_sup){
                            if(    K >= Nz_inf
                                && K < Nz_sup){

                                    size_t index_data = 
                                          I-Nx_inf 
                                        + data_new_size[0] * ( J-Ny_inf + (K-Nz_inf) * data_new_size[1] );
                                    if(index_data >= new_data.size()){
                                        printf("Out of bound for data.\n");
                                        printf("Index is %zu, size is %zu.\n",
                                            index_data,new_data.size());
                                        abort();
                                    }
                                    material[index] = new_data[index_data];

                            }
                        }
                    }else{
                        material[index] = 0;
                    }

            }
        }
    }

    // Save material:
    export_spoints_XML("materials", 0, grid, grid, Zip::ZIPPED);

    return EXIT_SUCCESS;
}
