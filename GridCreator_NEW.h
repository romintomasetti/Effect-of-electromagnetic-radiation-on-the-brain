#ifndef GRIDCREATOR_NEW_H
#define GRIDCREATOR_NEW_H

#include <vector>
#include <cstdio>

#include "Materials.h"
#include "MPI_Initializer.h"
#include "InputParser.h"

class GridCreator_NEW{
    public:
        ////////////////////////////////////////////////////////////
        /// VARIABLES:

        // Electric fields along X,Y,Z:
        double *E_x = NULL;
        double *E_y = NULL;
        double *E_z = NULL;
        std::vector<size_t> sizes_E = {0,0,0};

        // Magnetic fields along X,Y,Z:
        double *H_x = NULL;
        double *H_y = NULL;
        double *H_z = NULL;
        std::vector<size_t> sizes_H = {0,0,0};

        // Temperature field:
        double *temperature = NULL;
        std::vector<size_t> sizes_T = {0,0,0};

        // Spatial steps for electromagnetic fields:
        std::vector<double> delta_Electromagn = {-1.0,-1.0,-1.0};

        // Spatial steps for temperature:
        std::vector<double> delta_Temp = {-1.0,-1.0,-1.0};

        // Lengths of the grid along each direction:
        std::vector<double> lengths    = {-1.0,-1.0,-1.0};

        // Input parser:
        InputParser 	&input_parser;
        // Materials:
		Materials   	&materials;
        // MPI initializer:
		MPI_Initializer &MPI_communicator;

        // Origin of the indices:
        std::vector<size_t> originIndices = {0,0,0};

        ////////////////////////////////////////////////////////////
        /// FUNCTIONS:

        // Constructor:
        GridCreator_NEW(InputParser &input_parser,
					    Materials &materials,
					    MPI_Initializer &MPI_communicator);
        
        // Destructor:
        ~GridCreator_NEW(void);
};
#endif