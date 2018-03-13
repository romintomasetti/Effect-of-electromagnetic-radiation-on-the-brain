#ifndef GRIDCREATOR_NEW_H
#define GRIDCREATOR_NEW_H

#include <vector>
#include <cstdio>

#include "Materials.h"
#include "MPI_Initializer.h"
#include "InputParser.h"
#include "ProfilingClass.h"

#include "vtl.h"
#include "vtlVec3.h"
#include "vtlSPoints.h" 

class MPI_Initializer;
class GridCreator_NEW{
    public:

        const size_t DECALAGE_A_LA_MAIN = 0;

        /* VARIABLES */

        /* ELECTRIC AND MAGNETIC FIELDS - SIZES
         *   Ex of size (M − 1) × N × P
         *   Ey of size M × (N − 1) × P
         *   Ez of size M × N × (P − 1)
         *   Hx of size M × (N − 1) × (P − 1)
         *   Hy of size (M − 1) × N × (P − 1)
         *   Hz of size (M − 1) × (N − 1) × P
         */
        // Spatial steps for electromagnetic fields:
        std::vector<double> delta_Electromagn = {-1.0,-1.0,-1.0};
        // Number of nodes along each direction for the electromagnetic mesh, eqivalent to M,N,P:
        std::vector<size_t> sizes_EH = {0,0,0};
        // Electric fields along X,Y,Z, and the corresponding material 
        // + permittivity (eps) and electrical conductivity:
        double *E_x                 = NULL;
        unsigned char *E_x_material = NULL;
        double *E_x_eps             = NULL;
        double *E_x_electrical_cond = NULL;

        std::vector<size_t> size_Ex = {0,0,0}; //= (M − 1 + 2) * (N + 2) * (P + 2);

        double *E_y                 = NULL;
        unsigned char *E_y_material = NULL;
        double *E_y_eps             = NULL;
        double *E_y_electrical_cond = NULL;

        std::vector<size_t> size_Ey = {0,0,0}; //= (M + 2) * (N − 1 + 2) * (P + 2);

        double *E_z                 = NULL;
        unsigned char *E_z_material = NULL;
        double *E_z_eps             = NULL;
        double *E_z_electrical_cond = NULL;

        std::vector<size_t> size_Ez = {0,0,0}; //= (M + 2) * (N + 2) * (P − 1 + 2);

        // Magnetic fields along X,Y,Z, and the corresponding material 
        // + magnetic permeability (mu) and magnetic conductivity:
        double *H_x                 = NULL;
        unsigned char *H_x_material = NULL;
        double *H_x_mu              = NULL;
        double *H_x_magnetic_cond   = NULL;

        std::vector<size_t> size_Hx = {0,0,0}; //= (M+2) * (N − 1 +2) * (P − 1 +2)

        double *H_y                 = NULL;
        unsigned char *H_y_material = NULL;
        double *H_y_mu              = NULL;
        double *H_y_magnetic_cond   = NULL;

        std::vector<size_t> size_Hy = {0,0,0}; //= (M − 1 +2) * (N+2) * (P − 1 + 2)

        double *H_z                 = NULL;
        unsigned char *H_z_material = NULL;
        double *H_z_mu              = NULL;
        double *H_z_magnetic_cond   = NULL;

        std::vector<size_t> size_Hz = {0,0,0}; //= (M − 1 + 2) * (N − 1 +2) * (P+2)

        /*
         * Spatial step for the thermal grid, considered as homogeneous, i.e. the spatial step is the same
         * in every direction.
         */
        double delta_Thermal = -1.0;
        // Temperature field, and the corresponding material + conductivity and diffusivity:
        double *temperature                 = NULL;
        unsigned char *temperature_material = NULL;
        double *thermal_conductivity        = NULL;
        double *thermal_diffusivity         = NULL;
        std::vector<size_t> size_Thermal = {0,0,0};

        

        

        // Lengths of this grid along each direction:
        std::vector<double> lengths_ofSubgrid    = {-1.0,-1.0,-1.0};

        // Input parser:
        InputParser 	&input_parser;
        // Materials:
		Materials   	&materials;
        // MPI initializer:
		MPI_Initializer &MPI_communicator;
        // Profiler:
        ProfilingClass &profiler;

        // Origin of the indices of the grid:
        std::vector<size_t> originIndices_Electro = {0,0,0};
        std::vector<size_t> originIndices_Thermal = {0,0,0};

        // Origin of the whole simulation, for EM and TH grids:
		vtl::Vec3d originOfWholeSimulation_Electro = vtl::Vec3d(0.0,0.0,0.0);
        vtl::Vec3d originOfWholeSimulation_Thermal = vtl::Vec3d(0.0,0.0,0.0);

        /* FUNCTIONS */

        // Constructor:
        GridCreator_NEW(InputParser &input_parser,
					    Materials &materials,
					    MPI_Initializer &MPI_communicator,
                        ProfilingClass &profiler);
        
        // Destructor:
        ~GridCreator_NEW(void);

        // Grid initialization:
		void meshInitialization(void);

        // Assign a material to each node:
        void Assign_A_Material_To_Each_Node(void);

        // Assign to each node of the temperature field an initial temperature:
        void Assign_Init_Temperature_to_Temperature_nodes(void);

        // Assign to each electromagnetic node its properties as a function of the temperature:
        void Initialize_Electromagnetic_Properties(std::string whatToDo = string());

        // Get the global node number from the local node number, EM grid:
        void get_Global_from_Local_Electro(size_t *local,size_t* global);

        // Get the global node number from the local node number, Thermal grid:
        void get_Global_from_Local_Thermal(size_t *local,size_t* global);

        void Compute_nodes_inside_sources(
                std::vector<size_t>        &local_nodes_inside_source_NUMBER,
                std::vector<unsigned char> &ID_Source,
                std::vector<double>        &local_nodes_inside_source_FREQ,
                size_t                     *nbr_nodes_inside_source,
                const std::string          &TYPE_OF_FIELD
            );
};
#endif