#include "GridCreator_NEW.h"

#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include "omp.h"

/* CONSTRUCTOR */
GridCreator_NEW::GridCreator_NEW(InputParser &input_parser,
					    Materials &materials,
					    MPI_Initializer &MPI_communicator,
                        ProfilingClass &profiler):
                        input_parser(input_parser),
                        materials(materials),
                        MPI_communicator(MPI_communicator),
                        profiler(profiler){
    // Retrieve the spatial steps for the electromagnetic and thermal grids:
    this->delta_Electromagn[0] = this->input_parser.deltaX_Electro;
    this->delta_Electromagn[1] = this->input_parser.deltaY_Electro;
    this->delta_Electromagn[2] = this->input_parser.deltaZ_Electro;
    this->delta_Thermal        = this->input_parser.delta_Thermal;

    // Call the MPI division function, from the MPI_communicator field.
    // It retrieves the number of nodes for the electromagnetic grid along each direction:
	this->MPI_communicator.MPI_DIVISION(*this);
}

/* DESTRUCTOR */
GridCreator_NEW::~GridCreator_NEW(void){
    std::cout << "GridCreator_NEW::~GridCreator_NEW::IN" << std::endl;

    /* FREE ALLOCATED SPACE */

    // E_x:
    if(this->E_x != NULL){delete[] this->E_x;}
    // E_x_material:
    if(this->E_x_material !=NULL){delete[] this->E_x_material;}
    // E_y:
    if(this->E_y != NULL){delete[] this->E_y;}
    // E_y_material:
    if(this->E_y_material != NULL){delete[] this->E_y_material;}
    // E_z:
    if(this->E_z != NULL){delete[] this->E_z;}
    // E_z_material:
    if(this->E_z_material != NULL){delete[] this->E_z_material;}
    // H_x:
    if(this->H_x != NULL){delete[] this->H_x;}
    // H_x_material:
    if(this->H_x_material!= NULL){delete[] this->H_x_material;}
    // H_y:
    if(this->H_y != NULL){delete[] this->H_y;}
    // H_y_material:
    if(this->H_y_material != NULL){delete[] this->H_y_material;}
    // H_z:
    if(this->H_z != NULL){delete[] this->H_z;}
    // H_z_material:
    if(this->H_z_material != NULL){delete[] this->H_z_material;}
    // Temperature:
    if(this->temperature != NULL){delete[] this->temperature;}
    // Temperature_material:
    if(this->temperature_material != NULL){delete[] this->temperature_material;}

    std::cout << "GridCreator_NEW::~GridCreator_NEW::OUT" << std::endl;
}

/* GRID INITIALIZATION */
void GridCreator_NEW::meshInitialization(void){
    
    // Timing the grid initialization (in CPU time):
    std::clock_t start_grid_init_CPU_TIME;
    std::clock_t end___grid_init_CPU_TIME;
    this->profiler.addTimingInputToDictionnary("Grid_meshInit");
    start_grid_init_CPU_TIME = std::clock();

    double memory = 0.0;

    size_t M = this->sizes_EH[0];
    size_t N = this->sizes_EH[1];
    size_t P = this->sizes_EH[2];
    if(M == 0 || N == 0 || P == 0){
        fprintf(stderr,"GridCreator_NEW::meshInitialization::ERROR\n");
        fprintf(stderr,"\t>>> One of the quantities (M,N,P)=(%zu,%zu,%zu) is invalid.\nAborting.\n",M,N,P);
        fprintf(stderr,"\t>>> At line %d, in file %s.\n",__LINE__,__FILE__);
        abort();
    }

    /* ALLOCATE SPACE FOR THE ELECTRIC FIELDS */
    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing E_x" << std::endl;
    #endif
    // Size of E_x is  (M − 1) × N × P. Add 2 nodes in each direction for the neighboors.
    this->E_x          = new double[(M-1+2)*(N+2)*(P+2)];
    this->E_x_material = new unsigned char[(M-1+2)*(N+2)*(P+2)];

    memory = (8+1) * (M-1+2)*(N+2)*(P+2);
    this->profiler.addMemoryUsage("BYTES",memory);

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing E_y" << std::endl;
    #endif
    // Size of E_y is  M × (N − 1) × P. Add 2 nodes in each direction for the neighboors.
    this->E_y = new double[(M+2)*(N-1+2)*(P+2)];
    this->E_y_material = new unsigned char[(M+2)*(N-1+2)*(P+2)];

    memory = (8+1) * (M+2)*(N-1+2)*(P+2);
    this->profiler.addMemoryUsage("BYTES",memory);

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing E_y" << std::endl;
    #endif
    // Size of E_z is  M × N × (P − 1). Add 2 nodes in each direction for the neighboors.
    this->E_z          = new double[(M+2)*(N+2)*(P-1+2)];
    this->E_z_material = new unsigned char[(M+2)*(N+2)*(P-1+2)];

    memory = (8+1) * (M+2)*(N+2)*(P-1+2);
    this->profiler.addMemoryUsage("BYTES",memory);

    /* ALLOCATE SPACE FOR THE MAGNETIC FIELDS */

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing H_x" << std::endl;
    #endif
    // Size of H_x is  M × (N − 1) × (P − 1). Add 2 nodes in each direction for the neighboors.
    this->H_x          = new double[(M+2)*(N-1+2)*(P-1+2)];
    this->H_x_material = new unsigned char[(M+2)*(N-1+2)*(P-1+2)];

    memory = (8+1) * (M+2)*(N-1+2)*(P-1+2);
    this->profiler.addMemoryUsage("BYTES",memory);

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing H_y" << std::endl;
    #endif
    // Size of H_y is  (M − 1) × N × (P − 1). Add 2 nodes in each direction for the neighboors.
    this->H_y = new double[(M-1+2)*(N+2)*(P-1+2)];
    this->H_y_material = new unsigned char[(M-1+2)*(N+2)*(P-1+2)];

    memory = (8+1) * (M-1+2)*(N+2)*(P-1+2);
    this->profiler.addMemoryUsage("BYTES",memory);

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing H_z" << std::endl;
    #endif
    // Size of H_z is  (M − 1) × (N − 1) × P. Add 2 nodes in each direction for the nieghboors.
    this->H_z          = new double[(M-1+2)*(N-1+2)*(P+2)];
    this->H_z_material = new unsigned char[(M-1+2)*(N-1+2)*(P+2)];

    memory = (8+1) * (M-1+2)*(N-1+2)*(P+2);
    this->profiler.addMemoryUsage("BYTES",memory);

    /* ALLOCATE SPACE FOR THE TEMPERATURE FIELD */

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing temperature" << std::endl;
    #endif
    // The temperature grid is homogeneous. (this->size_Thermal^3).
    size_t T = this->size_Thermal * this->size_Thermal * this->size_Thermal;
    this->temperature          = new double[T];
    this->temperature_material = new unsigned char[T];

    memory = (8+1) * T;
    this->profiler.addMemoryUsage("BYTES",memory);

    /* INITIALIZATION OF THE NODES */
    this->Assign_A_Material_To_Each_Node();

    end___grid_init_CPU_TIME = std::clock();
    double elapsedTimeSec = (end___grid_init_CPU_TIME - start_grid_init_CPU_TIME)
                                 / (double)(CLOCKS_PER_SEC);
    this->profiler.incrementTimingInput("Grid_meshInit",elapsedTimeSec);
    std::cout << "GridInit => Time: " << elapsedTimeSec << " s" << std::endl;


}

/* Assign a material to each node */
void GridCreator_NEW::Assign_A_Material_To_Each_Node(){
    /*
     * This function fills in the vectors of material.
     */
    // Assign material as a function of the simulation type.

    unsigned int DEFAULT = 4;
    unsigned int nbr_omp_threads = 0;
    if((unsigned)omp_get_num_threads() > DEFAULT){
            nbr_omp_threads = omp_get_num_threads();
    }else{
        nbr_omp_threads = DEFAULT;
    }

    size_t M = this->sizes_EH[0];
    size_t N = this->sizes_EH[1];
    size_t P = this->sizes_EH[2];
    size_t T = this->size_Thermal;
    size_t index;

    if(this->input_parser.get_SimulationType() == "USE_AIR_EVERYWHERE"){
        
        unsigned char mat = this->materials.materialID_FromMaterialName["AIR"];

        double init_air_temp = this->input_parser.GetInitTemp_FromMaterialName["AIR"];

        #pragma omp parallel num_threads(nbr_omp_threads)
        {
            #pragma omp for collapse(3)\
                private(index)
            for(size_t K = 0 ; K < P+2 ; K ++){
                for(size_t J = 0 ; J < N+2 ; J ++ ){
                    for(size_t I = 0 ; I < M+2 ; I ++){

                        index = I + N * ( J + M * K );

                        // Fill in the electric field Ex of size (M − 1) × N × P:
                        if(I < (M-1)+2){
                            this->E_x_material[index] = mat;
                        }
                        // Fill in the electric field Ey of size M × (N − 1) × P:
                        if(J < (N-1)+2){
                            this->E_y_material[index] = mat;
                        }
                        // Fill in the electric field Ez of size M × N × (P − 1):
                        if(K < (P-1)+2){
                            this->E_z_material[index] = mat;
                        }
                        // Fill in the magnetic field Hx of size M × (N − 1) × (P − 1):
                        if(J < (N-1)+2 && K < (P-1)+2){
                            this->H_x_material[index] = mat;
                        }
                        // Fill in the magnetic field Hy of size (M − 1) × N × (P − 1):
                        if(I < (M-1)+2 && K < (P-1)+2){
                            this->H_y_material[index] = mat;
                        }
                        // Fill in the magnetic field Hz of size (M − 1) × (N − 1) × P:
                        if(I < (M-1)+2 && J < (N-1)+2){
                            this->H_z_material[index] = mat;
                        }
                    }
                    /* END OF for(size_t I = 0 ; I < M+2 ; I ++) */
                }
                /* END OF for(size_t J = 0 ; J < N+2 ; J ++ ) */
            }
            /* END OF for(size_t K = 0 ; K < P+2 ; K ++) */

            #pragma omp for collapse(3)\
                private(index)
            for(size_t I = 0 ; I < T ; I++){
                for(size_t J = 0 ; J < T ; J ++){
                    for(size_t K = 0 ; K < T ; K ++){
                        this->temperature[I + T * (J + T * K)] = init_air_temp;
                    }
                }
            }
        }
    }
    /* END OF     if(this->input_parser.get_SimulationType() == "USE_AIR_EVERYWHERE") */
    else if(this->input_parser.get_SimulationType() == "TEST_PARAVIEW"){

    }else{
        fprintf(stderr,"GridCreator_NEW::Assign_A_Material_To_Each_Node()::ERROR\n");
        fprintf(stderr,"\t>>> Simulation type doesn't correspond to any known type. Aborting.\n");
        fprintf(stderr,"\tFile %s:%d\n",__FILE__,__LINE__);
        abort();
    }
}