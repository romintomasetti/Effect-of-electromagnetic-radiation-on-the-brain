#include "GridCreator_NEW.h"

#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include "omp.h"

/*
 * Defines for the column number of the properties:
 */
#define COLUMN_PERMEABILITY 4
#define COLUMN_PERMITTIVITY 5
#define COLUMN_ELEC_CONDUC  6
#define COLUMN_MAGN_CONDUC  7

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
    // E_x_eps:
    if(this->E_x_eps != NULL){delete[] this->E_x_eps;}
    // E_x_electrical_cond
    if(this->E_x_electrical_cond != NULL){delete[] this->E_x_electrical_cond;}

    // E_y:
    if(this->E_y != NULL){delete[] this->E_y;}
    // E_y_material:
    if(this->E_y_material != NULL){delete[] this->E_y_material;}
    // E_y_eps:
    if(this->E_y_eps != NULL){delete[] this->E_y_eps;}
    // E_y_electrical_cond
    if(this->E_y_electrical_cond != NULL){delete[] this->E_y_electrical_cond;}

    // E_z:
    if(this->E_z != NULL){delete[] this->E_z;}
    // E_z_material:
    if(this->E_z_material != NULL){delete[] this->E_z_material;}
    // E_z_eps:
    if(this->E_z_eps != NULL){delete[] this->E_z_eps;}
    // E_z_electrical_cond
    if(this->E_z_electrical_cond != NULL){delete[] this->E_z_electrical_cond;}

    // H_x:
    if(this->H_x != NULL){delete[] this->H_x;}
    // H_x_material:
    if(this->H_x_material!= NULL){delete[] this->H_x_material;}
    // H_x_mu:
    if(this->H_x_mu != NULL){delete[] this->H_x_mu;}
    // H_x_magnetic_cond:
    if(this->H_x_magnetic_cond != NULL){delete[] this->H_x_magnetic_cond;}

    // H_y:
    if(this->H_y != NULL){delete[] this->H_y;}
    // H_y_material:
    if(this->H_y_material != NULL){delete[] this->H_y_material;}
    // H_y_mu:
    if(this->H_y_mu != NULL){delete[] this->H_y_mu;}
    // H_y_magnetic_cond:
    if(this->H_y_magnetic_cond != NULL){delete[] this->H_y_magnetic_cond;}

    // H_z:
    if(this->H_z != NULL){delete[] this->H_z;}
    // H_z_material:
    if(this->H_z_material != NULL){delete[] this->H_z_material;}
    // H_z_mu:
    if(this->H_z_mu != NULL){delete[] this->H_z_mu;}
    // H_z_magnetic_cond:
    if(this->H_z_magnetic_cond != NULL){delete[] this->H_z_magnetic_cond;}

    // Temperature:
    if(this->temperature != NULL){delete[] this->temperature;}
    // Temperature_material:
    if(this->temperature_material != NULL){delete[] this->temperature_material;}
    // Thermal conductivity:
    if(this->thermal_conductivity != NULL){delete[] this->thermal_conductivity;}
    // Thermal diffusivity:
    if(this->thermal_diffusivity != NULL){delete[] this->thermal_diffusivity;}


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

    size_t size;

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
    this->size_Ex[0] = M - 1 + 2;
    this->size_Ex[1] = N + 2;
    this->size_Ex[2] = P + 2;
    size = this->size_Ex[0] * this->size_Ex[1] * this->size_Ex[2];
    this->E_x                 = new double[size];
    this->E_x_material        = new unsigned char[size];
    this->E_x_eps             = new double[size];
    this->E_x_electrical_cond = new double [size];

    memory = (8+1+8+8) * size;
    this->profiler.addMemoryUsage("BYTES",memory);

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing E_y" << std::endl;
    #endif
    // Size of E_y is  M × (N − 1) × P. Add 2 nodes in each direction for the neighboors.
    this->size_Ey[0] = M + 2;
    this->size_Ey[1] = N - 1 + 2;
    this->size_Ey[2] = P + 2;
    size = this->size_Ey[0] * this->size_Ey[1] * this->size_Ey[2];
    this->E_y                 = new double[size];
    this->E_y_material        = new unsigned char[size];
    this->E_y_eps             = new double[size];
    this->E_y_electrical_cond = new double[size];

    memory = (8+1+8+8) * size;
    this->profiler.addMemoryUsage("BYTES",memory);

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing E_y" << std::endl;
    #endif
    // Size of E_z is  M × N × (P − 1). Add 2 nodes in each direction for the neighboors.
    this->size_Ez[0] = M + 2;
    this->size_Ez[1] = N + 2;
    this->size_Ez[2] = P - 1 + 2;
    size = this->size_Ez[0] * this->size_Ez[1] * this->size_Ez[2];
    this->E_z                 = new double[size];
    this->E_z_material        = new unsigned char[size];
    this->E_z_eps             = new double[size];
    this->E_z_electrical_cond = new double[size];

    memory = (8+1+8+8) * size;
    this->profiler.addMemoryUsage("BYTES",memory);

    /* ALLOCATE SPACE FOR THE MAGNETIC FIELDS */

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing H_x" << std::endl;
    #endif
    // Size of H_x is  M × (N − 1) × (P − 1). Add 2 nodes in each direction for the neighboors.
    this->size_Hx[0] = M + 2;
    this->size_Hx[1] = N - 1 + 2;
    this->size_Hx[2] = P - 1 + 2;
    size = this->size_Hx[0] * this->size_Hx[1] * this->size_Hx[2];
    this->H_x               = new double[size];
    this->H_x_material      = new unsigned char[size];
    this->H_x_magnetic_cond = new double[size];
    this->H_x_mu            = new double[size];

    memory = (8+1+8+8) * size;
    this->profiler.addMemoryUsage("BYTES",memory);

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing H_y" << std::endl;
    #endif
    // Size of H_y is  (M − 1) × N × (P − 1). Add 2 nodes in each direction for the neighboors.
    this->size_Hy[0] = M - 1 + 2;
    this->size_Hy[1] = N + 2;
    this->size_Hy[2] = P - 1 + 2;
    size = this->size_Hy[0] * this->size_Hy[1] * this->size_Hy[2];
    this->H_y               = new double[size];
    this->H_y_material      = new unsigned char[size];
    this->H_y_mu            = new double[size];
    this->H_y_magnetic_cond = new double[size];

    memory = (8+1+8+8) * size;
    this->profiler.addMemoryUsage("BYTES",memory);

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing H_z" << std::endl;
    #endif
    // Size of H_z is  (M − 1) × (N − 1) × P. Add 2 nodes in each direction for the nieghboors.
    this->size_Hz[0] = M - 1 + 2;
    this->size_Hz[1] = N - 1 + 2;
    this->size_Hz[2] = P + 2;
    size = this->size_Hz[0] * this->size_Hz[1] * this->size_Hz[2];
    this->H_z               = new double[size];
    this->H_z_material      = new unsigned char[size];
    this->H_z_mu            = new double[size];
    this->H_z_magnetic_cond = new double[size];

    memory = (8+1+8+8) * size;
    this->profiler.addMemoryUsage("BYTES",memory);

    /* ALLOCATE SPACE FOR THE TEMPERATURE FIELD */

    #if DEBUG > 2
    std::cout << "GridCreator_New::initializing temperature" << std::endl;
    #endif
    // The temperature grid is homogeneous. (this->size_Thermal^3).
    size_t T = this->size_Thermal * this->size_Thermal * this->size_Thermal;

    if(T == 0){
        fprintf(stderr,"GridCreator_NEW::meshInitialization::ERROR\n");
        fprintf(stderr,"Your temperature grid is empty (size T is 0).\n");
        fprintf(stderr,"Aborting.\nFile %s:%d\n",__FILE__,__LINE__);
        abort();
    }

    this->temperature          = new double[T];
    this->temperature_material = new unsigned char[T];
    this->thermal_conductivity = new double[T];
    this->thermal_diffusivity  = new double[T];

    memory = (8+1+8+8) * T;
    this->profiler.addMemoryUsage("BYTES",memory);

    /* INITIALIZATION OF THE NODES */
    this->Assign_A_Material_To_Each_Node();

    /* INITIALIZATION OF TEMPERATURE NODES (give a initial temperature) */
    this->Assign_Init_Temperature_to_Temperature_nodes();

    // Get elapsed CPU time:
    end___grid_init_CPU_TIME = std::clock();
    double elapsedTimeSec = (end___grid_init_CPU_TIME - start_grid_init_CPU_TIME)
                                 / (double)(CLOCKS_PER_SEC);
    this->profiler.incrementTimingInput("Grid_meshInit",elapsedTimeSec);


}

/* Assign a material to each node (both electromagnetic and thermal grid)*/
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

        #pragma omp parallel num_threads(nbr_omp_threads)
        {
            // EX field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Ex[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Ex[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Ex[0] ; I ++){

                        index = I + this->size_Ex[1] * ( J + this->size_Ex[0] * K );
                        this->E_x_material[index] = mat;
                        
                    }
                }
            }

            // EY field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Ey[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Ey[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Ey[0] ; I ++){

                        index = I + this->size_Ey[1] * ( J + this->size_Ey[0] * K );
                        this->E_y_material[index] = mat;
                        
                    }
                }
            }

            // EZ field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Ez[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Ez[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Ez[0] ; I ++){

                        index = I + this->size_Ez[1] * ( J + this->size_Ez[0] * K );
                        this->E_z_material[index] = mat;
                        
                    }
                }
            }

            // HX field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Hx[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Hx[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Hx[0] ; I ++){

                        index = I + this->size_Hx[1] * ( J + this->size_Hx[0] * K );
                        this->H_x_material[index] = mat;
                        
                    }
                }
            }

            // HY field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Hy[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Hy[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Hy[0] ; I ++){

                        index = I + this->size_Hy[1] * ( J + this->size_Hy[0] * K );
                        this->H_y_material[index] = mat;
                        
                    }
                }
            }

            // HZ field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Hz[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Hz[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Hz[0] ; I ++){

                        index = I + this->size_Hz[1] * ( J + this->size_Hz[0] * K );
                        this->H_z_material[index] = mat;
                        
                    }
                }
            }
                  
            // Temperature nodes:
            #pragma omp for collapse(3)
                    for(size_t I = 0 ; I < T ; I++){
                        for(size_t J = 0 ; J < T ; J ++){
                            for(size_t K = 0 ; K < T ; K ++){
                                this->temperature_material[I + T * (J + T * K)] = mat;
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

// Assign to each node of the temperature field an initial temperature:
void GridCreator_NEW::Assign_Init_Temperature_to_Temperature_nodes(void){

    unsigned int DEFAULT = 4;
    unsigned int nbr_omp_threads = 0;
    if((unsigned)omp_get_num_threads() > DEFAULT){
            nbr_omp_threads = omp_get_num_threads();
    }else{
        nbr_omp_threads = DEFAULT;
    }

    // Retrieve size of the temperature field:
    size_t T = this->size_Thermal;

    if(this->input_parser.get_SimulationType() == "USE_AIR_EVERYWHERE"){

        double init_air_temp = this->input_parser.GetInitTemp_FromMaterialName["AIR"];

        #pragma omp parallel num_threads(nbr_omp_threads)
        {
            #pragma omp for collapse(3)
                    for(size_t I = 0 ; I < T ; I++){
                        for(size_t J = 0 ; J < T ; J ++){
                            for(size_t K = 0 ; K < T ; K ++){
                                this->temperature[I + T * (J + T * K)] = init_air_temp;
                            }
                        }
                    }
        }
    }else{
        fprintf(stderr,"GridCreator_NEW::Assign_Temperature_to_Temperature_nodes::ERROR\n");
        fprintf(stderr,"\t>>> Simulation type doesn't correspond to any known type. Aborting.\n");
        fprintf(stderr,"\tFile %s:%d\n",__FILE__,__LINE__);
        abort();
    }
}

// Assign to each electromagnetic node its properties as a function of the temperature:
void GridCreator_NEW::Initialize_Electromagnetic_Properties(std::string whatToDo /*= string()*/){

    std::cout << "GridCreator_NEW::Initialize_Electromagnetic_Properties::IN" << std::endl;

    unsigned int DEFAULT = 4;
    unsigned int nbr_omp_threads = 0;
    if((unsigned)omp_get_num_threads() > DEFAULT){
            nbr_omp_threads = omp_get_num_threads();
    }else{
        nbr_omp_threads = DEFAULT;
    }

    // Decide what to do based on the argument string 'whatToDo':

    if(whatToDo == "AIR_AT_INIT_TEMP"){
        /*
         * The nodes properties (mu, eps, magn. and elec. cond.) are assigned by assuming initial
         * air temperature.
         */
        size_t M = this->sizes_EH[0];
        size_t N = this->sizes_EH[1];
        size_t P = this->sizes_EH[2];

        size_t index;

        // Retrieve the air temperature:
        double air_init_temp = this->input_parser.GetInitTemp_FromMaterialName["AIR"];
        unsigned char mat    = this->materials.materialID_FromMaterialName["AIR"];
        double eps           = this->materials.getProperty(air_init_temp,
                                                            mat,
                                                            COLUMN_PERMITTIVITY);
        double electric_cond = this->materials.getProperty(air_init_temp,
                                                            mat,
                                                            COLUMN_ELEC_CONDUC);
        double mu            = this->materials.getProperty(air_init_temp,
                                                            mat,
                                                            COLUMN_PERMEABILITY);
        double magnetic_cond = this->materials.getProperty(air_init_temp,
                                                            mat,
                                                            COLUMN_MAGN_CONDUC);

        #pragma omp parallel num_threads(nbr_omp_threads)
        {
            // EX field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Ex[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Ex[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Ex[0] ; I ++){

                        index = I + this->size_Ex[1] * ( J + this->size_Ex[0] * K );
                        this->E_x_eps[index] = eps;
                        this->E_x_electrical_cond[index] = electric_cond;
                        
                    }
                }
            }

            // EY field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Ey[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Ey[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Ey[0] ; I ++){

                        index = I + this->size_Ey[1] * ( J + this->size_Ey[0] * K );
                        this->E_y_eps[index] = eps;
                        this->E_y_electrical_cond[index] = electric_cond;
                        
                    }
                }
            }

            // EZ field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Ez[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Ez[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Ez[0] ; I ++){

                        index = I + this->size_Ez[1] * ( J + this->size_Ez[0] * K );
                        this->E_z_eps[index] = eps;
                        this->E_z_electrical_cond[index] = electric_cond;
                        
                    }
                }
            }

            // HX field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Hx[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Hx[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Hx[0] ; I ++){

                        index = I + this->size_Hx[1] * ( J + this->size_Hx[0] * K );
                        this->H_x_mu[index] = mu;
                        this->H_x_magnetic_cond[index] = magnetic_cond;
                        
                    }
                }
            }

            // HY field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Hy[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Hy[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Hy[0] ; I ++){

                        index = I + this->size_Hy[1] * ( J + this->size_Hy[0] * K );
                        this->H_y_mu[index] = mu;
                        this->H_y_magnetic_cond[index] = magnetic_cond;
                        
                    }
                }
            }

            // HZ field:
            #pragma omp for collapse(3) nowait\
                private(index)
            for(size_t K = 0 ; K < this->size_Hz[2] ; K ++){
                for(size_t J = 0 ; J < this->size_Hz[1] ; J ++ ){
                    for(size_t I = 0 ; I < this->size_Hz[0] ; I ++){

                        index = I + this->size_Hz[1] * ( J + this->size_Hz[0] * K );
                        this->H_z_mu[index] = mu;
                        this->H_z_magnetic_cond[index] = magnetic_cond;
                        
                    }
                }
            }

        }

    }else{
        fprintf(stderr,"GridCreator_NEW::Initialize_Electromagnetic_Properties::ERROR\n");
        fprintf(stderr,"No 'whatToDo' corresponding to %s. Aborting.\n",whatToDo.c_str());
        fprintf(stderr,"In file %s:%d\n",__FILE__,__LINE__);
        abort();
    }

}