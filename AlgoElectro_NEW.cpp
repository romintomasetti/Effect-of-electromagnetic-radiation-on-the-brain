#include "AlgoElectro_NEW.hpp"

#include <ctime>
#include <new>
#include "omp.h"
#include "mpi.h"
#include <algorithm>
#include <sys/time.h>

#include <assert.h>

#include "header_with_all_defines.hpp"

#define NBR_FACES_CUBE 6

#define DECALAGE_E_SUPP 1





void prepare_array_to_be_sent(
                double ** Electric_field_to_send,
                double ** Magnetic_field_to_send,
                std::vector<size_t> &electric_field_sizes,
                std::vector<size_t> &magnetic_field_sizes,
                double *E_x,
                double *E_y,
                double *E_z,
                double *H_x,
                double *H_y,
                double *H_z,
                int    *mpi_rank_neighboor,
                #ifndef NDEBUG
                    std::vector<size_t> &size_faces_electric,
                    std::vector<size_t> &size_faces_magnetic,
                #endif
                bool is_electric_to_prepare
            );
            
void use_received_array(
                double **Electric_field_to_recv,
                double **Magnetic_field_to_recv,
                std::vector<size_t> &electric_field_sizes,
                std::vector<size_t> &magnetic_field_sizes,
                double *E_x,
                double *E_y,
                double *E_z,
                double *H_x,
                double *H_y,
                double *H_z,
                int    *mpi_rank_neighboor,
                #ifndef NDEBUG
                    std::vector<size_t> &size_faces_electric,
                    std::vector<size_t> &size_faces_magnetic,
                #endif
                bool is_electric_to_use
            );   

void communicate_single_omp_thread(
                double **Electric_field_to_send,
                double **Electric_field_to_recv,
                double **Magnetic_field_to_send,
                double **Magnetic_field_to_recv,
                int *mpi_to_who,
                int  mpi_me,
                std::vector<size_t> size_faces_electric,
                std::vector<size_t> size_faces_magnetic,
                bool is_electric_to_communicate
);

void determine_size_face_based_on_direction(
        char direction,
        std::vector<size_t> &electric_field_sizes,
        std::vector<size_t> &magnetic_field_sizes,
        size_t *size_of_sent_vector_electric,
        size_t *size_of_sent_vector_magnetic){

    /// We don't want to send the neighboors to others, so just do 'minus two' everywhere.
    size_t SEND_ALL = 2;

    /// Reset the sizes of the vectors to be sent to zero:
    *size_of_sent_vector_electric = 0;
    *size_of_sent_vector_magnetic = 0;

    if(direction == 'S' || direction == 'N'){

        /// Compute size of the data to be sent for the electric field:

        /// size += size_ex[1]*size_ex[2] (without neighboors):
        *size_of_sent_vector_electric += (electric_field_sizes[1]-SEND_ALL)    *(electric_field_sizes[2]-SEND_ALL);
        /// size += size_ey[1]*size_ey[2] (without neighboors):
        *size_of_sent_vector_electric += (electric_field_sizes[1+3]-SEND_ALL)  *(electric_field_sizes[2+3]-SEND_ALL);
        /// size += size_ez[1]*size_ez[2] (without neighboors):
        *size_of_sent_vector_electric += (electric_field_sizes[1+2*3]-SEND_ALL)*(electric_field_sizes[2+2*3]-SEND_ALL);

        /// Compute size of the data to be sent for the magnetic field:

        /// size += size_hx[1]*size_hx[2] (without neighboors):
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[1]-SEND_ALL)    *(magnetic_field_sizes[2]-SEND_ALL);
        /// size += size_hy[1]*size_hy[2] (without neighboors):
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[1+3]-SEND_ALL)  *(magnetic_field_sizes[2+3]-SEND_ALL);
        /// size += size_hz[1]*size_hz[2] (without neighboors):
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[1+2*3]-SEND_ALL)*(magnetic_field_sizes[2+2*3]-SEND_ALL);


    }else if(direction == 'W' || direction == 'E'){

        /// Compute size of the data to be sent for the electric field:

        /// size += size_ex[0]*size_ex[2] (without neighboors):
        *size_of_sent_vector_electric += (electric_field_sizes[0]-SEND_ALL) * (electric_field_sizes[2]-SEND_ALL);
        /// size += size_ey[0]*size_ey[2] (without neighboors):
        *size_of_sent_vector_electric += (electric_field_sizes[3]-SEND_ALL) * (electric_field_sizes[5]-SEND_ALL);
        /// size += size_ez[0]*size_ez[2] (without neighboors):
        *size_of_sent_vector_electric += (electric_field_sizes[6]-SEND_ALL) * (electric_field_sizes[8]-SEND_ALL);

        /// Compute size of the data to be sent for the magnetic field:

        /// size += size_hx[0]*size_hx[2] (without neighboors):
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[0]-SEND_ALL)*(magnetic_field_sizes[2]-SEND_ALL);
        /// size += size_hy[0]*size_hy[2] (without neighboors):
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[3]-SEND_ALL)*(magnetic_field_sizes[5]-SEND_ALL);
        /// size += size_z[0]*size_hz[2] (without neighboors):
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[6]-SEND_ALL)*(magnetic_field_sizes[8]-SEND_ALL);

    }else if(direction == 'U' || direction == 'D'){

        /// Compute size of the data to be sent for the electric field:

        /// size += size_ex[0]*size_ex[1]:
        *size_of_sent_vector_electric += (electric_field_sizes[0]-SEND_ALL)*(electric_field_sizes[1]-SEND_ALL);
        /// size += size_ey[0]*size_ey[1]:
        *size_of_sent_vector_electric += (electric_field_sizes[3]-SEND_ALL)*(electric_field_sizes[4]-SEND_ALL);
        /// size += size_ez[0]*size_ez[1]:
        *size_of_sent_vector_electric += (electric_field_sizes[6]-SEND_ALL)*(electric_field_sizes[7]-SEND_ALL);

        /// size += size_ex[0]*size_ex[1]:
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[0]-SEND_ALL)*(magnetic_field_sizes[1]-SEND_ALL);
        /// size += size_ey[0]*size_ey[1]:
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[3]-SEND_ALL)*(magnetic_field_sizes[4]-SEND_ALL);
        /// size += size_ez[0]*size_ez[1]:
        *size_of_sent_vector_magnetic += (magnetic_field_sizes[6]-SEND_ALL)*(magnetic_field_sizes[7]-SEND_ALL);

    }else{
        fprintf(stderr,"In function %s :: unknown direction ! (has %c) -> ABORTING\n",
            __FUNCTION__,direction);
        fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
        #ifdef MPI_COMM_WORLD
        MPI_Abort(MPI_COMM_WORLD,-1);
        #else
        abort();
        #endif
    }

}


/**
 * @brief Compute the smallest time step required for the algorithm statibility.
 */
double AlgoElectro_NEW::Compute_dt(GridCreator_NEW &mesh){
    // Retrieve the spatial step in each direction:
    double dx = mesh.delta_Electromagn[0];
    double dy = mesh.delta_Electromagn[1];
    double dz = mesh.delta_Electromagn[2];
    // Time step:
    double dt = 0.0;
    // Temporary variable:
    double tmp= 0.0;
    // Iterator:
    unsigned char i=0;                                      

    // Iterate over the number of materials. For each material, compute the required time step.
    // At the end, the smallest time step is chosen.                                          
    for (i = 0 ; i < mesh.materials.numberOfMaterials ; i++ ){

            // Get material:
            string material = mesh.materials.materialName_FromMaterialID[i];
            // Get permeability:
            double mu_material = mesh.materials.getProperty(
                    mesh.input_parser.GetInitTemp_FromMaterialName[material],
                    i,4);    

            // Get permittivity:
            double epsilon_material = mesh.materials.getProperty(
                   mesh.input_parser.GetInitTemp_FromMaterialName[material],
                    i,5);     
            // Compute speed of light:
            double c = 1/(sqrt(mu_material*epsilon_material));
            // Take the smallest time step:
            if( i == 0 ){
                dt = 1 / ( c * sqrt( 1 / ( dx * dx ) + 1 / ( dy * dy ) + 1 / (dz *dz) ) );
            }
            else{
                tmp = 1/(c*sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));
                if( tmp < dt ){
                    dt = tmp;
                }
            }
    }
    return dt;
}

/**
 * @brief This is the electromagnetic algorithm (FDTD scheme).
 */
void AlgoElectro_NEW::update(
    GridCreator_NEW &grid,
    InterfaceToParaviewer &interfaceParaview)
{

    /// Start monitoring the time taken for the algorithm to compute:
    struct timeval start, end;
    gettimeofday(&start, NULL);

    grid.profiler.addTimingInputToDictionnary("AlgoElectro_NEW_UPDATE_gettimeofday");

    // Retrieve the time step:
    double dt = this->Compute_dt(grid);
    double current_time = 0.0;
    std::cout << "AlgoElectro_NEW :: dt is " << dt << std::endl;

    /*
     *  TEMPERATURE WILL NEVER CHANGE IN THIS ALGORITHM.
     *  INITIALIZE ALL THE COEFFICIENTS NEEDED FOR THE UPDATE EQUATIONS.
     */
    /* ELECTRIC AND MAGNETIC FIELDS - SIZES
    *   Ex of size (M − 1) × N × P
    *   Ey of size M × (N − 1) × P
    *   Ez of size M × N × (P − 1)
    *   Hx of size M × (N − 1) × (P − 1)
    *   Hy of size (M − 1) × N × (P − 1)
    *   Hz of size (M − 1) × (N − 1) × P
    */

    // In the object grid, set the properties mu, eps, magnetic cond. and electric cond. for each node:
    grid.Initialize_Electromagnetic_Properties("AIR_AT_INIT_TEMP");

    /* Set the coefficients for the electromagnetic update algorithm */

    size_t size;
    
    //size_t M = grid.sizes_EH[0];
    //size_t N = grid.sizes_EH[1];
    //size_t P = grid.sizes_EH[2];

    // Magnetic field Hx:
    size = grid.size_Hx[0] * grid.size_Hx[1] * grid.size_Hx[2];
    double *C_hxh   = new double[size]();
    double *C_hxe_1 = new double[size]();
    double *C_hxe_2 = new double[size]();

    // Magnetic field Hy:
    size = grid.size_Hy[0] * grid.size_Hy[1] * grid.size_Hy[2];
    double *C_hyh   = new double[size]();
    double *C_hye_1 = new double[size]();
    double *C_hye_2 = new double[size]();

    // Magnetic field Hz:
    size = grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2];
    double *C_hzh   = new double[size]();
    double *C_hze_1 = new double[size]();
    double *C_hze_2 = new double[size]();

    // Electric field Ex:
    size = grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2];
    double *C_exe   = new double[size]();
    double *C_exh_1 = new double[size]();
    double *C_exh_2 = new double[size]();

    // Electric field Ey:
    size = grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]; 
    double *C_eye   = new double[size]();
    double *C_eyh_1 = new double[size]();
    double *C_eyh_2 = new double[size]();

    // Electric field Ez:
    size = grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2];
    double *C_eze   = new double[size]();
    double *C_ezh_1 = new double[size]();
    double *C_ezh_2 = new double[size]();

    /* COMPUTING COEFFICIENTS */
    #pragma omp parallel
    {
        size_t index;
        /* Coefficients for Ex */
        // Ex of size (M − 1) × N × P
        #pragma omp for collapse(3) 
            for(size_t K = 0 ; K < grid.size_Ex[2] ; K ++){
                for(size_t J = 0 ; J < grid.size_Ex[1] ; J ++){
                    for(size_t I = 0 ; I < grid.size_Ex[0] ; I ++){

                        index = I + grid.size_Ex[0] * ( J + grid.size_Ex[1] * K);

                        double COEF_E = grid.E_x_electrical_cond[index] * dt
                            / (2.0 * grid.E_x_eps[index]);

                        // Coefficient C_exe:
                        C_exe[index] = (1-COEF_E) / (1+COEF_E);

                        // Coefficient C_exh_1:
                        C_exh_1[index] = 1 / ( 1 + COEF_E) * dt 
                            / (grid.E_x_eps[index] * grid.delta_Electromagn[1]);

                        // Coefficient C_exh_2:
                        C_exh_2[index] = 1 / ( 1 + COEF_E) * dt 
                            / (grid.E_x_eps[index] * grid.delta_Electromagn[2]);

                    }
                }
            }

        /* Coefficients of Ey */
        // Ey is of size M × (N − 1) × P.
        #pragma omp for collapse(3) 
            for(size_t K = 0 ; K < grid.size_Ey[2] ; K ++){
                for(size_t J = 0 ; J < grid.size_Ey[1] ; J ++){
                    for(size_t I = 0 ; I < grid.size_Ey[0] ; I ++){

                        index = I + grid.size_Ey[0] * ( J + grid.size_Ey[1] * K);

                        double COEF_E = grid.E_y_electrical_cond[index] * dt
                            / (2.0 * grid.E_y_eps[index]);

                        // Coefficient C_eye:
                        C_eye[index] = (1-COEF_E) / (1+COEF_E);

                        // Coefficient C_eyh_1:
                        C_eyh_1[index] = 1 / ( 1 + COEF_E) * dt 
                            / (grid.E_y_eps[index] * grid.delta_Electromagn[2]);

                        // Coefficient C_eyh_2:
                        C_eyh_2[index] = 1 / ( 1 + COEF_E) * dt 
                            / (grid.E_y_eps[index] * grid.delta_Electromagn[0]);
                    }
                }
            }

        /* Coefficients of Ez */
        // Ez is of size  M × N × (P − 1)
        #pragma omp for collapse(3) 
            for(size_t K = 0 ; K < grid.size_Ez[2] ; K ++){
                for(size_t J = 0 ; J < grid.size_Ez[1] ; J ++){
                    for(size_t I = 0 ; I < grid.size_Ez[0] ; I ++){

                        index = I + grid.size_Ez[0] * ( J + grid.size_Ez[1] * K);

                        double COEF_E = grid.E_z_electrical_cond[index] * dt
                            / (2.0 * grid.E_z_eps[index]);

                        // Coefficient C_eze:
                        C_eze[index] = (1-COEF_E) / (1+COEF_E);

                        // Coefficient C_ezh_1:
                        C_ezh_1[index] = 1 / ( 1 + COEF_E) * dt 
                            / (grid.E_z_eps[index] * grid.delta_Electromagn[0]);

                        // Coefficient C_ezh_2:
                        C_ezh_2[index] = 1 / ( 1 + COEF_E) * dt 
                            / (grid.E_z_eps[index] * grid.delta_Electromagn[1]);
                    }
                }
            }

        /* Coefficients of Hx */
        // Hx is of size M × (N − 1) × (P − 1):
        #pragma omp for collapse(3) 
            for(size_t K = 0 ; K < grid.size_Hx[2] ; K ++){
                for(size_t J = 0 ; J < grid.size_Hx[1] ; J ++){
                    for(size_t I = 0 ; I < grid.size_Hx[0] ; I ++){

                        index = I + grid.size_Hx[0] * ( J + grid.size_Hx[1] * K);

                        double COEF_H = grid.H_x_magnetic_cond[index] * dt
                            / (2.0 * grid.H_x_mu[index]);

                        // Coefficient C_hxh:
                        C_hxh[index] = (1-COEF_H) / (1+COEF_H);

                        // Coefficient C_hxe_1:
                        C_hxe_1[index] = 1 / ( 1 + COEF_H) * dt 
                            / (grid.H_x_mu[index] * grid.delta_Electromagn[2]);

                        // Coefficient C_hxe_2:
                        C_hxe_2[index] = 1 / ( 1 + COEF_H) * dt 
                            / (grid.H_x_mu[index] * grid.delta_Electromagn[1]);
                    }
                }
            }

        /* Coefficients of Hy */
        // Hy is of size (M − 1) × N × (P − 1)
        #pragma omp for collapse(3) 
            for(size_t K = 0 ; K < grid.size_Hy[2] ; K ++){
                for(size_t J = 0 ; J < grid.size_Hy[1] ; J ++){
                    for(size_t I = 0 ; I < grid.size_Hy[0] ; I ++){

                        index = I + grid.size_Hy[0] * ( J + grid.size_Hy[1] * K);

                        double COEF_H = grid.H_y_magnetic_cond[index] * dt
                            / (2.0 * grid.H_y_mu[index]);

                        // Coefficient C_hxh:
                        C_hyh[index] = (1-COEF_H) / (1+COEF_H);

                        // Coefficient C_hxe_1:
                        C_hye_1[index] = 1 / ( 1 + COEF_H) * dt 
                            / (grid.H_y_mu[index] * grid.delta_Electromagn[0]);

                        // Coefficient C_hxe_2:
                        C_hye_2[index] = 1 / ( 1 + COEF_H) * dt 
                            / (grid.H_y_mu[index] * grid.delta_Electromagn[2]);
                    }
                }
            }

        /* Coefficients of Hz */
        // Hz is of size (M − 1) × (N − 1) × P
        #pragma omp for collapse(3) 
            for(size_t K = 0 ; K < grid.size_Hz[2] ; K ++){
                for(size_t J = 0 ; J < grid.size_Hz[1] ; J ++){
                    for(size_t I = 0 ; I < grid.size_Hz[0] ; I ++){

                        index = I + grid.size_Hz[0] * ( J + grid.size_Hz[1] * K);

                        double COEF_H = grid.H_z_magnetic_cond[index] * dt
                            / (2.0 * grid.H_z_mu[index]);

                        // Coefficient C_hxh:
                        C_hzh[index] = (1-COEF_H) / (1+COEF_H);

                        // Coefficient C_hxe_1:
                        C_hze_1[index] = 1 / ( 1 + COEF_H) * dt 
                            / (grid.H_z_mu[index] * grid.delta_Electromagn[1]);

                        // Coefficient C_hxe_2:
                        C_hze_2[index] = 1 / ( 1 + COEF_H) * dt 
                            / (grid.H_z_mu[index] * grid.delta_Electromagn[0]);
                    }
                }
            }

    }
    // END OF #pragma omp parallel

    /////////////////////////////////////////////////////
    // COMPUTE NODES INSIDE THE SOURCE:                //
    //      1) Call GridCreator_NEW                    //
    //      2) GridCreator_NEW calls its source object //
    /////////////////////////////////////////////////////
    
    // Size 6 because 3 E and 3 H components:
    std::vector<size_t>        *local_nodes_inside_source_NUMBER ;
    local_nodes_inside_source_NUMBER = new std::vector<size_t>[6];

    // Size 6 because 3 E and 3 H components:
    std::vector<unsigned char> *ID_Source                         ;
    ID_Source = new std::vector<unsigned char>[6];

    std::vector<double>        local_nodes_inside_source_FREQ    ;

    std::vector<std::string> TYPE = {"Ex","Ey","Ez","Hx","Hy","Hz"};

    for(unsigned int i = 0 ; i < TYPE.size() ; i ++){
        grid.Compute_nodes_inside_sources(
            local_nodes_inside_source_NUMBER[i],
            ID_Source[i],
            local_nodes_inside_source_FREQ,
            TYPE[i]
        );
        if(local_nodes_inside_source_NUMBER[i].size() != ID_Source[i].size()){
            fprintf(stderr,"In function %s :: wrong sizes !\n",__FUNCTION__);
            fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
            #ifdef MPI_COMM_WORLD
            MPI_Abort(MPI_COMM_WORLD,-1);
            #else
            abort();
            #endif
        }
    }

    /// Assign frequencies:
    local_nodes_inside_source_FREQ = grid.input_parser.source.frequency;
    
    /// Clean the output:
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    if(omp_get_thread_num() == 0 
        && grid.MPI_communicator.isRootProcess() != INT_MIN)
        {
            printf(">>> FDTD scheme started with time step of %lf seconds.\n",
                dt);
        }


    ///////////////////////////////////////////////
    // UPDATE WHILE LOOP - PARALLELIZED WITH     //
    // OPENMP THREADS    - MINIMUM 6 OPENMP      //
    // THREADS ARE REQUIRED FOR MPI COMMUICATION //
    // TO WORK.                                  //
    ///////////////////////////////////////////////

    /// Set OMP_DYNAMIC=false:
    this->check_OMP_DYNAMIC_envVar();

    /**
     * The following arrays contain the electric field nodes to be sent and received.
     * 6 is for the number of faces.
     */
    double **Electric_field_to_send = (double**) calloc(NBR_FACES_CUBE,sizeof(double*));
    double **Electric_field_to_recv = (double**) calloc(NBR_FACES_CUBE,sizeof(double*));

    double **Magnetic_field_to_send = (double**) calloc(NBR_FACES_CUBE,sizeof(double*));
    double **Magnetic_field_to_recv = (double**) calloc(NBR_FACES_CUBE,sizeof(double*));

    /**
     * Initialize Electric_field_to_send/recv when it is appropriate:
     */
    std::vector<char> DIRECTIONS = {'S','N','W','E','D','U'};

    std::vector<size_t> electric_field_sizes = {
            grid.size_Ex[0],
            grid.size_Ex[1],
            grid.size_Ex[2],
            grid.size_Ey[0],
            grid.size_Ey[1],
            grid.size_Ey[2],
            grid.size_Ez[0],
            grid.size_Ez[1],
            grid.size_Ez[2]
        };
    std::vector<size_t> magnetic_field_sizes = {
            grid.size_Hx[0],
            grid.size_Hx[1],
            grid.size_Hx[2],
            grid.size_Hy[0],
            grid.size_Hy[1],
            grid.size_Hy[2],
            grid.size_Hz[0],
            grid.size_Hz[1],
            grid.size_Hz[2]
        };

    std::vector<size_t> size_faces_electric(6);
    std::vector<size_t> size_faces_magnetic(6);

    for(unsigned int i = 0 ; i < NBR_FACES_CUBE ; i ++){
        if(grid.MPI_communicator.RankNeighbour[i] != -1){
            char dir = DIRECTIONS[i];

            determine_size_face_based_on_direction(
                                dir,
                                electric_field_sizes,
                                magnetic_field_sizes,
                                &size_faces_electric[i],
                                &size_faces_magnetic[i]);

            Electric_field_to_send[i] = (double*) calloc(size_faces_electric[i],sizeof(double));
            Electric_field_to_recv[i] = (double*) calloc(size_faces_electric[i],sizeof(double));
            Magnetic_field_to_send[i] = (double*) calloc(size_faces_magnetic[i],sizeof(double));
            Magnetic_field_to_recv[i] = (double*) calloc(size_faces_magnetic[i],sizeof(double));
        }
    }


    ////////////////////////////////////
    /// BEGINNING OF PARALLEL REGION ///
    ////////////////////////////////////
    /**
     * Here is a list of shared and private variables in he following parallel region:
     *      1) grid         : instance of the class GridCreator_NEW;
     *      2) dt           : time step of the electromagnetic algorithm;
     *      3) current_time : current simulation time, including thermal solver;
     *      4) currentStep  : current step of the electromagnetic solver;
     *      5) local_nodes_inside_source_NUMBER : index of all nodes inside the sources;
     *      6) interfaceParaviewer : to write the output..
     * 
     *      TO DO
     */
    #pragma omp parallel num_threads(omp_get_max_threads()) default(none)\
        shared(grid,current_time,end,start)\
        shared(local_nodes_inside_source_NUMBER)\
        shared(local_nodes_inside_source_FREQ,ID_Source)\
        shared(interfaceParaview)\
        shared(C_hxh,C_hxe_1,C_hxe_2)\
        shared(C_hyh,C_hye_1,C_hye_2)\
        shared(C_hzh,C_hze_1,C_hze_2)\
        shared(C_exe,C_exh_1,C_exh_2)\
        shared(C_eye,C_eyh_1,C_eyh_2)\
        shared(C_eze,C_ezh_1,C_ezh_2)\
        shared(ompi_mpi_comm_world,ompi_mpi_int)\
        shared(Electric_field_to_send,Electric_field_to_recv)\
        shared(Magnetic_field_to_send,Magnetic_field_to_recv)\
        shared(electric_field_sizes,magnetic_field_sizes,dt)\
        shared(size_faces_electric,size_faces_magnetic)
    {

        /*
        shared(H_x_tmp,H_y_tmp,H_z_tmp)\
        shared(E_x_tmp,E_y_tmp,E_z_tmp)\
        */

        // Temporary pointers, to avoid doing grid.sthg !
        double *H_x_tmp = grid.H_x;
        double *H_y_tmp = grid.H_y;
        double *H_z_tmp = grid.H_z;
        double *E_x_tmp = grid.E_x;
        double *E_y_tmp = grid.E_y;
        double *E_z_tmp = grid.E_z;

        size_t index;
        size_t index_1Plus;
        size_t index_1Moins;
        size_t index_2Plus;
        size_t index_2Moins;
        size_t size_x;
        size_t size_y;
        size_t size_x_1;
        size_t size_y_1;
        size_t size_x_2;
        size_t size_y_2;

        size_t currentStep = 0;

        /**
         * Important for the electric field update !
         */
        size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX = 0;
        size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY = 0;
        size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ = 0;

        size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX = 0;
        size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY = 0;
        size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ = 0;

        if(grid.MPI_communicator.MPI_POSITION[0] == grid.MPI_communicator.MPI_MAX_POSI[0]){
            IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX = 1;
        }

        if(grid.MPI_communicator.MPI_POSITION[1] == grid.MPI_communicator.MPI_MAX_POSI[1]){
            IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY = 1;
        }

        if(grid.MPI_communicator.MPI_POSITION[2] == grid.MPI_communicator.MPI_MAX_POSI[2]){
            IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ = 1;
        }

        if(grid.MPI_communicator.MPI_POSITION[0] == 0){
            IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX = 1;
        }

        if(grid.MPI_communicator.MPI_POSITION[1] == 0){
            IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY = 1;
        }

        if(grid.MPI_communicator.MPI_POSITION[2] == 0){
            IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ = 1;
        }

        // Some indexing variables:
        size_t I,J,K;

        /// Variables to monitore the time spent communicating:
        struct timeval start_mpi_comm;
        struct timeval end___mpi_comm;
        double         total_mpi_comm = 0.0;

        /// Variables to compute the time taken by each iteration:
        struct timeval start_while_iter;
        struct timeval end___while_iter;
        double         total_while_iter = 0.0;

        while(current_time < grid.input_parser.get_stopTime()
                && currentStep < grid.input_parser.maxStepsForOneCycleOfElectro){

            gettimeofday( &start_while_iter , NULL);

            // Updating the magnetic field Hx.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Hx[0];
            size_y = grid.size_Hx[1];

            size_x_1 = grid.size_Ey[0];
            size_y_1 = grid.size_Ey[1];

            size_x_2 = grid.size_Ez[0];
            size_y_2 = grid.size_Ez[1];

            #ifndef NDEBUG
                #pragma omp master
                printf("%s>>> %s!!! WARNING !!!%s 'NDEBUG' is not defied. You are in debug mode."
                       " Be aware that the code is subsequently much slower. As an example,"
                       " a lot of asserts and printf's are performed.%s\n",
                        ANSI_COLOR_RED,
                        ANSI_COLOR_YELLOW,
                        ANSI_COLOR_GREEN,
                        ANSI_COLOR_RESET);
            #endif

            #pragma omp for schedule(static) collapse(3) 
            for (K = 1 ; K < grid.size_Hx[2]-1 ; K++){
                for(J = 1 ; J < grid.size_Hx[1]-1 ; J ++){
                    for(I = 1 ; I < grid.size_Hx[0]-1 ; I++){

                        // Hx(mm, nn, pp):
                        index        = I + size_x   * ( J     + size_y   * K);
                        // Ey(mm, nn, pp + 1):
                        index_1Plus  = I + size_x_1 * ( J     + size_y_1 * (K+1));
                        // Ey(mm, nn, pp):
                        index_1Moins = I + size_x_1 * ( J     + size_y_1 * K);
                        // Ez(mm, nn + 1, pp):
                        index_2Plus  = I + size_x_2 * ( (J+1) + size_y_2 * K);
                        // Ez(mm, nn, pp):
                        index_2Moins = I + size_x_2 * ( J     + size_y_2 * K);
                        

                        ASSERT(index,<,grid.size_Hx[0]*grid.size_Hx[1]*grid.size_Hx[2]);
                        ASSERT(index_1Plus,<,grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]);
                        ASSERT(index_1Moins,<,grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]);
                        ASSERT(index_2Plus,<,grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2]);
                        ASSERT(index_2Moins,<,grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2]);

                        H_x_tmp[index] = C_hxh[index] * H_x_tmp[index]
                                + C_hxe_1[index] * (E_y_tmp[index_1Plus] - E_y_tmp[index_1Moins])
                                - C_hxe_2[index] * (E_z_tmp[index_2Plus] - E_z_tmp[index_2Moins]);

                    }
                }
            }

            // Updating the magnetic field Hy.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Hy[0];
            size_y = grid.size_Hy[1];

            size_x_1 = grid.size_Ez[0];
            size_y_1 = grid.size_Ez[1];

            size_x_2 = grid.size_Ex[0];
            size_y_2 = grid.size_Ex[1];

            #pragma omp for schedule(static) collapse(3) 
            for(K = 1 ; K < grid.size_Hy[2]-1 ; K ++){
                for(J = 1 ; J < grid.size_Hy[1]-1 ; J ++){
                    for(I = 1; I < grid.size_Hy[0]-1 ; I ++){

                        index        = I   + size_x   * ( J  + size_y   * K);
                        // Ez(mm + 1, nn, pp):
                        index_1Plus  = I+1 + size_x_1 * ( J  + size_y_1 * K);
                        // Ez(mm, nn, pp)
                        index_1Moins = I   + size_x_1 * ( J  + size_y_1 * K);
                        // Ex(mm, nn, pp + 1):
                        index_2Plus  = I   + size_x_2 * ( J  + size_y_2 * (K+1));
                        // Ex(mm, nn, pp):
                        index_2Moins = I   + size_x_2 * ( J  + size_y_2 * K);

                        ASSERT(index,<,grid.size_Hy[0]*grid.size_Hy[1]*grid.size_Hy[2]);
                        ASSERT(index_2Plus,<,grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2]);
                        ASSERT(index_2Moins,<,grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2]);
                        ASSERT(index_1Plus,<,grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2]);
                        ASSERT(index_1Moins,<,grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2]);

                        H_y_tmp[index] = C_hyh[index] * H_y_tmp[index]
                                + C_hye_1[index] * (E_z_tmp[index_1Plus] - E_z_tmp[index_1Moins])
                                - C_hye_2[index] * (E_x_tmp[index_2Plus] - E_x_tmp[index_2Moins]);

                    }
                }
            }

            // Updating the magnetic field Hz.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Hz[0];
            size_y = grid.size_Hz[1];

            size_x_1 = grid.size_Ex[0];
            size_y_1 = grid.size_Ex[1];

            size_x_2 = grid.size_Ey[0];
            size_y_2 = grid.size_Ey[1];

            #pragma omp for schedule(static) collapse(3) 
            for(K = 1 ; K < grid.size_Hz[2]-1  ; K ++){
                for(J = 1 ; J < grid.size_Hz[1]-1 ; J ++){
                    for(I = 1 ; I < grid.size_Hz[0]-1 ; I ++){

                        index        = I   + size_x   * ( J     + size_y   * K);
                        // Ex(mm, nn + 1, pp)
                        index_1Plus  = I   + size_x_1 * ( (J+1) + size_y_1 * K);
                        // Ex(mm, nn, pp)
                        index_1Moins = I   + size_x_1 * ( J     + size_y_1 * K);
                        // Ey(mm + 1, nn, pp)
                        index_2Plus  = I+1 + size_x_2 * ( J     + size_y_2 * K);
                        // Ey(mm, nn, pp)
                        index_2Moins = I   + size_x_2 * ( J     + size_y_2 * K);

                        ASSERT(index,<,grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2]);
                        ASSERT(index_1Plus,<,grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2]);
                        ASSERT(index_1Moins,<,grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2]);
                        ASSERT(index_2Plus,<,grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]);
                        ASSERT(index_2Moins,<,grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]);

                        H_z_tmp[index] = C_hzh[index] * H_z_tmp[index]
                                + C_hze_1[index] * (E_x_tmp[index_1Plus] - E_x_tmp[index_1Moins])
                                - C_hze_2[index] * (E_y_tmp[index_2Plus] - E_y_tmp[index_2Moins]);

                    }
                }
            }

            /////////////////////////////////////////////////////
            /// OPENMP barrier because we must ensure all the ///
            /// magnetic fields have been updated.            ///
            /////////////////////////////////////////////////////
            #pragma omp barrier


            /////////////////////////
            /// MPI COMMUNICATION ///
            /////////////////////////
            gettimeofday( &start_mpi_comm, NULL);
            /// Prepare the array to send:
            prepare_array_to_be_sent(
                Electric_field_to_send,
                Magnetic_field_to_send,
                electric_field_sizes,
                magnetic_field_sizes,
                E_x_tmp,
                E_y_tmp,
                E_z_tmp,
                H_x_tmp,
                H_y_tmp,
                H_z_tmp,
                grid.MPI_communicator.RankNeighbour,
                #ifndef NDEBUG
                    size_faces_electric,
                    size_faces_magnetic,
                #endif
                false /* false : tells the function we want to deal with magnetic field only */
            );

            /// Wait all omp threads:
            #pragma omp barrier

            /// Only the master thread communicates:
            #pragma omp master
            {
                communicate_single_omp_thread(
                    Electric_field_to_send,
                    Electric_field_to_recv,
                    Magnetic_field_to_send,
                    Magnetic_field_to_recv,
                    grid.MPI_communicator.RankNeighbour,
                    grid.MPI_communicator.getRank(),
                    size_faces_electric,
                    size_faces_magnetic,
                    false /* false : tells the function we want to deal with magnetic field only */
                );
            }

            /// Other threads wait for the communication to be done:
            #pragma omp barrier

            /// Fill in the matrix of electric field with what was received:
            use_received_array(
                Electric_field_to_recv,
                Magnetic_field_to_recv,
                electric_field_sizes,
                magnetic_field_sizes,
                E_x_tmp,
                E_y_tmp,
                E_z_tmp,
                H_x_tmp,
                H_y_tmp,
                H_z_tmp,
                grid.MPI_communicator.RankNeighbour,
                #ifndef NDEBUG
                    size_faces_electric,
                    size_faces_magnetic,
                #endif
                false /* false : tells the function we want to deal with magnetic field only */
            );           

            #pragma omp barrier
            gettimeofday( &end___mpi_comm , NULL);
            total_mpi_comm += end___mpi_comm.tv_sec  - start_mpi_comm.tv_sec + 
                                (end___mpi_comm.tv_usec - start_mpi_comm.tv_usec) / 1.e6;
            /////////////////////////
            ///      END OF       ///
            /// MPI COMMUNICATION ///
            /////////////////////////



            // Updating the electric field Ex.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Ex[0];
            size_y = grid.size_Ex[1];

            size_x_1 = grid.size_Hz[0];
            size_y_1 = grid.size_Hz[1];

            size_x_2 = grid.size_Hy[0];
            size_y_2 = grid.size_Hy[1];


            #pragma omp for schedule(static) collapse(3) 
            for(K = 1 + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ ;
                    K < grid.size_Ex[2]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ ; 
                    K ++){
                for(J = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY ; 
                        J < grid.size_Ex[1]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY ; 
                        J ++){
                    for(I = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX; 
                            I < grid.size_Ex[0]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX; 
                            I ++){

                        index        = I   + size_x   * ( J     + size_y   * K);
                        // Hz(mm, nn, pp)
                        index_1Plus  = I
                                   + size_x_1 * ( J
                                        + size_y_1 * (K));
                        // Hz(mm, nn - 1, pp)
                        index_1Moins = I
                                   + size_x_1 * ( J-1
                                      + size_y_1 * (K));
                        // Hy(mm, nn, pp)
                        index_2Plus  = I 
                                   + size_x_2 * ( J
                                        + size_y_2 * (K));
                        // Hy(mm, nn, pp - 1)
                        index_2Moins = I 
                                   + size_x_2 * ( J 
                                        + size_y_2 * (K-1 ));

                        ASSERT(index,<,grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2]);
                        ASSERT(index_1Plus,<,grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2]);
                        ASSERT(index_1Moins,<,grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2]);
                        ASSERT(index_2Plus,<,grid.size_Hy[0]*grid.size_Hy[1]*grid.size_Hy[2]);
                        ASSERT(index_2Moins,<,grid.size_Hy[0]*grid.size_Hy[1]*grid.size_Hy[2]);

                        E_x_tmp[index] = C_exe[index] * E_x_tmp[index]
                                + C_exh_1[index] * (H_z_tmp[index_1Plus] - H_z_tmp[index_1Moins])
                                - C_exh_2[index] * (H_y_tmp[index_2Plus] - H_y_tmp[index_2Moins]);

                    }
                }
            }

            // Updating the electric field Ey.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Ey[0];
            size_y = grid.size_Ey[1];

            size_x_1 = grid.size_Hx[0];
            size_y_1 = grid.size_Hx[1];

            size_x_2 = grid.size_Hz[0];
            size_y_2 = grid.size_Hz[1];

            #pragma omp for schedule(static) collapse(3) 
            for(K = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ ; K < grid.size_Ey[2]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ ; K ++){
                for(J = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY; J < grid.size_Ey[1]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY ; J ++){
                    for(I = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX ; I < grid.size_Ey[0]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX ; I ++){

                        index        = I   + size_x   * ( J     + size_y   * K);
                        // Hx(mm, nn, pp)
                        index_1Plus  = I 
                               + size_x_1 * ( J 
                                    + size_y_1 * (K ));
                        // Hx(mm, nn, pp - 1)
                        index_1Moins = I 
                                    + size_x_1 * ( J 
                                         + size_y_1 * (K-1 ));
                        // Hz(mm, nn, pp)
                        index_2Plus  = I
                                    + size_x_2 * ( J 
                                         + size_y_2 * (K ));
                        // Hz(mm - 1, nn, pp)
                        index_2Moins = I-1 
                                    + size_x_2 * ( J 
                                         + size_y_2 * (K ));

                        ASSERT(index,<,grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]);
                        ASSERT(index_2Plus,<,grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2]);
                        ASSERT(index_2Moins,<,grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2]);
                        ASSERT(index_1Plus,<,grid.size_Hx[0]*grid.size_Hx[1]*grid.size_Hx[2]);
                        ASSERT(index_1Moins,<,grid.size_Hx[0]*grid.size_Hx[1]*grid.size_Hx[2]);

                        E_y_tmp[index] = C_eye[index] * E_y_tmp[index]
                                + C_eyh_1[index] * (H_x_tmp[index_1Plus] - H_x_tmp[index_1Moins])
                                - C_eyh_2[index] * (H_z_tmp[index_2Plus] - H_z_tmp[index_2Moins]);

                    }
                }
            }

            // Updating the electric field Ez.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Ez[0];
            size_y = grid.size_Ez[1];

            size_x_1 = grid.size_Hy[0];
            size_y_1 = grid.size_Hy[1];

            size_x_2 = grid.size_Hx[0];
            size_y_2 = grid.size_Hx[1];

            #pragma omp for schedule(static) collapse(3) 
            for(K = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ; K < grid.size_Ez[2]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ ; K ++){
                for(J = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY ; J < grid.size_Ez[1]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY ; J ++){
                    for(I = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX ; I < grid.size_Ez[0]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX ; I ++){

                        index        = I   + size_x   * ( J     + size_y   * K);
                        // Hy(mm, nn, pp)
                        index_1Plus  = I 
                                + size_x_1 * ( J 
                                        + size_y_1 * (K ));
                        // Hy(mm - 1, nn, pp)
                        index_1Moins = I-1 
                                + size_x_1 * ( J 
                                        + size_y_1 * (K ));
                        // Hx(mm, nn, pp)
                        index_2Plus  = I 
                                + size_x_2 * ( J 
                                        + size_y_2 * (K));
                        // Hx(mm, nn - 1, pp)
                        index_2Moins = I 
                                + size_x_2 * ( J-1 
                                        + size_y_2 * (K ));

                        ASSERT(index,<,grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2]);
                        ASSERT(index_1Plus,<,grid.size_Hy[0]*grid.size_Hy[1]*grid.size_Hy[2]);
                        ASSERT(index_1Moins,<,grid.size_Hy[0]*grid.size_Hy[1]*grid.size_Hy[2]);
                        ASSERT(index_2Plus,<,grid.size_Hx[0]*grid.size_Hx[1]*grid.size_Hx[2]);
                        ASSERT(index_2Moins,<,grid.size_Hx[0]*grid.size_Hx[1]*grid.size_Hx[2]);

                        E_z_tmp[index] = C_eze[index] * E_z_tmp[index]
                                + C_ezh_1[index] * (H_y_tmp[index_1Plus] - H_y_tmp[index_1Moins])
                                - C_ezh_2[index] * (H_x_tmp[index_2Plus] - H_x_tmp[index_2Moins]);
                    }
                }
            }
            #pragma omp barrier

            ////////////////////////////
            /// IMPOSING THE SOURCES ///
            ////////////////////////////
            #pragma omp for schedule(static)
            for(size_t it = 0 ; it < local_nodes_inside_source_NUMBER[2].size() ; it ++){

                index = local_nodes_inside_source_NUMBER[2][it];
                ASSERT(index,<,grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2]);

                double frequency = local_nodes_inside_source_FREQ[ID_Source[2][it]];

                E_z_tmp[index] = sin(2*M_PI*frequency*current_time);
            }

            #pragma omp for schedule(static)
            for(size_t it = 0 ; it < local_nodes_inside_source_NUMBER[1].size() ; it ++){

                index = local_nodes_inside_source_NUMBER[1][it];
                ASSERT(index,<,grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]);

                E_y_tmp[index] = 0;
            }

            #pragma omp for schedule(static)
            for(size_t it = 0 ; it < local_nodes_inside_source_NUMBER[0].size() ; it ++){

                index = local_nodes_inside_source_NUMBER[0][it];
                ASSERT(index,<,grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2]);

                E_x_tmp[index] = 0;
            }

            
            
            /////////////////////////
            /// MPI COMMUNICATION ///
            /////////////////////////

            /// Wait all OPENMP threads to be sure computations are done for this step:
            #pragma omp barrier
            gettimeofday( &start_mpi_comm , NULL);
            /// Prepare the array to send:
            prepare_array_to_be_sent(
                Electric_field_to_send,
                Magnetic_field_to_send,
                electric_field_sizes,
                magnetic_field_sizes,
                E_x_tmp,
                E_y_tmp,
                E_z_tmp,
                H_x_tmp,
                H_y_tmp,
                H_z_tmp,
                grid.MPI_communicator.RankNeighbour,
                #ifndef NDEBUG
                    size_faces_electric,
                    size_faces_magnetic,
                #endif
                true // True : telling the function we communicate only electric field
            );

            /// Wait all omp threads:
            #pragma omp barrier

            /// Only the master thread communicates:
            #pragma omp master
            {
                communicate_single_omp_thread(
                    Electric_field_to_send,
                    Electric_field_to_recv,
                    Magnetic_field_to_send,
                    Magnetic_field_to_recv,
                    grid.MPI_communicator.RankNeighbour,
                    grid.MPI_communicator.getRank(),
                    size_faces_electric,
                    size_faces_magnetic,
                    true
                );
            }

            /// Other threads wait for the communication to be done:
            #pragma omp barrier

            /// Fill in the matrix of electric field with what was received:
            use_received_array(
                Electric_field_to_recv,
                Magnetic_field_to_recv,
                electric_field_sizes,
                magnetic_field_sizes,
                E_x_tmp,
                E_y_tmp,
                E_z_tmp,
                H_x_tmp,
                H_y_tmp,
                H_z_tmp,
                grid.MPI_communicator.RankNeighbour,
                #ifndef NDEBUG
                    size_faces_electric,
                    size_faces_magnetic,
                #endif
                true // True : telling the function that we communicate electric field only.
            );           

            #pragma omp barrier
            gettimeofday( &end___mpi_comm , NULL);
            total_mpi_comm += end___mpi_comm.tv_sec  - start_mpi_comm.tv_sec + 
                                (end___mpi_comm.tv_usec - start_mpi_comm.tv_usec) / 1.e6;

            /////////////////////////
            ///      END OF       ///
            /// MPI COMMUNICATION ///
            /////////////////////////
            
            gettimeofday( &end___while_iter , NULL);
            total_while_iter += end___while_iter.tv_sec  - start_while_iter.tv_sec + 
                                (end___while_iter.tv_usec - start_while_iter.tv_usec) / 1.e6;

            currentStep ++;

            #pragma omp master
            {
                /// If this is the first step, add some inputs to the profiler:
                if(currentStep == 1){
                    grid.profiler.addTimingInputToDictionnary("ELECTRO_WRITING_OUTPUTS",true);
                    grid.profiler.addTimingInputToDictionnary("ELECTRO_MPI_COMM",true);
                }

                /* WRITE TO OUTPUT FILES */
                if( (currentStep%grid.input_parser.SAMPLING_FREQ_ELECTRO) == 0 ){

                    /// Probe time spent writing output:
                    double timeWriting = omp_get_wtime();

                    interfaceParaview.convertAndWriteData(
                        currentStep,
                        "ELECTRO"
                    );

                    timeWriting = omp_get_wtime()-timeWriting;
                    grid.profiler.incrementTimingInput("ELECTRO_WRITING_OUTPUTS",timeWriting);
                    if(grid.MPI_communicator.isRootProcess() != INT_MIN){
                        printf("%s[MPI %d - Electro - Update - step %zu]%s"
                               " Time for writing : %lf seconds.\n",
                               ANSI_COLOR_GREEN,
                               grid.MPI_communicator.getRank(),
                               currentStep,
                               ANSI_COLOR_RESET,
                               timeWriting);
                    }
                }

                grid.profiler.incrementTimingInput("ELECTRO_MPI_COMM",total_mpi_comm);

                if(grid.MPI_communicator.isRootProcess() != INT_MIN){
                        printf("%s[MPI %d - Electro - Update - step %zu]%s\n"
                               "\t> Current simulation time is   %.12lf seconds (over %.12lf).\n"
                               "\t> Current step is              %zu over %zu.\n"
                               "\t> Time elapsed inside while is %.6lf seconds (TOTAL).\n"
                               "\t> Time elaspsed per iter. is   %.6lf seconds (ITER).\n"
                               "\t> Time spent in MPI comm. is   %.6lf seconds (TOTAL).\n"
                               "\t> Time spent in MPI comm. is   %.6lf seconds (ITER).\n"
                               "\t> Time spent in writing is     %.6lf seconds (TOTAL).\n"
                               "\t> Using %d MPI process(es) and %d OMP thread(s).\n\n",
                               ANSI_COLOR_GREEN,
                               grid.MPI_communicator.getRank(),
                               currentStep,
                               ANSI_COLOR_RESET,
                               current_time,
                               grid.input_parser.get_stopTime(),
                               currentStep,
                               grid.input_parser.maxStepsForOneCycleOfElectro,
                               total_while_iter,
                               total_while_iter/currentStep,
                               grid.profiler.getTimingInput("ELECTRO_MPI_COMM"),
                               total_mpi_comm,
                               grid.profiler.getTimingInput("ELECTRO_WRITING_OUTPUTS"),
                               grid.MPI_communicator.getNumberOfMPIProcesses(),
                               omp_get_num_threads()
                               );
                    }

                total_mpi_comm = 0;

                current_time += dt;
                
            }
            #pragma omp barrier
                        

        } /* END OF WHILE LOOP */

    }/* END OF PARALLEL REGION */


    /* FREE MEMORY */
    
    delete[] local_nodes_inside_source_NUMBER;
    delete[] ID_Source;
    
    // Free H_x coefficients:
    size = grid.size_Hx[0]*grid.size_Hx[1]*grid.size_Hx[2];
    delete[] C_hxh;
    delete[] C_hxe_1;
    delete[] C_hxe_2;


    // Free H_y coefficents:
    size = grid.size_Hy[0]*grid.size_Hy[1]*grid.size_Hy[2];
    delete[] C_hyh;
    delete[] C_hye_1;
    delete[] C_hye_2;

    // Free H_z coefficients:
    size = grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2];
    delete[] C_hzh;
    delete[] C_hze_1;
    delete[] C_hze_2;

    // Free E_x coefficients:
    size = grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2];
    delete[] C_exe;
    delete[] C_exh_1;
    delete[] C_exh_2;

    // Free E_y coefficients:
    size = grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2];
    delete[] C_eye;
    delete[] C_eyh_1;
    delete[] C_eyh_2;

    // Free E_z coefficients:
    size = grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2];
    delete[] C_eze;
    delete[] C_ezh_1;
    delete[] C_ezh_2;

    /// Free electric and magnetic fields:
    for(unsigned int i = 0 ; i < NBR_FACES_CUBE ; i ++){
        if(Electric_field_to_send[i] != NULL)
            free(Electric_field_to_send[i]);
        if(Electric_field_to_recv[i] != NULL)
            free(Electric_field_to_recv[i]);
        if(Magnetic_field_to_recv[i] != NULL)
            free(Magnetic_field_to_recv[i]);
        if(Magnetic_field_to_send[i] != NULL)
            free(Magnetic_field_to_send[i]);
    }
    if(Electric_field_to_send != NULL){free(Electric_field_to_send);}
    if(Electric_field_to_recv != NULL){free(Electric_field_to_recv);}
    if(Magnetic_field_to_send != NULL){free(Magnetic_field_to_send);}
    if(Magnetic_field_to_recv != NULL){free(Magnetic_field_to_recv);}

    /// Compute total elapsed time inside UPDATE:
    gettimeofday(&end, NULL);

    double delta = end.tv_sec  - start.tv_sec + 
                        (end.tv_usec - start.tv_usec) / 1.e6;

    grid.profiler.incrementTimingInput("AlgoElectro_NEW_UPDATE_gettimeofday",delta);
    std::cout << "AlgoElectro_NEW_UPDATE => Time: " << delta << " s" << std::endl;
}

void AlgoElectro_NEW::check_OMP_DYNAMIC_envVar(void){
    /* SET OMP_DYNAMIC to false */
	if(const char *omp_dynamic_env = std::getenv("OMP_DYNAMIC")){
		// Already declared. Check it is false.
		if(std::strcmp(omp_dynamic_env,"false") == 0){
			//printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
		}else{
			std::string set_env = "OMP_DYNAMIC=false";
			putenv(&set_env[0]);
			//printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
		}
	}else{
		// OMP_DYNAMIC was not declared. Declare it.
		std::string set_env = "OMP_DYNAMIC=false";
		putenv(&set_env[0]);
		//printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
	}
}





void prepare_array_to_be_sent(
                double ** Electric_field_to_send,
                double ** Magnetic_field_to_send,
                std::vector<size_t> &electric_field_sizes,
                std::vector<size_t> &magnetic_field_sizes,
                double *E_x,
                double *E_y,
                double *E_z,
                double *H_x,
                double *H_y,
                double *H_z,
                int    *mpi_rank_neighboor,
                #ifndef NDEBUG
                    std::vector<size_t> &size_faces_electric,
                    std::vector<size_t> &size_faces_magnetic,
                #endif
                bool is_electric_to_prepare
            )
{
    /////////////////////////////////////////////////////////
    /// CONVENTION FOR COMMUNICATION                      ///
    /// 1) OMP thread(0) communicates with SOUTH. (face 0)///
    /// 2) OMP thread(1) communicates with NORTH. (face 1)///
    /// 3) OMP thread(2) communicates with WEST.  (face 2)///
    /// 3) OMP thread(3) communicates with EAST.  (face 3)///
    /// 4) OMP thread(4) communicates with DOWN.  (face 4)///
    /// 5) OMP thread(5) communicates with UP.    (face 5)///
    /////////////////////////////////////////////////////////

    size_t DECAL = DECALAGE_E_SUPP;

    //if(direction == 'W'){
    /// The W direction corresponds to (-y) axis. Denoted by face number 2 (numbering starting from 0).
    if(mpi_rank_neighboor[2] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_prepare){
            /// Put E_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0]-DECAL ; i++){
                    index   = i + electric_field_sizes[0] * ( 1 + electric_field_sizes[1] * k );
                    
                    counter = i-DECAL + (k-DECAL) * (electric_field_sizes[0]-2*DECAL);

                    ASSERT( counter, <, size_faces_electric[2]);

                    Electric_field_to_send[2][counter] = E_x[index];    
                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0]-DECAL ; i++){
                    index   = i + magnetic_field_sizes[0] * ( 1 + magnetic_field_sizes[1] * k );
                    
                    counter = i-DECAL + (k-DECAL) * (magnetic_field_sizes[0]-2*DECAL);

                    ASSERT( counter, <, size_faces_magnetic[2]);

                    Magnetic_field_to_send[2][counter] = H_x[index];

                    //printf("Electric_field_to_send[2][%zu] = %lf\n",counter,Electric_field_to_send[2][counter]);
                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_y:
            /// Offset due to already put inside the face[2]:
            counter_prev_elec = (electric_field_sizes[0]-2*DECAL) * (electric_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+3]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+3] * ( 1 + electric_field_sizes[1+3] * k );

                    counter = counter_prev_elec + (i-DECAL) + (k-DECAL) * (electric_field_sizes[0+3]-2*DECAL);
                    
                    ASSERT( counter ,<, size_faces_electric[2]);

                    Electric_field_to_send[2][counter] = E_y[index];

                    //printf("Electric_field_to_send[2][%zu] = %lf\n",counter,Electric_field_to_send[2][counter]);

                }
            }
        }else{
            /// Put H_y:
            /// Offset due to already put inside the face[2]:
            counter_prev_magn = (magnetic_field_sizes[0]-2*DECAL) * (magnetic_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+3]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+3] * ( 1 + magnetic_field_sizes[1+3] * k );

                    counter = counter_prev_magn + (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0+3]-2*DECAL);
                    
                    ASSERT( counter ,<, size_faces_magnetic[2]);

                    Magnetic_field_to_send[2][counter] = H_y[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[0+3]-2*DECAL) * (electric_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+2*3] * ( 1 + electric_field_sizes[1+2*3] * k);

                    counter = counter_prev_elec + (i-DECAL) + (k-DECAL) * (electric_field_sizes[0+2*3]-2*DECAL);

                    ASSERT( counter ,<, size_faces_electric[2]);

                    Electric_field_to_send[2][counter] = E_z[index];

                    //printf("Electric_field_to_send[2][%zu] = %lf\n",counter,Electric_field_to_send[2][counter]);

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[0+3]-2*DECAL) * (magnetic_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+2*3] * ( 1 + magnetic_field_sizes[1+2*3] * k);

                    counter = counter_prev_magn + (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0+2*3]-2*DECAL);

                    ASSERT( counter ,<, size_faces_magnetic[2]);

                    Magnetic_field_to_send[2][counter] = H_z[index];

                }
            }
        }
    }
    //if(direction == 'E'){
    /// The E direction corresponds to (+y) axis. Denoted by face number 3 (numbering starting from 0).
    if(mpi_rank_neighboor[3] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_prepare){
            /// Put E_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                     * En fait, si on avait itér sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                     */
                    index      = i + electric_field_sizes[0] * ( electric_field_sizes[1]-2 + electric_field_sizes[1] * k );

                    counter = (i-DECAL) + (k-DECAL) * (electric_field_sizes[0]-2*DECAL);

                    ASSERT( counter ,<, size_faces_electric[3]);

                    Electric_field_to_send[3][counter] = E_x[index];

                }
            }
        }else{
            /// Put = H_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                     * En fait, si on avait itér sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                     */
                    index      = i + magnetic_field_sizes[0] * ( magnetic_field_sizes[1]-2 +
                                    magnetic_field_sizes[1] * k );

                    counter = (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0]-2*DECAL);

                    ASSERT( counter ,<, size_faces_magnetic[3]);

                    Magnetic_field_to_send[3][counter] = H_x[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[0]-2*DECAL) * (electric_field_sizes[2]-2*DECAL);
            #pragma omp for
            for(size_t k = DECAL ; k < electric_field_sizes[2+3]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+3]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                     * En fait, si on avait itér sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                     */
                    index = i + electric_field_sizes[0+3] * ( electric_field_sizes[1+3]-2 + electric_field_sizes[1+3] * k );

                    counter = counter_prev_elec + (i-DECAL) + (k-DECAL) * (electric_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter ,< ,size_faces_electric[3]);

                    Electric_field_to_send[3][counter] = E_y[index];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[0]-2*DECAL) * (magnetic_field_sizes[2]-2*DECAL);
            #pragma omp for
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+3]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+3]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                     * En fait, si on avait itér sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                     */
                    index = i + magnetic_field_sizes[0+3] * ( magnetic_field_sizes[1+3]-2
                            + magnetic_field_sizes[1+3] * k );

                    counter = counter_prev_magn + (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter ,< ,size_faces_magnetic[3]);

                    Magnetic_field_to_send[3][counter] = H_y[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[0+3]-2*DECAL) * (electric_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+2*3]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                     * En fait, si on avait itér sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                     */
                    index = i + electric_field_sizes[0+2*3] * ( electric_field_sizes[1+2*3]-2 + electric_field_sizes[1+2*3] * k);

                    counter = counter_prev_elec + (i-DECAL) + (k-DECAL) * (electric_field_sizes[0+2*3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[3]);

                    Electric_field_to_send[3][counter] = E_z[index];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[0+3]-2*DECAL) * (magnetic_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+2*3]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                     * En fait, si on avait itér sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                     */
                    index = i + magnetic_field_sizes[0+2*3] * ( magnetic_field_sizes[1+2*3]-2 + 
                                magnetic_field_sizes[1+2*3] * k);

                    counter = counter_prev_magn + (i-DECAL) +
                                    (k-DECAL) * (magnetic_field_sizes[0+2*3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[3]);

                    Magnetic_field_to_send[3][counter] = H_z[index];

                }
            }
        }

    }
    //if(direction == 'S'){
    /// The S direction corresponds to (+x) axis. Denoted by face number 0 (numbering starting from 0).
    if(mpi_rank_neighboor[0] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_prepare){
            /// Put E_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1]-DECAL ; j++){
                    index      = electric_field_sizes[0]-2 + electric_field_sizes[0] * ( j + electric_field_sizes[1] * k );

                    counter = (j-DECAL) + (k-DECAL) * (electric_field_sizes[1]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[0]);

                    Electric_field_to_send[0][counter] = E_x[index];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1]-DECAL ; j++){
                    index      = magnetic_field_sizes[0]-2 + magnetic_field_sizes[0] * ( j +
                                        magnetic_field_sizes[1] * k );

                    counter = (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[0]);

                    Magnetic_field_to_send[0][counter] = H_x[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[1]-2*DECAL) * (electric_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+3]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1+3]-DECAL ; j++){
                    index = electric_field_sizes[0+3]-2 + electric_field_sizes[0+3] * ( j + electric_field_sizes[1+3] * k );

                    counter = counter_prev_elec + (j-DECAL) + (k-DECAL) * (electric_field_sizes[1+3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[0]);

                    Electric_field_to_send[0][counter] = E_y[index];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[1]-2*DECAL) * (magnetic_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+3]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1+3]-DECAL ; j++){
                    index = magnetic_field_sizes[0+3]-2 + magnetic_field_sizes[0+3] * ( j + 
                                magnetic_field_sizes[1+3] * k );

                    counter = counter_prev_magn + (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1+3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[0]);

                    Magnetic_field_to_send[0][counter] = H_y[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[1+3]-2*DECAL) * (electric_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1+2*3]-DECAL ; j++){
                    index = electric_field_sizes[0+2*3]-2 
                                + electric_field_sizes[0+2*3] * ( j + electric_field_sizes[1+2*3] * k);

                    counter = counter_prev_elec + (j-DECAL) + (k-DECAL) * (electric_field_sizes[1+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_electric[0]);

                    Electric_field_to_send[0][counter] = E_z[index];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[1+3]-2*DECAL) * (magnetic_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1+2*3]-DECAL ; j++){
                    index = magnetic_field_sizes[0+2*3]-2 
                                + magnetic_field_sizes[0+2*3] * ( j + magnetic_field_sizes[1+2*3] * k);

                    counter = counter_prev_magn + (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_magnetic[0]);

                    Magnetic_field_to_send[0][counter] = H_z[index];

                }
            }
        }

    }
    //if(direction == 'N'){
    /// The N direction corresponds to (-x) axis. Denoted by face number 1 (numbering starting from 0).
    if(mpi_rank_neighboor[1] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_prepare){
            /// Put E_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1]-DECAL ; j++){
                    index      = 1 + electric_field_sizes[0] * ( j + electric_field_sizes[1] * k );

                    counter = (j-DECAL) + (k-DECAL) * (electric_field_sizes[1]-2*DECAL);

                    ASSERT( counter ,<, size_faces_electric[1]);
                    
                    Electric_field_to_send[1][counter] = E_x[index];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1]-DECAL ; j++){
                    index      = 1 + magnetic_field_sizes[0] * ( j + magnetic_field_sizes[1] * k );

                    counter = (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1]-2*DECAL);

                    ASSERT( counter ,<, size_faces_magnetic[1]);
                    
                    Magnetic_field_to_send[1][counter] = H_x[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[1]-2*DECAL) * (electric_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+3]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1+3]-DECAL ; j++){
                    index = 1 + electric_field_sizes[0+3] * ( j + electric_field_sizes[1+3] * k );

                    counter = counter_prev_elec + (j-DECAL) + (k-DECAL) * (electric_field_sizes[1+3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[1]);

                    Electric_field_to_send[1][counter] = E_y[index];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[1]-2*DECAL) * (magnetic_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+3]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1+3]-DECAL ; j++){
                    index = 1 + magnetic_field_sizes[0+3] * ( j + magnetic_field_sizes[1+3] * k );

                    counter = counter_prev_magn + (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1+3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[1]);

                    Magnetic_field_to_send[1][counter] = H_y[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[1+3]-2*DECAL) * (electric_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1+2*3]-DECAL ; j++){
                    index = 1 + electric_field_sizes[0+2*3] * ( j + electric_field_sizes[1+2*3] * k);

                    counter = counter_prev_elec + (j-DECAL) + (k-DECAL) * (electric_field_sizes[1+6]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[1]);

                    Electric_field_to_send[1][counter] = E_z[index];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[1+3]-2*DECAL) * (magnetic_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1+2*3]-DECAL ; j++){
                    index = 1 + magnetic_field_sizes[0+2*3] * ( j + magnetic_field_sizes[1+2*3] * k);

                    counter = counter_prev_magn + (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1+6]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[1]);

                    Magnetic_field_to_send[1][counter] = H_z[index];

                }
            }
        }

    }
    //if(direction == 'U'){
    /// The N direction corresponds to (+z) axis. Denoted by face number 5 (numbering starting from 0).
    if(mpi_rank_neighboor[5] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_prepare){
            /// Put E_x:
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0]-DECAL ; i++){
                    index      = i + electric_field_sizes[0] * ( j 
                            + electric_field_sizes[1] * (electric_field_sizes[2]-2) );
                    
                    counter = (i-DECAL) + (j-DECAL) * (electric_field_sizes[0]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[5]);

                    Electric_field_to_send[5][counter] = E_x[index];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0]-DECAL ; i++){
                    index      = i + magnetic_field_sizes[0] * ( j 
                            + magnetic_field_sizes[1] * (magnetic_field_sizes[2]-2) );
                    
                    counter = (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[5]);

                    Magnetic_field_to_send[5][counter] = H_x[index];

                }
            }
        }


        if(is_electric_to_prepare){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[0]-2*DECAL) * (electric_field_sizes[1]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1+3]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+3] * ( j 
                                + electric_field_sizes[1+3] * (electric_field_sizes[2+3]-2) );

                    counter = counter_prev_elec + (i-DECAL) + (j-DECAL) * (electric_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[5]);

                    Electric_field_to_send[5][counter] = E_y[index];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[0]-2*DECAL) * (magnetic_field_sizes[1]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1+3]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+3] * ( j 
                                + magnetic_field_sizes[1+3] * (magnetic_field_sizes[2+3]-2) );

                    counter = counter_prev_magn + (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[5]);

                    Magnetic_field_to_send[5][counter] = H_y[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[0+3]-2*DECAL) * (electric_field_sizes[1+3]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1+2*3]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+2*3] * ( j 
                                + electric_field_sizes[1+2*3] * (electric_field_sizes[2+2*3]-2) );

                    counter = counter_prev_elec + (i-DECAL) + (j-DECAL) * (electric_field_sizes[0+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_electric[5]);

                    Electric_field_to_send[5][counter] = E_z[index];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[0+3]-2*DECAL) * (magnetic_field_sizes[1+3]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1+2*3]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+2*3] * ( j 
                                + magnetic_field_sizes[1+2*3] * (magnetic_field_sizes[2+2*3]-2) );

                    counter = counter_prev_magn + (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_magnetic[5]);

                    Magnetic_field_to_send[5][counter] = H_z[index];

                }
            }
        }
    }
    //if(direction == 'D'){
    /// The D direction corresponds to (-z) axis. Denoted by face number 4 (numbering starting from 0).
    if(mpi_rank_neighboor[4] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_prepare){
            /// Put E_x:
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0]-DECAL ; i++){
                    index      = i + electric_field_sizes[0] * ( j + electric_field_sizes[1] * (1) );
                    
                    counter = (i-DECAL) + (j-DECAL) * (electric_field_sizes[0]-2*DECAL);

                    ASSERT( counter ,< ,size_faces_electric[4]);

                    Electric_field_to_send[4][counter] = E_x[index];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0]-DECAL ; i++){
                    index      = i + magnetic_field_sizes[0] * ( j + magnetic_field_sizes[1] * (1) );
                    
                    counter = (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0]-2*DECAL);

                    ASSERT( counter ,< ,size_faces_magnetic[4]);

                    Magnetic_field_to_send[4][counter] = H_x[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[0]-2*DECAL) * (electric_field_sizes[1]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1+3]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+3] * ( j + electric_field_sizes[1+3] * (1) );

                    counter = counter_prev_elec + (i-DECAL) + (j-DECAL) * (electric_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter ,<, size_faces_electric[4]);

                    Electric_field_to_send[4][counter] = E_y[index];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[0]-2*DECAL) * (magnetic_field_sizes[1]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1+3]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+3] * ( j + magnetic_field_sizes[1+3] * (1) );

                    counter = counter_prev_magn + (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter ,<, size_faces_magnetic[4]);

                    Magnetic_field_to_send[4][counter] = H_y[index];

                }
            }
        }

        if(is_electric_to_prepare){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[0+3]-2*DECAL) * (electric_field_sizes[1+3]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1+2*3]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+2*3] * ( j + electric_field_sizes[1+2*3] * (1) );

                    counter = counter_prev_elec + (i-DECAL) + (j-DECAL) * (electric_field_sizes[0+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_electric[4]);

                    Electric_field_to_send[4][counter] = E_z[index];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[0+3]-2*DECAL) * (magnetic_field_sizes[1+3]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1+2*3]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+2*3] * ( j + magnetic_field_sizes[1+2*3] * (1) );

                    counter = counter_prev_magn + (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_magnetic[4]);

                    Magnetic_field_to_send[4][counter] = H_z[index];

                }
            }
        }
    }
}
            
void use_received_array(
                double **Electric_field_to_recv,
                double **Magnetic_field_to_recv,
                std::vector<size_t> &electric_field_sizes,
                std::vector<size_t> &magnetic_field_sizes,
                double *E_x,
                double *E_y,
                double *E_z,
                double *H_x,
                double *H_y,
                double *H_z,
                int    *mpi_rank_neighboor,
                #ifndef NDEBUG
                    std::vector<size_t> &size_faces_electric,
                    std::vector<size_t> &size_faces_magnetic,
                #endif
                bool is_electric_to_use
            )
{
    /**
     * This function uses the received electric fields to push it inside the electric field matrix.
     */
    /////////////////////////////////////////////////////////
    /// CONVENTION FOR COMMUNICATION                      ///
    /// 1) OMP thread(0) communicates with SOUTH. (face 0)///
    /// 2) OMP thread(1) communicates with NORTH. (face 1)///
    /// 3) OMP thread(2) communicates with WEST.  (face 2)///
    /// 3) OMP thread(3) communicates with EAST.  (face 3)///
    /// 4) OMP thread(4) communicates with DOWN.  (face 4)///
    /// 5) OMP thread(5) communicates with UP.    (face 5)///
    /////////////////////////////////////////////////////////
    
    size_t DECAL = DECALAGE_E_SUPP;

    /// Depending on the direction:

    //if(direction == 'W'){
    /// The W direction corresponds to (-y) axis. Denoted by face number 2 (numbering starting from 0).
    if(mpi_rank_neighboor[2] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        
        if(is_electric_to_use){
            /// Put E_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0]-DECAL ; i++){
                    index      = i + electric_field_sizes[0] * ( 0 + electric_field_sizes[1] * k );
                    
                    counter = (i-DECAL) + (k-DECAL) * (electric_field_sizes[0]-2*DECAL);

                    ASSERT(counter, < ,size_faces_electric[2]);
                    
                    E_x[index] = Electric_field_to_recv[2][counter];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0]-DECAL ; i++){
                    index      = i + magnetic_field_sizes[0] * ( 0 + magnetic_field_sizes[1] * k );
                    
                    counter = (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0]-2*DECAL);

                    ASSERT(counter, < ,size_faces_magnetic[2]);
                    
                    H_x[index] = Magnetic_field_to_recv[2][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[0]-2*DECAL) * (electric_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+3]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+3] * ( 0 + electric_field_sizes[1+3] * k );

                    counter = counter_prev_elec + (i-DECAL) + (k-DECAL) * (electric_field_sizes[0+3]-2*DECAL);

                    ASSERT(counter ,<, size_faces_electric[2]);

                    E_y[index] = Electric_field_to_recv[2][counter];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[0]-2*DECAL) * (magnetic_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+3]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+3] * ( 0 + magnetic_field_sizes[1+3] * k );

                    counter = counter_prev_magn + (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0+3]-2*DECAL);

                    ASSERT(counter ,<, size_faces_magnetic[2]);

                    H_y[index] = Magnetic_field_to_recv[2][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[0+3]-2*DECAL) * (electric_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+2*3] * ( 0 + electric_field_sizes[1+2*3] * k);

                    counter = counter_prev_elec + (i-DECAL) + (k-DECAL) * (electric_field_sizes[0+6]-2*DECAL);

                    ASSERT(counter, < ,size_faces_electric[2]);

                    E_z[index] = Electric_field_to_recv[2][counter];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[0+3]-2*DECAL) * (magnetic_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+2*3] * ( 0 + magnetic_field_sizes[1+2*3] * k);

                    counter = counter_prev_magn + (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0+6]-2*DECAL);

                    ASSERT(counter, < ,size_faces_magnetic[2]);

                    H_z[index] = Magnetic_field_to_recv[2][counter];

                }
            }
        }
    }

    

    //if(direction == 'E'){
    /// The E direction corresponds to (+y) axis. Denoted by face number 3 (numbering starting from 0).
    if(mpi_rank_neighboor[3] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_use){
            /// Put E_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index      = i + 
                        electric_field_sizes[0] * ( electric_field_sizes[1]-1 + 
                                electric_field_sizes[1] * k );

                    counter = (i-DECAL) + (k-DECAL) * (electric_field_sizes[0]-2*DECAL);

                    ASSERT(counter ,<, size_faces_electric[3]);
                    
                    E_x[index] = Electric_field_to_recv[3][counter];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index      = i + 
                        magnetic_field_sizes[0] * ( magnetic_field_sizes[1]-1 + 
                                magnetic_field_sizes[1] * k );

                    counter = (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0]-2*DECAL);

                    ASSERT(counter ,<, size_faces_magnetic[3]);
                    
                    H_x[index] = Magnetic_field_to_recv[3][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[0]-2*DECAL) * (electric_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+3]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+3]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = i + electric_field_sizes[0+3] 
                            * ( electric_field_sizes[1+3]-1 + electric_field_sizes[1+3] * k );

                    counter = counter_prev_elec + (i-DECAL) + (k-DECAL) * (electric_field_sizes[0+3]-2*DECAL);

                    ASSERT(counter, < ,size_faces_electric[3]);

                    E_y[index] = Electric_field_to_recv[3][counter];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[0]-2*DECAL) * (magnetic_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+3]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+3]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = i + magnetic_field_sizes[0+3] 
                            * ( magnetic_field_sizes[1+3]-1 + magnetic_field_sizes[1+3] * k );

                    counter = counter_prev_magn + (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0+3]-2*DECAL);

                    ASSERT(counter, < ,size_faces_magnetic[3]);

                    H_y[index] = Magnetic_field_to_recv[3][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[0+3]-2*DECAL) * (electric_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+2*3]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = i + electric_field_sizes[0+2*3] 
                            * ( electric_field_sizes[1+2*3]-1 + electric_field_sizes[1+2*3] * k);

                    counter = counter_prev_elec + (i-DECAL) + (k-DECAL) * (electric_field_sizes[0+6]-2*DECAL);

                    ASSERT(counter, < ,size_faces_electric[3]);

                    E_z[index] = Electric_field_to_recv[3][counter];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[0+3]-2*DECAL) * (magnetic_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+2*3]-DECAL ; i++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = i + magnetic_field_sizes[0+2*3] 
                            * ( magnetic_field_sizes[1+2*3]-1 + magnetic_field_sizes[1+2*3] * k);

                    counter = counter_prev_magn + (i-DECAL) + (k-DECAL) * (magnetic_field_sizes[0+6]-2*DECAL);

                    ASSERT(counter, < ,size_faces_magnetic[3]);

                    H_z[index] = Magnetic_field_to_recv[3][counter];

                }
            }
        }
    }
    /**
     * I receive information with tag 'S'.
     */
    //if(direction == 'S'){
    /// The S direction corresponds to (+x) axis. Denoted by face number 0 (numbering starting from 0).
    if(mpi_rank_neighboor[0] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_use){
            /// Put E_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index      = electric_field_sizes[0]-1 
                            + electric_field_sizes[0] * ( j + electric_field_sizes[1] * k );

                    counter = (j-DECAL) + (k-DECAL) * (electric_field_sizes[1]-2*DECAL);

                    ASSERT(counter, <, size_faces_electric[0]);
                    
                    E_x[index] = Electric_field_to_recv[0][counter];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index      = magnetic_field_sizes[0]-1 
                            + magnetic_field_sizes[0] * ( j + magnetic_field_sizes[1] * k );

                    counter = (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1]-2*DECAL);

                    ASSERT(counter, <, size_faces_magnetic[0]);
                    
                    H_x[index] = Magnetic_field_to_recv[0][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[1]-2*DECAL) * (electric_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+3]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1+3]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = electric_field_sizes[0+3]-1 
                            + electric_field_sizes[0+3] * ( j + electric_field_sizes[1+3] * k );

                    counter = counter_prev_elec + (j-DECAL) + (k-DECAL) * (electric_field_sizes[1+3]-2*DECAL);

                    ASSERT(counter ,< ,size_faces_electric[0]);

                    E_y[index] = Electric_field_to_recv[0][counter];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[1]-2*DECAL) * (magnetic_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+3]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1+3]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = magnetic_field_sizes[0+3]-1 
                            + magnetic_field_sizes[0+3] * ( j + magnetic_field_sizes[1+3] * k );

                    counter = counter_prev_magn + (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1+3]-2*DECAL);

                    ASSERT(counter ,< ,size_faces_magnetic[0]);

                    H_y[index] = Magnetic_field_to_recv[0][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[1+3]-2*DECAL) * (electric_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1+2*3]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = electric_field_sizes[0+2*3]-1 
                            + electric_field_sizes[0+2*3] * ( j + electric_field_sizes[1+2*3] * k);

                    counter = counter_prev_elec + (j-DECAL) + (k-DECAL) * (electric_field_sizes[1+6]-2*DECAL);

                    ASSERT(counter ,<, size_faces_electric[0]);

                    E_z[index] = Electric_field_to_recv[0][counter];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[1+3]-2*DECAL) * (magnetic_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1+2*3]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = magnetic_field_sizes[0+2*3]-1 
                            + magnetic_field_sizes[0+2*3] * ( j + magnetic_field_sizes[1+2*3] * k);

                    counter = counter_prev_magn + (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1+6]-2*DECAL);

                    ASSERT(counter ,<, size_faces_magnetic[0]);

                    H_z[index] = Magnetic_field_to_recv[0][counter];

                }
            }
        }
    }
    //if(direction == 'N'){
    /// The N direction corresponds to (-x) axis. Denoted by face number 1 (numbering starting from 0).
    if(mpi_rank_neighboor[1] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_use){
            /// Put E_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index      = 0 + electric_field_sizes[0] * ( j + electric_field_sizes[1] * k );

                    counter = (j-DECAL) + (k-DECAL) * (electric_field_sizes[1]-2*DECAL);

                    ASSERT(counter, <, size_faces_electric[1]);
                    
                    E_x[index] = Electric_field_to_recv[1][counter];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index      = 0 + magnetic_field_sizes[0] * ( j + magnetic_field_sizes[1] * k );

                    counter = (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1]-2*DECAL);

                    ASSERT(counter, <, size_faces_magnetic[1]);
                    
                    H_x[index] = Magnetic_field_to_recv[1][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[1]-2*DECAL) * (electric_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+3]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1+3]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = 0 + electric_field_sizes[0+3] * ( j + electric_field_sizes[1+3] * k );

                    counter = counter_prev_elec + (j-DECAL) + (k-DECAL) * (electric_field_sizes[1+3]-2*DECAL);

                    ASSERT(counter, < ,size_faces_electric[1]);

                    E_y[index] = Electric_field_to_recv[1][counter];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[1]-2*DECAL) * (magnetic_field_sizes[2]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+3]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1+3]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = 0 + magnetic_field_sizes[0+3] * ( j + magnetic_field_sizes[1+3] * k );

                    counter = counter_prev_magn + (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1+3]-2*DECAL);

                    ASSERT(counter, < ,size_faces_magnetic[1]);

                    H_y[index] = Magnetic_field_to_recv[1][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[1+3]-2*DECAL) * (electric_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < electric_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t j = DECAL ; j < electric_field_sizes[1+2*3]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = 0 
                            + electric_field_sizes[0+2*3] * ( j + electric_field_sizes[1+2*3] * k);

                    counter = counter_prev_elec + (j-DECAL) + (k-DECAL) * (electric_field_sizes[1+6]-2*DECAL);

                    ASSERT(counter ,< ,size_faces_electric[1]);

                    E_z[index] = Electric_field_to_recv[1][counter];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[1+3]-2*DECAL) * (magnetic_field_sizes[2+3]-2*DECAL);
            #pragma omp for 
            for(size_t k = DECAL ; k < magnetic_field_sizes[2+2*3]-DECAL ; k++){
                for(size_t j = DECAL ; j < magnetic_field_sizes[1+2*3]-DECAL ; j++){
                    /**
                     * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                     * En fait, si on avait itéré sur les j, on aurait fait:
                     *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                     *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                     */
                    index = 0 
                            + magnetic_field_sizes[0+2*3] * ( j + magnetic_field_sizes[1+2*3] * k);

                    counter = counter_prev_magn + (j-DECAL) + (k-DECAL) * (magnetic_field_sizes[1+6]-2*DECAL);

                    ASSERT(counter ,< ,size_faces_magnetic[1]);

                    H_z[index] = Magnetic_field_to_recv[1][counter];

                }
            }
        }
    }
    //if(direction == 'U'){
    /// The N direction corresponds to (+z) axis. Denoted by face number 5 (numbering starting from 0).
    if(mpi_rank_neighboor[5] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_use){
            /// Put E_x:
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0]-DECAL ; i++){
                    index      = i + electric_field_sizes[0] 
                            * ( j + electric_field_sizes[1] * (electric_field_sizes[2]-1) );
                    
                    counter = (i-DECAL) + (j-DECAL) * (electric_field_sizes[0]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[5]);

                    E_x[index] = Electric_field_to_recv[5][counter];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0]-DECAL ; i++){
                    index      = i + magnetic_field_sizes[0] 
                            * ( j + magnetic_field_sizes[1] * (magnetic_field_sizes[2]-1) );
                    
                    counter = (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[5]);

                    H_x[index] = Magnetic_field_to_recv[5][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[0]-2*DECAL) * (electric_field_sizes[1]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1+3]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+3] 
                        * ( j + electric_field_sizes[1+3] * (electric_field_sizes[2+3]-1) );

                    counter = counter_prev_elec + (i-DECAL) + (j-DECAL) * (electric_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_electric[5]);

                    E_y[index] = Electric_field_to_recv[5][counter];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[0]-2*DECAL) * (magnetic_field_sizes[1]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1+3]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+3] 
                        * ( j + magnetic_field_sizes[1+3] * (magnetic_field_sizes[2+3]-1) );

                    counter = counter_prev_magn + (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter, < ,size_faces_magnetic[5]);

                    H_y[index] = Magnetic_field_to_recv[5][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[0+3]-2*DECAL) * (electric_field_sizes[1+3]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1+2*3]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + 
                        electric_field_sizes[0+2*3] * ( j + 
                            electric_field_sizes[1+2*3] * (electric_field_sizes[2+2*3]-1) );

                    counter = counter_prev_elec + (i-DECAL) + (j-DECAL) * (electric_field_sizes[0+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_electric[5]);

                    E_z[index] = Electric_field_to_recv[5][counter];

                }
            }
        }else{
            /// Put E_z:
            counter_prev_magn += (magnetic_field_sizes[0+3]-2*DECAL) * (magnetic_field_sizes[1+3]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1+2*3]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + 
                        magnetic_field_sizes[0+2*3] * ( j + 
                            magnetic_field_sizes[1+2*3] * (magnetic_field_sizes[2+2*3]-1) );

                    counter = counter_prev_magn + (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_magnetic[5]);

                    H_z[index] = Magnetic_field_to_recv[5][counter];

                }
            }
        }
    }
    //if(direction == 'D'){
    /// The D direction corresponds to (-z) axis. Denoted by face number 4 (numbering starting from 0).
    if(mpi_rank_neighboor[4] != -1){
        size_t index   = 0;
        size_t counter = 0;
        size_t counter_prev_elec = 0;
        size_t counter_prev_magn = 0;

        if(is_electric_to_use){
            /// Put E_x:
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0]-DECAL ; i++){
                    index      = i + electric_field_sizes[0] * ( j + electric_field_sizes[1] * (0) );
                    
                    counter = (i-DECAL) + (j-DECAL) * (electric_field_sizes[0]-2*DECAL);

                    ASSERT( counter ,< ,size_faces_electric[4]);

                    E_x[index] = Electric_field_to_recv[4][counter];

                }
            }
        }else{
            /// Put H_x:
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0]-DECAL ; i++){
                    index      = i + magnetic_field_sizes[0] * ( j + magnetic_field_sizes[1] * (0) );
                    
                    counter = (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0]-2*DECAL);

                    ASSERT( counter ,< ,size_faces_electric[4]);

                    H_x[index] = Magnetic_field_to_recv[4][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_y:
            counter_prev_elec = (electric_field_sizes[0]-2*DECAL) * (electric_field_sizes[1]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1+3]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+3] * ( j + electric_field_sizes[1+3] * (0) );

                    counter = counter_prev_elec + (i-DECAL) + (j-DECAL) * (electric_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter ,<, size_faces_electric[4]);

                    E_y[index] = Electric_field_to_recv[4][counter];

                }
            }
        }else{
            /// Put H_y:
            counter_prev_magn = (magnetic_field_sizes[0]-2*DECAL) * (magnetic_field_sizes[1]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < magnetic_field_sizes[1+3]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+3] * ( j + magnetic_field_sizes[1+3] * (0) );

                    counter = counter_prev_magn + (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0+3]-2*DECAL);

                    ASSERT( counter ,<, size_faces_magnetic[4]);

                    H_y[index] = Magnetic_field_to_recv[4][counter];

                }
            }
        }

        if(is_electric_to_use){
            /// Put E_z:
            counter_prev_elec += (electric_field_sizes[0+3]-2*DECAL) * (electric_field_sizes[1+3]-2*DECAL);
            #pragma omp for 
            for(size_t j = DECAL ; j < electric_field_sizes[1+2*3]-DECAL ; j++){
                for(size_t i = DECAL ; i < electric_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + electric_field_sizes[0+2*3] * ( j + electric_field_sizes[1+2*3] * (0) );

                    counter = counter_prev_elec + (i-DECAL) + (j-DECAL) * (electric_field_sizes[0+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_electric[4]);

                    E_z[index] = Electric_field_to_recv[4][counter];

                }
            }
        }else{
            /// Put H_z:
            counter_prev_magn += (magnetic_field_sizes[0+3]-2*DECAL) * (magnetic_field_sizes[1+3]-2*DECAL);
            #pragma omp for
            for(size_t j = DECAL ; j < magnetic_field_sizes[1+2*3]-DECAL ; j++){
                for(size_t i = DECAL ; i < magnetic_field_sizes[0+2*3]-DECAL ; i++){
                    index = i + magnetic_field_sizes[0+2*3] * ( j + magnetic_field_sizes[1+2*3] * (0) );

                    counter = counter_prev_magn + (i-DECAL) + (j-DECAL) * (magnetic_field_sizes[0+6]-2*DECAL);

                    ASSERT( counter, <, size_faces_magnetic[4]);

                    H_z[index] = Magnetic_field_to_recv[4][counter];

                }
            }
        }
    }
}

/**
 * Communicates both electric and magnetic fields:
 */
void communicate_single_omp_thread(
                double **Electric_field_to_send,
                double **Electric_field_to_recv,
                double **Magnetic_field_to_send,
                double **Magnetic_field_to_recv,
                int *mpi_to_who,
                int  mpi_me,
                std::vector<size_t> size_faces_electric,
                std::vector<size_t> size_faces_magnetic,
                bool is_electric_to_communicate
            )
{
    /// Only the master OPENMP thread can access this fuction !
    if(omp_get_thread_num() != 0){
        fprintf(stderr,"In %s :: only the master OMP thread can accessthis function ! Aborting.\n",
            __FUNCTION__);
        fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
        #ifdef MPI_COMM_WORLD
        MPI_Abort(MPI_COMM_WORLD,-1);
        #else
        abort();
        #endif
    }

    #ifndef NDEBUG
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        printf("[MPI %d] - NEIGHBOORS [%d,%d,%d,%d,%d,%d]\n",
            mpi_me,
            mpi_to_who[0],mpi_to_who[1],mpi_to_who[2],mpi_to_who[3],
            mpi_to_who[4],mpi_to_who[5]);
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
    #endif

    /// LOOP OVER THE 6 FACES
    for(unsigned int FACE = 0 ; FACE < NBR_FACES_CUBE ; FACE ++){

        // If it is -1, then no need to communicate ! Just continue.
        if(mpi_to_who[FACE] == -1){continue;}

        /// Define a tag for communication:
        int neighboorComm = -1;

        if( FACE == 0 ){ neighboorComm = 1; }
        if( FACE == 1 ){ neighboorComm = 0; }
        if( FACE == 2 ){ neighboorComm = 3; }
        if( FACE == 3 ){ neighboorComm = 2; }
        if( FACE == 4 ){ neighboorComm = 5; }
        if( FACE == 5 ){ neighboorComm = 4; }

        /// The order of send/recv operations depends on the MPI rank of current and neighboor MPI processes.

        /* NEIGHBOOR IS EVEN, I AM ODD */
        if(mpi_me%2 != 0 && mpi_to_who[FACE]%2 == 0){
            #ifndef NDEBUG
                printf("[MPI %d - EVEN - FACE %d -OMP %d] send to   [MPI %d] | sendTag %d , recvTag %d\n",
                        mpi_me,
                        FACE,
                        omp_get_thread_num(),
                        mpi_to_who[FACE],
                        FACE,
                        neighboorComm);
            #endif
            if(is_electric_to_communicate){
                MPI_Send(
                        Electric_field_to_send[FACE],
                        size_faces_electric[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        FACE,
                        MPI_COMM_WORLD
                );
            }else{
                MPI_Send(
                        Magnetic_field_to_send[FACE],
                        size_faces_magnetic[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        FACE,
                        MPI_COMM_WORLD
                );
            }
            #ifndef NDEBUG
                printf("[MPI %d - EVEN - FACE %d -OMP %d] recv from [MPI %d] | sendTag %d , recvTag %d\n",
                        mpi_me,
                        FACE,
                        omp_get_thread_num(),
                        mpi_to_who[FACE],
                        FACE,
                        neighboorComm);
            #endif
            if(is_electric_to_communicate){
                MPI_Recv(
                        Electric_field_to_recv[FACE],
                        size_faces_electric[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        neighboorComm,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
            }else{
                MPI_Recv(
                        Magnetic_field_to_recv[FACE],
                        size_faces_magnetic[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        neighboorComm,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
            }
            #ifndef NDEBUG
                printf("[MPI %d - ODD ] to [MPI %d] :: Done\n",
                        mpi_me,
                        mpi_to_who[FACE]);
            #endif
                
        /* NEIGHBOOR IS ODD, I AM EVEN */
        }else if(mpi_me%2 == 0 && mpi_to_who[FACE]%2 != 0){
            #ifndef NDEBUG
                printf("[MPI %d - ODD - FACE %d - OMP %d] recv from [MPI %d] | sendTag %d , recvTag %d\n",
                        mpi_me,
                        FACE,
                        omp_get_thread_num(),
                        mpi_to_who[FACE],
                        FACE,
                        neighboorComm);
            #endif
            if(is_electric_to_communicate){
                MPI_Recv(
                        Electric_field_to_recv[FACE],
                        size_faces_electric[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        neighboorComm,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
            }else{
                MPI_Recv(
                        Magnetic_field_to_recv[FACE],
                        size_faces_magnetic[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        neighboorComm,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
            }
            #ifndef NDEBUG
                printf("[MPI %d - ODD - FACE %d - OMP %d] send to   [MPI %d] | sendTag %d , recvTag %d\n",
                        mpi_me,
                        FACE,
                        omp_get_thread_num(),
                        mpi_to_who[FACE],
                        FACE,
                        neighboorComm);
            #endif
            if(is_electric_to_communicate){
                MPI_Send(
                        Electric_field_to_send[FACE],
                        size_faces_electric[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        FACE,
                        MPI_COMM_WORLD
                );
            }else{
                MPI_Send(
                        Magnetic_field_to_send[FACE],
                        size_faces_magnetic[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        FACE,
                        MPI_COMM_WORLD
                );
            }
            #ifndef NDEBUG
                printf("[MPI %d - EVEN] to [MPI %d] :: DONE\n",
                        mpi_me,
                        mpi_to_who[FACE]);
            #endif
                
        /* NEIGHBOOR LARGER THAN ME */
        }else if( mpi_to_who[FACE] > mpi_me ){
            #ifndef NDEBUG
                printf("[MPI %d - LARGER - FACE %d -OMP %d] to [MPI %d] | sendTag %d , recvTag %d\n",
                        mpi_me,
                        FACE,
                        omp_get_thread_num(),
                        mpi_to_who[FACE],
                        FACE,
                        neighboorComm);
            #endif
            if(is_electric_to_communicate){
                MPI_Recv(
                        Electric_field_to_recv[FACE],
                        size_faces_electric[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        neighboorComm,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
            }else{
                MPI_Recv(
                        Magnetic_field_to_recv[FACE],
                        size_faces_magnetic[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        neighboorComm,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
            }
            if(is_electric_to_communicate){   
                MPI_Send(
                        Electric_field_to_send[FACE],
                        size_faces_electric[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        FACE,
                        MPI_COMM_WORLD
                );
            }else{
                MPI_Send(
                        Magnetic_field_to_send[FACE],
                        size_faces_magnetic[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        FACE,
                        MPI_COMM_WORLD
                );
            } 
            #ifndef NDEBUG
                printf("[MPI %d - LARGER] to [MPI %d] :: DONE\n",
                        mpi_me,
                        mpi_to_who[FACE]);
            #endif
                
        }else if( mpi_to_who[FACE] < mpi_me ){
            #ifndef NDEBUG
                printf("[MPI %d - SMALLER - FACE %d - OMP %d] to [MPI %d] | sendTag %d , recvTag %d\n",
                        mpi_me,
                        FACE,
                        omp_get_thread_num(),
                        mpi_to_who[FACE],
                        FACE,
                        neighboorComm);
            #endif
            if(is_electric_to_communicate){
                MPI_Send(
                        Electric_field_to_send[FACE],
                        size_faces_electric[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        FACE,
                        MPI_COMM_WORLD
                );
            }else{
                MPI_Send(
                        Magnetic_field_to_send[FACE],
                        size_faces_magnetic[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        FACE,
                        MPI_COMM_WORLD
                );
            }
            if(is_electric_to_communicate){  
                MPI_Recv(
                        Electric_field_to_recv[FACE],
                        size_faces_electric[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        neighboorComm,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
            }else{
                MPI_Recv(
                        Magnetic_field_to_recv[FACE],
                        size_faces_magnetic[FACE],
                        MPI_DOUBLE,
                        mpi_to_who[FACE],
                        neighboorComm,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
            }  
            #ifndef NDEBUG
                printf("[MPI %d - SMALLER] to [MPI %d] :: DONE\n",
                        mpi_me,
                        mpi_to_who[FACE]);
            #endif
                
        }else{
                fprintf(stderr,"In function %s :: no way to communicate between MPI %d and"
                                " MPI %d, on face %d. Aborting.\n",
                                __FUNCTION__,mpi_me,mpi_to_who[FACE],FACE);
                fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
                #ifdef MPI_COMM_WORLD
                    MPI_Abort(MPI_COMM_WORLD,-1);
                #else
                    abort();
                #endif
        }
    }
    
}