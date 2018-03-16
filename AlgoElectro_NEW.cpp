#include "AlgoElectro_NEW.hpp"

#include <ctime>
#include <new>
#include "omp.h"
#include "mpi.h"
#include <algorithm>


void MPI_communication_with_neighboors(
    int  omp_thread_id,
    int  mpi_me,
    int  mpi_from_to_who,
    char direction,
    double *ElectricFieldToSend,
    size_t size_to_send,
    double *ElectricFieldToRecv,
    size_t size_to_recv,
    double *E_x,
    double *E_y,
    double *E_z,
    std::vector<size_t> &omp_sizes
);


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

void AlgoElectro_NEW::update(
    GridCreator_NEW &grid,
    InterfaceToParaviewer &interfaceParaview){

    // Start clock for monitoring CPU time:
    double start_algo_update;
    double end___algo_update;
    grid.profiler.addTimingInputToDictionnary("AlgoElectro_NEW_UPDATE_omp_get_wtime");
    start_algo_update = omp_get_wtime();

    // Retrieve the time step:
    double dt = this->Compute_dt(grid);
    double current_time = 0.0;
    size_t currentStep = 0;
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

    // In the object grid, set the properties mu, eps, magnetic cond. and electric cond. f each node:
    grid.Initialize_Electromagnetic_Properties("AIR_AT_INIT_TEMP");

    /* Set the coefficients for the magnetic update */

    size_t size;
    size_t memory;
    
    //size_t M = grid.sizes_EH[0];
    //size_t N = grid.sizes_EH[1];
    //size_t P = grid.sizes_EH[2];

    // Magnetic field Hx of size M × (N − 1) × (P − 1):
    size = grid.size_Hx[0] * grid.size_Hx[1] * grid.size_Hx[2];
    double *C_hxh   = new double[size];
    double *C_hxe_1 = new double[size];
    double *C_hxe_2 = new double[size];

    memory = (8+8+8) * size;
    grid.profiler.addMemoryUsage("BYTES",memory);

    // Magnetic field Hy of size (M − 1) × N × (P − 1):
    size = grid.size_Hy[0] * grid.size_Hy[1] * grid.size_Hy[2];
    double *C_hyh   = new double[size];
    double *C_hye_1 = new double[size];
    double *C_hye_2 = new double[size];

    memory = (8+8+8) * size;
    grid.profiler.addMemoryUsage("BYTES",memory);

    // Magnetic field Hz of size (M − 1) × (N − 1) × P:
    size = grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2];
    double *C_hzh   = new double[size];
    double *C_hze_1 = new double[size];
    double *C_hze_2 = new double[size];

    memory = (8+8+8) * size;
    grid.profiler.addMemoryUsage("BYTES",memory);

    // Electric field Ex of size (M − 1) × N × P:
    size = grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2];
    double *C_exe   = new double[size];
    double *C_exh_1 = new double[size];
    double *C_exh_2 = new double[size];

    memory = (8+8+8) * size;
    grid.profiler.addMemoryUsage("BYTES",memory);

    // Electric field Ey of size M × (N − 1) × P:
    size = grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]; 
    double *C_eye   = new double[size];
    double *C_eyh_1 = new double[size];
    double *C_eyh_2 = new double[size];

    memory = (8+8+8) * size;
    grid.profiler.addMemoryUsage("BYTES",memory);

    // Electric field Ez of size M × N × (P − 1):
    size = grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2];
    double *C_eze   = new double[size];
    double *C_ezh_1 = new double[size];
    double *C_ezh_2 = new double[size];

    memory = (8+8+8) * size;
    grid.profiler.addMemoryUsage("BYTES",memory);

    /* COMPUTING COEFFICIENTS */
    if(omp_get_max_threads() < 6 ){
        omp_set_num_threads(6);
    }
    #pragma omp parallel
    {
        if(omp_get_num_threads() < 6){
            fprintf(stderr,"AlgoElectro_NEW::UPDATE::ERROR\n");
            fprintf(stderr,"Not enough OMP threads. Needs 6 but has %d.\n",omp_get_num_threads());
            fprintf(stderr,"Aborting in file %s:%d.\n",__FILE__,__LINE__);
            abort();
        }
    }

    #pragma omp parallel num_threads(1)
    {
        size_t index;
        /* Coefficients for Ex */
        // Ex of size (M − 1) × N × P
        #pragma omp for collapse(3) nowait
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
        #pragma omp for collapse(3) nowait
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
        #pragma omp for collapse(3) nowait
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
        #pragma omp for collapse(3) nowait
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
        #pragma omp for collapse(3) nowait
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
        #pragma omp for collapse(3) nowait
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
    // END OF #pragma omp parallel num_threads(6)

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

    // Size is 6 because 3 E et 3 H components:
    std::vector<size_t>        nbr_nodes_inside_source = {0,0,0,0,0,0};

    std::vector<std::string> TYPE = {"Ex","Ey","Ez","Hx","Hy","Hz"};

    for(unsigned int i = 0 ; i < TYPE.size() ; i ++){
        grid.Compute_nodes_inside_sources(
            local_nodes_inside_source_NUMBER[i],
            ID_Source[i],
            local_nodes_inside_source_FREQ,
            &nbr_nodes_inside_source[i],
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
        grid.profiler.addMemoryUsage("BYTES",sizeof(size_t)*local_nodes_inside_source_NUMBER[i].size());
        grid.profiler.addMemoryUsage("BYTES",sizeof(unsigned char)*ID_Source[i].size());

        printf("For %s, found %zu nodess !\n",TYPE[i].c_str(),ID_Source[i].size());

    }

    fprintf(stderr,"\n\t>>> Ready to compute sa mère !!!\n\n");

    ///////////////////////////////////////////////
    // UPDATE WHILE LOOP - PARALLELIZED WITH     //
    // OPENMP THREADS    - MINIMUM 6 OPENMP      //
    // THREADS ARE REQUIRED FOR MPI COMMUICATION //
    // TO WORK.                                  //
    ///////////////////////////////////////////////

    // Temporary pointers, to avoid doing grid.sthg !
    double *H_x_tmp = grid.H_x;
    double *H_y_tmp = grid.H_y;
    double *H_z_tmp = grid.H_z;
    double *E_x_tmp = grid.E_x;
    double *E_y_tmp = grid.E_y;
    double *E_z_tmp = grid.E_z;

    /// Set OMP_DYNAMIC=false:
    this->check_OMP_DYNAMIC_envVar();


    double parallelRegionStartingTime = omp_get_wtime();

    #pragma omp parallel num_threads(omp_get_max_threads()) default(none)\
        shared(grid,dt,current_time,currentStep)\
        shared(local_nodes_inside_source_NUMBER)\
        shared(interfaceParaview)\
        shared(parallelRegionStartingTime)\
        shared(H_x_tmp,H_y_tmp,H_z_tmp)\
        shared(E_x_tmp,E_y_tmp,E_z_tmp)\
        shared(C_hxh,C_hxe_1,C_hxe_2)\
        shared(C_hyh,C_hye_1,C_hye_2)\
        shared(C_hzh,C_hze_1,C_hze_2)\
        shared(C_exe,C_exh_1,C_exh_2)\
        shared(C_eye,C_eyh_1,C_eyh_2)\
        shared(C_eze,C_ezh_1,C_ezh_2)\
        shared(ompi_mpi_comm_world,ompi_mpi_int)
    {
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

        /**
         * Important for the electric field update !
         */
        size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX = 0;
        size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY = 0;
        size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ = 0;

        if(grid.MPI_communicator.MPI_POSITION[0] == grid.MPI_communicator.MPI_MAX_POSI[0])
            IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX = 1;

        if(grid.MPI_communicator.MPI_POSITION[1] == grid.MPI_communicator.MPI_MAX_POSI[1])
            IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY = 1;

        if(grid.MPI_communicator.MPI_POSITION[2] == grid.MPI_communicator.MPI_MAX_POSI[2])
            IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ = 1;

            

        size_t I,J,K;

        ////////////////////////
        // MPI INITIALIZATION //
        ////////////////////////
        bool OMP_thread_has_neighboor = false;
        char direction = '\n';

        /// Put all the electric field sizes inside one vector for convenience.
        std::vector<size_t> omp_sizes = {
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

        /// Vectors to store data to be sent or received when communicating with
        /// other MPI processes.
        size_t size_to_be_sent = 0;
        size_t size_to_recv    = 0;
        double *ElectricFieldToSend = NULL;
        double *ElectricFieldToRecv = NULL;

        /// Determine therole of each OMP thread in the MPI communication process:
        this->determine_OMP_thread_role_in_MPI_communication(
            omp_get_thread_num(),
            &OMP_thread_has_neighboor,
            &direction,
            grid,
            omp_sizes,
            &ElectricFieldToSend,
            &ElectricFieldToRecv,
            &size_to_be_sent,
            &size_to_recv
        );

        while(current_time < grid.input_parser.get_stopTime()
                && currentStep < grid.input_parser.maxStepsForOneCycleOfElectro){


            // Updating the magnetic field Hx.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Hx[0];
            size_y = grid.size_Hx[1];

            size_x_1 = grid.size_Ey[0];
            size_y_1 = grid.size_Ey[1];

            size_x_2 = grid.size_Ez[0];
            size_y_2 = grid.size_Ez[1];


            #pragma omp for schedule(static) collapse(3) nowait
            for (K = 0; K < grid.size_Hx[2] ; K++){
                for(J = 0; J < grid.size_Hx[1] ; J ++){
                    for(I = 0; I < grid.size_Hx[0] ; I++){

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

            #pragma omp for schedule(static) collapse(3) nowait
            for(K = 0 ; K < grid.size_Hy[2] ; K ++){
                for(J = 0 ; J < grid.size_Hy[1] ; J ++){
                    for(I = 0; I < grid.size_Hy[0] ; I ++){

                        index        = I   + size_x   * ( J  + size_y   * K);
                        // Ez(mm + 1, nn, pp):
                        index_1Plus  = I+1 + size_x_1 * ( J  + size_y_1 * K);
                        // Ez(mm, nn, pp)
                        index_1Moins = I   + size_x_1 * ( J  + size_y_1 * K);
                        // Ex(mm, nn, pp + 1):
                        index_2Plus  = I   + size_x_2 * ( J  + size_y_2 * (K+1));
                        // Ex(mm, nn, pp):
                        index_2Moins = I   + size_x_2 * ( J  + size_y_2 * K);

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

            #pragma omp for schedule(static) collapse(3) nowait
            for(K = 0 ; K < grid.size_Hz[2]  ; K ++){
                for(J = 0 ; J < grid.size_Hz[1] ; J ++){
                    for(I = 0 ; I < grid.size_Hz[0] ; I ++){

                        index        = I   + size_x   * ( J     + size_y   * K);
                        // Ex(mm, nn + 1, pp)
                        index_1Plus  = I   + size_x_1 * ( (J+1) + size_y_1 * K);
                        // Ex(mm, nn, pp)
                        index_1Moins = I   + size_x_1 * ( J     + size_y_1 * K);
                        // Ey(mm + 1, nn, pp)
                        index_2Plus  = I+1 + size_x_2 * ( J     + size_y_2 * K);
                        // Ey(mm, nn, pp)
                        index_2Moins = I   + size_x_2 * ( J     + size_y_2 * K);

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

            // Updating the electric field Ex.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Ex[0];
            size_y = grid.size_Ex[1];

            size_x_1 = grid.size_Hz[0];
            size_y_1 = grid.size_Hz[1];

            size_x_2 = grid.size_Hy[0];
            size_y_2 = grid.size_Hy[1];


            #pragma omp for schedule(static) collapse(3) nowait
            for(K = 1  ; K < grid.size_Ex[2]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ ; K ++){
                for(J = 1  ; J < grid.size_Ex[1]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY ; J ++){
                    for(I = 1 ; I < grid.size_Ex[0]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX; I ++){

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

            #pragma omp for schedule(static) collapse(3) nowait
            for(K = 1  ; K < grid.size_Ey[2]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ ; K ++){
                for(J = 1 ; J < grid.size_Ey[1]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY ; J ++){
                    for(I = 1  ; I < grid.size_Ey[0]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX ; I ++){

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
                        index_2Plus  = I/*-DECALAGE_H_in_E   */
                                    + size_x_2 * ( J 
                                         + size_y_2 * (K ));
                        // Hz(mm - 1, nn, pp)
                        index_2Moins = I-1 
                                    + size_x_2 * ( J 
                                         + size_y_2 * (K ));

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

            #pragma omp for schedule(static) collapse(3) nowait
            for(K = 1 ; K < grid.size_Ez[2]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ ; K ++){
                for(J = 1  ; J < grid.size_Ez[1]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY ; J ++){
                    for(I = 1  ; I < grid.size_Ez[0]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX ; I ++){

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

                E_z_tmp[index] = sin(2*M_PI*900E6*current_time);
            }

            #pragma omp for schedule(static)
            for(size_t it = 0 ; it < local_nodes_inside_source_NUMBER[1].size() ; it ++){

                index = local_nodes_inside_source_NUMBER[1][it];

                E_y_tmp[index] = 0;
            }

            #pragma omp for schedule(static)
            for(size_t it = 0 ; it < local_nodes_inside_source_NUMBER[0].size() ; it ++){

                index = local_nodes_inside_source_NUMBER[0][it];

                E_x_tmp[index] = 0;
            }

            #pragma omp barrier

            /////////////////////////
            /// MPI COMMUNICATION ///
            /////////////////////////
            #pragma omp single
            {
                /// Synchronize MPI processes:
                MPI_Barrier(MPI_COMM_WORLD);
            }
            if(omp_get_thread_num() >=0 
                && omp_get_thread_num() <= 5
                && OMP_thread_has_neighboor == true
                && true)
            {
                /// Call the MPI communication function:
                MPI_communication_with_neighboors(
                    omp_get_thread_num(), /*int  omp_thread_id,*/
                    grid.MPI_communicator.getRank(),
                    grid.MPI_communicator.RankNeighbour[omp_get_thread_num()],
                    direction, /*char direction,*/
                    ElectricFieldToSend, /*double *ElectricFieldToSend,*/
                    size_to_be_sent, /*size_t size_to_send,*/
                    ElectricFieldToRecv, /*double *ElectricFieldToRecv,*/
                    size_to_recv, /*size_t size_to_recv,*/
                    E_x_tmp, /*double *E_x,*/
                    E_y_tmp, /*double *E_y,*/
                    E_z_tmp, /*double *E_z,*/
                    omp_sizes /*std::vector<size_t> &omp_sizes*/
                );
            }
            /////////////////////////
            ///      END OF       ///
            /// MPI COMMUNICATION ///
            /////////////////////////
            
            
            #pragma omp single
            {
                if(currentStep == 0){
                    grid.profiler.addTimingInputToDictionnary("ELECTRO_WRITING_OUTPUTS");
                }

                current_time += dt;
                currentStep ++;

                // Monitor time spent in writing !
                double timeWriting = omp_get_wtime();
                interfaceParaview.convertAndWriteData(
                    currentStep,
                    "ELECTRO"
                );
                double time_taken_writing = omp_get_wtime()-timeWriting;
                
                // Add time for writing:
                grid.profiler.incrementTimingInput("ELECTRO_WRITING_OUTPUTS",time_taken_writing);


                if(currentStep > 100){
                    printf("\n\t%zu ITER ABORTING\n\n",currentStep);
                    abort();
                }

                double elapsedTot = omp_get_wtime() - parallelRegionStartingTime;

                printf("AlgoElectro_NEW : iter %zu, time %.10f, time per iter : %.10f (using %d OMP and %d MPI / on MPI %d).\n",
                        currentStep,current_time,elapsedTot/currentStep,omp_get_num_threads(),
                        grid.MPI_communicator.getNumberOfMPIProcesses(),
                        grid.MPI_communicator.getRank());
                printf("Time for writing : %f seconds.\n",time_taken_writing);
            }
            #pragma omp barrier
            #pragma omp master
            {
                /*if(currentStep > 6){
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Abort(MPI_COMM_WORLD,-1);
                }*/
            }

        } /* END OF WHILE */

        /// Free memory of each thread:
        if(ElectricFieldToRecv != NULL){
            delete[] ElectricFieldToRecv;
        }
        if(ElectricFieldToSend != NULL){
            delete[] ElectricFieldToSend;
        }

    }/* END OF PARALLEL REGION */


    /* FREE MEMORY */

    for(int i = 0 ; i < 6 ; i ++){
        grid.profiler.removeMemoryUsage(
            "BYTES",
            sizeof(size_t)*local_nodes_inside_source_NUMBER[i].size(),
            "Free_nodes_inside_source_NUMBER[]");
        grid.profiler.removeMemoryUsage(
            "BYTES",
            sizeof(unsigned char)*ID_Source[i].size(),
            "Free_ID_Source[]"
        );
    }
    
    delete[] local_nodes_inside_source_NUMBER;
    delete[] ID_Source;
    
    // Free H_x coefficients:
    size = grid.size_Hx[0]*grid.size_Hx[1]*grid.size_Hx[2];
    delete[] C_hxh;
    delete[] C_hxe_1;
    delete[] C_hxe_2;

    memory = (8+8+8) * size;
    grid.profiler.removeMemoryUsage("BYTES",memory,"Algo_Free_H_x_coeffs");


    // Free H_y coefficents:
    size = grid.size_Hy[0]*grid.size_Hy[1]*grid.size_Hy[2];
    delete[] C_hyh;
    delete[] C_hye_1;
    delete[] C_hye_2;

    memory = (8+8+8) * size;
    grid.profiler.removeMemoryUsage("BYTES",memory,"Algo_Free_H_y_coeffs");

    // Free H_z coefficients:
    size = grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2];
    delete[] C_hzh;
    delete[] C_hze_1;
    delete[] C_hze_2;

    memory = (8+8+8) * size;
    grid.profiler.removeMemoryUsage("BYTES",memory,"Algo_Free_H_z_coeffs");

    // Free E_x coefficients:
    size = grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2];
    delete[] C_exe;
    delete[] C_exh_1;
    delete[] C_exh_2;

    memory = (8+8+8) * size;
    grid.profiler.removeMemoryUsage("BYTES",memory,"Algo_Free_E_x_coeffs");

    // Free E_y coefficients:
    size = grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2];
    delete[] C_eye;
    delete[] C_eyh_1;
    delete[] C_eyh_2;

    memory = (8+8+8) * size;
    grid.profiler.removeMemoryUsage("BYTES",memory,"Algo_Free_E_y_coeffs");

    // Free E_z coefficients:
    size = grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2];
    delete[] C_eze;
    delete[] C_ezh_1;
    delete[] C_ezh_2;

    memory = (8+8+8) * size;
    grid.profiler.removeMemoryUsage("BYTES",memory,"Algo_Free_E_z_coeffs");

    // End clock for monitoring CPU time
    end___algo_update = omp_get_wtime();
    double elapsedTimeSec = end___algo_update - start_algo_update;
    grid.profiler.incrementTimingInput("AlgoElectro_NEW_UPDATE_omp_get_wtime",elapsedTimeSec);
    std::cout << "AlgoElectro_NEW_UPDATE => Time: " << elapsedTimeSec << " s" << std::endl;
}

void AlgoElectro_NEW::check_OMP_DYNAMIC_envVar(void){
    /* SET OMP_DYNAMIC to false */
	if(const char *omp_dynamic_env = std::getenv("OMP_DYNAMIC")){
		// Already declared. Check it is false.
		if(std::strcmp(omp_dynamic_env,"false") == 0){
			printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
		}else{
			std::string set_env = "OMP_DYNAMIC=false";
			putenv(&set_env[0]);
			printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
		}
	}else{
		// OMP_DYNAMIC was not declared. Declare it.
		std::string set_env = "OMP_DYNAMIC=false";
		putenv(&set_env[0]);
		printf("OMP_DYNAMIC=%s.\n",std::getenv("OMP_DYNAMIC"));
	}
}

/**
 * @brief This function initialize some usefull variables for each OMP thread, in order to set MPI comminucation.
 */
void AlgoElectro_NEW::determine_OMP_thread_role_in_MPI_communication(
            int omp_thread_id /* omp_get_thread_num()*/,
            bool *OMP_thread_has_neighboor,
            char *direction,
            GridCreator_NEW &grid,
            std::vector<size_t> &omp_sizes,
            double **ElectricFieldToSend,
            double **ElectricFieldToRecv,
            size_t *size_of_sent_vector,
            size_t *size_of_recv_vector
        )
{
    /// We only need 6 OMP threads to handle the communication.
    if(omp_thread_id >= 0 && omp_thread_id <= 5){
        /// The omp thread is needed for MPI communication.

        /// Set the direction of communication.
        std::vector<char> DIRECTIONS = {'S','N','W','E','D','U'};

        *direction = DIRECTIONS[omp_thread_id];

        /// Determine if the MPI process has a neighboor is the direction:
        /////////////////////////////////////////////////
        /// CONVENTION FOR COMMUNICATION              ///
        /// 1) OMP thread(0) communicates with SOUTH. ///
        /// 2) OMP thread(1) communicates with NORTH. ///
        /// 3) OMP thread(2) communicates with WEST.  ///
        /// 3) OMP thread(3) communicates with EAST.  ///
        /// 4) OMP thread(4) communicates with DOWN.  ///
        /// 5) OMP thread(5) communicates with UP.    ///
        /////////////////////////////////////////////////
        if(grid.MPI_communicator.RankNeighbour[omp_thread_id] != -1){
            *OMP_thread_has_neighboor = true;
        }

        /// If the OMP thread has an associated neighboor, initialize
        /// the vectors used for sending and recieving.
        if(*OMP_thread_has_neighboor){

            /// Use omp_sizes to determine the size of those vectors:

            if(*direction == 'S' || *direction == 'N'){

                /// Determine the total size to be sent:
                *size_of_sent_vector  = 0;
                /// size += size_ex[1]*size_ex[2]:
                *size_of_sent_vector += (omp_sizes[1]-2)*(omp_sizes[2]-2);
                /// size += size_ey[1]*size_ey[2]:
                *size_of_sent_vector += (omp_sizes[4]-2)*(omp_sizes[5]-2);
                /// size += size_ez[1]*size_ez[2]:
                *size_of_sent_vector += (omp_sizes[7]-2)*(omp_sizes[8]-2);

            }else if(*direction == 'W' || *direction == 'E'){

                /// Determine the total size to be sent:
                *size_of_sent_vector  = 0;
                /// size += size_ex[0]*size_ex[2] (without neighboors)
                *size_of_sent_vector += (omp_sizes[0]-2)*(omp_sizes[2]-2);
                /// size += size_ey[0]*size_ey[2] (without neighboors)
                *size_of_sent_vector += (omp_sizes[3]-2)*(omp_sizes[5]-2);
                /// size += size_ez[0]*size_ez[2] (without neighboors)
                *size_of_sent_vector += (omp_sizes[6]-2)*(omp_sizes[8]-2);

            }else if(*direction == 'U' || *direction == 'D'){

                /// Determine the total size to be sent:
                *size_of_recv_vector  = 0;
                /// size += size_ex[0]*size_ex[1]:
                *size_of_sent_vector += (omp_sizes[0]-2)*(omp_sizes[1]-2);
                /// size += size_ey[0]*size_ey[1]:
                *size_of_sent_vector += (omp_sizes[3]-2)*(omp_sizes[4]-2);
                /// size += size_ez[0]*size_ez[1]:
                *size_of_sent_vector += (omp_sizes[6]-2)*(omp_sizes[7]-2);

            }else{

                fprintf(stderr,"In %s :: fatal error. Direction is %c, unknow to the system.\n",
                    __FUNCTION__,*direction);
                fprintf(stderr,"In file %s:%d\n",__FILE__,__LINE__);
                #ifdef MPI_COMM_WORLD
                MPI_Abort(MPI_COMM_WORLD,-1);
                #else
                abort();
                #endif

            }

            /// Increment memory usage, omp safe:
            *size_of_recv_vector = *size_of_sent_vector;
            
            *ElectricFieldToSend = new double[*size_of_sent_vector];
            *ElectricFieldToRecv = new double[*size_of_recv_vector];

            #pragma omp critical
            {
                grid.profiler.addMemoryUsage(
                    "BYTES",
                    sizeof(double)* (*size_of_sent_vector) * (*size_of_recv_vector)
                );
            }
        }

    }else{
        /// The thread has an ID larger than 5, just do nothing.
        return;
    }
}

/**
 * @brief Fills the electric field vectors with the received vector.
 */
void fill_electric_field_with_recv_vector(
        double *E_x,
        double *E_y,
        double *E_z,
        double *ElectricFieldToRecv,
        size_t size_to_recv,
        char   direction,
        std::vector<size_t> &omp_sizes,
        int mpi_me
){
    /**
     * This function takes the received vector and place it at the suited indices
     * of the electric field.
     */

    /// Depending on the direction:
    if(direction == 'W'){
        size_t index   = 0;
        size_t counter = 0;

        /// Put E_x:
        for(size_t k = 1 ; k < omp_sizes[2]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0]-1 ; i++){
                index      = i + omp_sizes[0] * ( 0 + omp_sizes[1] * k );
                
                E_x[index] = ElectricFieldToRecv[3*counter];

                counter++;
            }
        }
        /// Put E_y:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+3]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0+3]-1 ; i++){
                index = i + omp_sizes[0+3] * ( 0 + omp_sizes[1+3] * k );

                E_y[index] = ElectricFieldToRecv[3*counter+1];

                /*printf("Recv : E_y(%zu,%zu,%zu) = %f (MPI %d)\n",
                    i,omp_sizes[1+3],k,
                    E_y[index],mpi_me);*/

                counter++;
            }
        }
        /// Put E_z:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+2*3]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0+2*3]-1 ; i++){
                index = i + omp_sizes[0+2*3] * ( 0 + omp_sizes[1+2*3] * k);

                E_z[index] = ElectricFieldToRecv[3*counter+2];

                counter++;
            }
        }
        return;
    }
    if(direction == 'E'){
        size_t index   = 0;
        size_t counter = 0;

        /// Put E_x:
        for(size_t k = 1 ; k < omp_sizes[2]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0]-1 ; i++){
                /**
                 * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                 * En fait, si on avait itéré sur les j, on aurait fait:
                 *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                 *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                 */
                index      = i + omp_sizes[0] * ( omp_sizes[1]-1 + omp_sizes[1] * k );
                
                E_x[index] = ElectricFieldToRecv[3*counter];

                counter++;
            }
        }
        /// Put E_y:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+3]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0+3]-1 ; i++){
                /**
                 * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                 * En fait, si on avait itéré sur les j, on aurait fait:
                 *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                 *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                 */
                index = i + omp_sizes[0+3] * ( omp_sizes[1+3]-1 + omp_sizes[1+3] * k );

                E_y[index] = ElectricFieldToRecv[3*counter+1];

                /*printf("Recv : E_y(%zu,%zu,%zu) = %f (MPI %d)\n",
                    i,omp_sizes[1+3],k,
                    E_y[index],mpi_me);*/

                counter++;
            }
        }
        /// Put E_z:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+2*3]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0+2*3]-1 ; i++){
                /**
                 * Attention, il faut absolument faire omp_sizes[1+3]-1 !!!!!
                 * En fait, si on avait itéré sur les j, on aurait fait:
                 *      for(j = 1 ; j < omp_sizes[1+3] ; j ++)
                 *  donc j va de 1 à OMP_SIZES[1+3]-1 !!!
                 */
                index = i + omp_sizes[0+2*3] * ( omp_sizes[1+2*3]-1 + omp_sizes[1+2*3] * k);

                E_z[index] = ElectricFieldToRecv[3*counter+2];

                counter ++;
            }
        }
        return;
    }
    if(direction == 'S'){
        return;
    }
    if(direction == 'N'){
        return;
    }
    if(direction == 'U'){
        return;
    }
    if(direction == 'D'){
        return;
    }
}

/**
 * @brief Fills the vector to send with the electric field vectors:
 */
void fill_vector_to_send(
        double *E_x,
        double *E_y,
        double *E_z,
        double *ElectricFieldToSend,
        size_t size_to_send,
        char   direction,
        std::vector<size_t> &omp_sizes,
        int  mpi_me
){
    /**
     * This function takes the suited electric field components and put them inside the 
     * vector to be sent to the neighboors.
     */

    /// Depending on the direction:
    if(direction == 'W'){
        size_t index   = 0;
        size_t counter = 0;

        /// Put E_x:
        for(size_t k = 1 ; k < omp_sizes[2]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0]-1 ; i++){
                index      = i + omp_sizes[0] * ( 1 + omp_sizes[1] * k );
                
                ElectricFieldToSend[3*counter] = E_x[index];

                counter++;
            }
        }
        /// Put E_y:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+3]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0+3]-1 ; i++){
                index = i + omp_sizes[0+3] * ( 1 + omp_sizes[1+3] * k );

                ElectricFieldToSend[3*counter+1] = E_y[index];

                /*printf("Send : E_y(%zu,%zu,%zu) = %f (MPI %d)\n",
                    i,(size_t)1,k,E_y[index],mpi_me);*/

                counter++;
            }
        }
        /// Put E_z:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+2*3]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0+2*3]-1 ; i++){
                index = i + omp_sizes[0+2*3] * ( 1 + omp_sizes[1+2*3] * k);

                ElectricFieldToSend[3*counter+2] = E_z[index];

                counter++;
            }
        }
        return;
    }
    if(direction == 'E'){
        size_t index   = 0;
        size_t counter = 0;

        /// Put E_x:
        for(size_t k = 1 ; k < omp_sizes[2]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0]-1 ; i++){
                /**
                 * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                 * En fait, si on avait itér sur les j, on aurait fait:
                 *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                 *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                 */
                index      = i + omp_sizes[0] * ( omp_sizes[1]-2 + omp_sizes[1] * k );
                
                ElectricFieldToSend[3*counter] = E_x[index];

                counter++;
            }
        }
        /// Put E_y:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+3]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0+3]-1 ; i++){
                /**
                 * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                 * En fait, si on avait itér sur les j, on aurait fait:
                 *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                 *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                 */
                index = i + omp_sizes[0+3] * ( omp_sizes[1+3]-2 + omp_sizes[1+3] * k );

                ElectricFieldToSend[3*counter+1] = E_y[index];

                /*printf("SEND : E_y(%zu,%zu,%zu) = %f (MPI %d) - OMP_SIZES(1+3) = %zu\n",
                    i,omp_sizes[1+3],k,
                    E_y[index],mpi_me,
                    omp_sizes[1+3]);*/

                counter++;
            }
        }
        /// Put E_z:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+2*3]-1 ; k++){
            for(size_t i = 1 ; i < omp_sizes[0+2*3]-1 ; i++){
                /**
                 * Attention, il faut absolument faire omp_sizes[1+3]-2 !!!!!
                 * En fait, si on avait itér sur les j, on aurait fait:
                 *      for(j = 1 ; j < omp_sizes[1+3]-1 ; j ++)
                 *  donc j va de 1 à OMP_SIZES[1+3]-2 !!!
                 */
                index = i + omp_sizes[0+2*3] * ( omp_sizes[1+2*3]-2 + omp_sizes[1+2*3] * k);

                ElectricFieldToSend[3*counter+2] = E_z[index];

                counter ++;
            }
        }
        return;
    }
    if(direction == 'S'){
        size_t index   = 0;
        size_t counter = 0;

        /// Put E_x:
        for(size_t k = 1 ; k < omp_sizes[2]-1 ; k++){
            for(size_t j = 1 ; j < omp_sizes[1]-1 ; j++){
                index      = omp_sizes[0]-2 + omp_sizes[0] * ( j + omp_sizes[1] * k );
                
                ElectricFieldToSend[3*counter] = E_x[index];

                counter++;
            }
        }
        /// Put E_y:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+3]-1 ; k++){
            for(size_t j = 1 ; j < omp_sizes[1+3]-1 ; j++){
                index = omp_sizes[0+3]-2 + omp_sizes[0+3] * ( j + omp_sizes[1+3] * k );

                ElectricFieldToSend[3*counter+1] = E_y[index];

                /*printf("Send : E_y(%zu,%zu,%zu) = %f (MPI %d)\n",
                    i,(size_t)1,k,E_y[index],mpi_me);*/

                counter++;
            }
        }
        /// Put E_z:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+2*3]-1 ; k++){
            for(size_t j = 1 ; j < omp_sizes[1+2*3]-1 ; j++){
                index = omp_sizes[0+2*3]-2 + omp_sizes[0+2*3] * ( j + omp_sizes[1+2*3] * k);

                ElectricFieldToSend[3*counter+2] = E_z[index];

                counter++;
            }
        }
        return;
    }
    if(direction == 'N'){
        size_t index   = 0;
        size_t counter = 0;

        /// Put E_x:
        for(size_t k = 1 ; k < omp_sizes[2]-1 ; k++){
            for(size_t j = 1 ; j < omp_sizes[1]-1 ; j++){
                index      = 1 + omp_sizes[0] * ( j + omp_sizes[1] * k );
                
                ElectricFieldToSend[3*counter] = E_x[index];

                counter++;
            }
        }
        /// Put E_y:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+3]-1 ; k++){
            for(size_t j = 1 ; j < omp_sizes[1+3]-1 ; j++){
                index = 1 + omp_sizes[0+3] * ( j + omp_sizes[1+3] * k );

                ElectricFieldToSend[3*counter+1] = E_y[index];

                /*printf("Send : E_y(%zu,%zu,%zu) = %f (MPI %d)\n",
                    i,(size_t)1,k,E_y[index],mpi_me);*/

                counter++;
            }
        }
        /// Put E_z:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+2*3]-1 ; k++){
            for(size_t j = 1 ; j < omp_sizes[1+2*3]-1 ; j++){
                index = 1 + omp_sizes[0+2*3] * ( j + omp_sizes[1+2*3] * k);

                ElectricFieldToSend[3*counter+2] = E_z[index];

                counter++;
            }
        }
        return;
    }
    if(direction == 'U'){
        size_t index   = 0;
        size_t counter = 0;

        /// Put E_x:
        for(size_t j = 1 ; j < omp_sizes[2]-1 ; j++){
            for(size_t i = 1 ; i < omp_sizes[1]-1 ; i++){
                index      = omp_sizes[0]-2 + omp_sizes[0] * ( j + omp_sizes[1] * k );
                
                ElectricFieldToSend[3*counter] = E_x[index];

                counter++;
            }
        }
        /// Put E_y:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+3]-1 ; k++){
            for(size_t j = 1 ; j < omp_sizes[1+3]-1 ; j++){
                index = omp_sizes[0+3]-2 + omp_sizes[0+3] * ( j + omp_sizes[1+3] * k );

                ElectricFieldToSend[3*counter+1] = E_y[index];

                /*printf("Send : E_y(%zu,%zu,%zu) = %f (MPI %d)\n",
                    i,(size_t)1,k,E_y[index],mpi_me);*/

                counter++;
            }
        }
        /// Put E_z:
        counter = 0;
        for(size_t k = 1 ; k < omp_sizes[2+2*3]-1 ; k++){
            for(size_t j = 1 ; j < omp_sizes[0+2*3]-1 ; j++){
                index = omp_sizes[0+2*3]-2 + omp_sizes[0+2*3] * ( j + omp_sizes[1+2*3] * k);

                ElectricFieldToSend[3*counter+2] = E_z[index];

                counter++;
            }
        }
        return;
    }
    if(direction == 'D'){
        return;
    }

}

void MPI_communication_custom_function(
        size_t size_to_send,
        double *ElectricFieldToSend,
        size_t size_to_recv,
        double *ElectricFieldToRecv,
        int    mpi_from_to_who,
        int    mpi_me,
        int    omp_thread_id
    )
{
    /// Determine the omp thread ID with which we'll communicate:
    int omp_Neighboor = INT_MIN;

    if(omp_thread_id == 0){
		omp_Neighboor=1;

	}else if(omp_thread_id == 1){
		omp_Neighboor = 0;

	}else if(omp_thread_id == 2){
		omp_Neighboor =3;

	}else if(omp_thread_id == 3){
		omp_Neighboor = 2;

	}else if(omp_thread_id == 4){
		omp_Neighboor = 5;

	}else if(omp_thread_id == 5){
		omp_Neighboor = 4;

	}else{
		printf("In %s ::ERROR\n",__FUNCTION__);
		printf("\tIn MPI %d, at line %d, in file %s\n",mpi_me,
			__LINE__,__FILE__);
		std::abort();
	}

    /// Do som MPI send/recv in a logic order:

    /* NEIGHBOOR IS EVEN, I AM ODD */
	if(( mpi_from_to_who%2 ) == 0 &&
	   ( mpi_me%2          ) != 0){
		   	printf(">>>>>>>>>>>>>>>>>>>NeigEven :: OMP"
               " %d , RANK(%d)=%d, MYRANK=%d\n",
                omp_thread_id,omp_thread_id,
				mpi_from_to_who,
				mpi_me);

			/* We first send and then receive */
			MPI_Send(
                ElectricFieldToSend,
                size_to_send,
				MPI_DOUBLE,
				mpi_from_to_who,
                omp_thread_id,
				MPI_COMM_WORLD);

			MPI_Recv(
                ElectricFieldToRecv,
                size_to_recv,
				MPI_DOUBLE,
                mpi_from_to_who,
				omp_Neighboor,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);

	/* NEIGHBOOR IS ODD, I AM EVEN */
	}else if( 
        ( mpi_from_to_who%2 ) != 0 &&
		( mpi_me%2 )          == 0)
    {
			printf(">>>>>>>>>>>>>>>>>>>NeighOdd :: OMP"
            " %d , RANK(%d)=%d, MYRANK=%d \n",
                omp_thread_id,
                omp_thread_id,
				mpi_from_to_who,
				mpi_me);
                		
        	/* We first receive and then send */

			MPI_Recv(
                ElectricFieldToRecv,
                size_to_recv,
				MPI_DOUBLE,
                mpi_from_to_who,
				omp_Neighboor,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);

			MPI_Send(
                ElectricFieldToSend,
                size_to_send,
                MPI_DOUBLE,
				mpi_from_to_who,
                omp_thread_id,
				MPI_COMM_WORLD);

	/* NEIGHBOOR LARGER THAN ME */
	}else if(
        mpi_from_to_who > mpi_me){

			printf("Larger :: OMP %d , RANK(%d)=%d, MYRANK=%d\n",
                omp_thread_id,omp_thread_id,
				mpi_from_to_who,
				mpi_me);

			/* First receive and then send */
			
			MPI_Recv(
                ElectricFieldToRecv,
                size_to_recv,
				MPI_DOUBLE,
                mpi_from_to_who,
				omp_Neighboor,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);

			MPI_Send(
                ElectricFieldToSend,
                size_to_send,
				MPI_DOUBLE,
				mpi_from_to_who,
                omp_thread_id,
				MPI_COMM_WORLD);

	/* NEIGHBOOR SMALLER THAN ME */
	}else if( mpi_from_to_who < mpi_me ){

			printf("Smaller :: OMP %d , RANK(%d)=%d, MYRANK=%d\n",
                omp_thread_id,omp_thread_id,
				mpi_from_to_who,
				mpi_me);
			/* First send and then receive */
			
			MPI_Send(
                ElectricFieldToSend,
                size_to_send,
                MPI_DOUBLE,
				mpi_from_to_who,
                omp_thread_id,
				MPI_COMM_WORLD);

			MPI_Recv(
                ElectricFieldToRecv,
                size_to_recv,
				MPI_DOUBLE,
                mpi_from_to_who,
				omp_Neighboor,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);

	}else{
		printf("comminucateWithNeighboor not yet implemented, abort!\n");
		std::abort();
	}

}

void MPI_communication_with_neighboors(
    int  omp_thread_id,
    int  mpi_me,
    int  mpi_from_to_who,
    char direction,
    double *ElectricFieldToSend,
    size_t size_to_send,
    double *ElectricFieldToRecv,
    size_t size_to_recv,
    double *E_x,
    double *E_y,
    double *E_z,
    std::vector<size_t> &omp_sizes
){
    /// Check that the OMP thread ID is in the good range [0-5]:
    if(omp_thread_id < 0 || omp_thread_id > 5){
        fprintf(stderr,"In %s :: an OMP thread with ID not "
            "inside the range [0-5] should not end up here (has %d). Aborting.\n",
            __FUNCTION__,omp_thread_id);
        fprintf(stderr,"In file %s:%d\n",__FILE__,__LINE__);
        #ifdef MPI_COMM_WORLD
        MPI_Abort(MPI_COMM_WORLD,-1);
        #else
        abort();
        #endif
    }

    /// Verify the vectors send/recv are not null pointers:
    if( ElectricFieldToRecv == NULL || ElectricFieldToSend == NULL){
        fprintf(stderr,"In %s :: The vector to recv/send is a NULL pointer. Aborting.\n",
            __FUNCTION__);
        fprintf(stderr,"In file %s:%d\n",__FILE__,__LINE__);
        #ifdef MPI_COMM_WORLD
        MPI_Abort(MPI_COMM_WORLD,-1);
        #else
        abort();
        #endif
    }

    /// Verify the electric field pointers:
    if( E_x == NULL || E_y == NULL || E_z == NULL){
        fprintf(stderr,"In %s :: The electric field vectors are NULL. Aborting.\n",
            __FUNCTION__);
        fprintf(stderr,"In file %s:%d\n",__FILE__,__LINE__);
        #ifdef MPI_COMM_WORLD
        MPI_Abort(MPI_COMM_WORLD,-1);
        #else
        abort();
        #endif
    }

    /// Fill in the vector to send:
    printf("\n\n>>> MPI %d sends to %c (MPI %d)\n\n",
        mpi_me,direction,mpi_from_to_who);

    
    

    fill_vector_to_send(
        E_x,
        E_y,
        E_z,
        ElectricFieldToSend,
        size_to_send,
        direction,
        omp_sizes,
        mpi_me);

    /// Perform the communication:
    MPI_communication_custom_function(
        size_to_send,
        ElectricFieldToSend,
        size_to_recv,
        ElectricFieldToRecv,
        mpi_from_to_who,
        mpi_me,
        omp_thread_id
    );

    /// Fill in the electric field vectors with the received vector:
    fill_electric_field_with_recv_vector(
        E_x,
        E_y,
        E_z,
        ElectricFieldToRecv,
        size_to_recv,
        direction,
        omp_sizes,
        mpi_me
    );

}