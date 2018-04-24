#include "AlgoElectro_NEW.hpp"

#include <ctime>
#include <new>
#include "omp.h"
#include "mpi.h"
#include <algorithm>
#include <sys/time.h>

#include "UTILS/vector_utilities.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <limits.h>
#include <unistd.h>

#include <cmath>

#include "DISC_INTEGR/discrete_integration_util.hpp"

#include <assert.h>

#include "header_with_all_defines.hpp"

#define NBR_FACES_CUBE 6

#define DECALAGE_E_SUPP 1

#include <sys/file.h>
 #define   LOCK_SH   1    /* shared lock */
 #define   LOCK_EX   2    /* exclusive lock */
 #define   LOCK_NB   4    /* don't block when locking */
 #define   LOCK_UN   8    /* unlock */

void fflush_stdout(void){fflush(stdout);}

void probe_a_field(
    GridCreator_NEW &grid,
    std::string &which_field,
    std::string &filename,
    std::string &which_form_to_probe,
    std::vector<double> &infoOnForm,
    std::vector<double> &electro_deltas,
    double current_time,
    double dt
);


int tryGetLock( char const *lockName );
void releaseLock( int fd);

void testlock(void) {
  # pragma omp parallel num_threads(16)
  {    
    int fd = -1; char ln[] = "testlock.lock";
    while (fd == -1) fd = tryGetLock(ln);

    cout << omp_get_thread_num() << ": got the lock!\n";
    cout << omp_get_thread_num() << ": removing the lock\n";
    FILE *file = fopen(ln,"a");
    if( file != NULL){
        fprintf(file,"Coucou de OMP %d.\n",omp_get_thread_num());
        fclose(file);
    }else{
        printf("Fail opening file !\n");
        abort();
    }

    system("cat testlock.lock");

    releaseLock(fd);
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
    double dt = DBL_MAX;
    // Temporary variable:
    double tmp= 0.0;
    // Iterator:
    unsigned char i=0;                                      

    // Iterate over the number of materials. For each material, compute the required time step.
    // At the end, the smallest time step is chosen. 
	size_t nbr_mat = mesh.materials.unified_material_list.size();
    for (i = 0 ; i < nbr_mat ; i++ ){
			/// Get material name:
			std::string mat_name = mesh.materials.unified_material_list[i].name;
			/// Get permeability:
			double rel_permeability = 0;
			std::map<std::string,double>::iterator it;
			it = mesh.materials.unified_material_list[i].properties.find("RELATIVEPERMEABILITY");
			if(it == mesh.materials.unified_material_list[i].properties.end()){
				/// Use default relative permeability
				rel_permeability = 1;
			}else{
				rel_permeability = mesh
									.materials
									.unified_material_list[i].properties["RELATIVEPERMEABILITY"];
				if(rel_permeability == nan("") || rel_permeability == 0)
					continue;
			}
            double mu_material = rel_permeability * VACUUM_PERMEABILITY;

            // Get permittivity:
			double rel_permittivity = 0;
			it = mesh.materials.unified_material_list[i].properties.find("RELATIVEPERMITTIVITY");
			if(it == mesh.materials.unified_material_list[i].properties.end()){
				/// Use default relative permeability
				rel_permittivity = 1;
				DISPLAY_WARNING("using def rel per");
			}else{
				rel_permittivity = mesh
									.materials
									.unified_material_list[i].properties["RELATIVEPERMITTIVITY"];
				if(rel_permittivity == nan("") || rel_permittivity == 0)
					continue;
			}
            double epsilon_material = rel_permittivity * VACUUM_PERMITTIVITY;
			
			/*printf("%20s - [mu::%.9g | eps::%.9g | rel_mu::%.9g | rel_eps::%.9g]\n",
					mat_name.c_str(),
					mu_material,
					epsilon_material,
					rel_permeability,
					rel_permittivity);*/
			
            // Compute speed of light:
            double c = 1/(sqrt(mu_material*epsilon_material));
            // Take the smallest time step:

			tmp = 1/(c*sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));
			if( tmp < dt )
				dt = tmp;
    }
	/*printf(">>> Chosen dt is %.9g.\n",dt);*/
	
	if(dt == DBL_MAX){
		DISPLAY_ERROR_ABORT(
			"There was an error because dt == DBL_MAX is true."
		);
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
    //std::cout << "AlgoElectro_NEW :: dt is " << dt << std::endl;

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
	if(this->VERBOSITY >= 1)
		printf("\t> [MPI %d] - Initializing electromagnetic properties...\n",
				grid.MPI_communicator.getRank());
    if(grid.input_parser.get_SimulationType() == "USE_AIR_EVERYWHERE"){
		DISPLAY_ERROR_ABORT("USE_AIR_EVERYWHERE is depreciated.");
        grid.Initialize_Electromagnetic_Properties("AIR_AT_INIT_TEMP");
    }else{
        grid.Initialize_Electromagnetic_Properties("INIT_TEMP");
    }

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

    fill_double_vector_with_zeros(C_hxh,size);
    fill_double_vector_with_zeros(C_hxe_1,size);
    fill_double_vector_with_zeros(C_hxe_2,size);

    // Magnetic field Hy:
    size = grid.size_Hy[0] * grid.size_Hy[1] * grid.size_Hy[2];
    double *C_hyh   = new double[size]();
    double *C_hye_1 = new double[size]();
    double *C_hye_2 = new double[size]();

    fill_double_vector_with_zeros(C_hyh,size);
    fill_double_vector_with_zeros(C_hye_1,size);
    fill_double_vector_with_zeros(C_hye_2,size);

    // Magnetic field Hz:
    size = grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2];
    double *C_hzh   = new double[size]();
    double *C_hze_1 = new double[size]();
    double *C_hze_2 = new double[size]();

    fill_double_vector_with_zeros(C_hzh,size);
    fill_double_vector_with_zeros(C_hze_1,size);
    fill_double_vector_with_zeros(C_hze_2,size);

    // Electric field Ex:
    size = grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2];
    double *C_exe   = new double[size]();
    double *C_exh_1 = new double[size]();
    double *C_exh_2 = new double[size]();

    fill_double_vector_with_zeros(C_exe,size);
    fill_double_vector_with_zeros(C_exh_1,size);
    fill_double_vector_with_zeros(C_exh_2,size);

    // Electric field Ey:
    size = grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]; 
    double *C_eye   = new double[size]();
    double *C_eyh_1 = new double[size]();
    double *C_eyh_2 = new double[size]();

    fill_double_vector_with_zeros(C_eye,size);
    fill_double_vector_with_zeros(C_eyh_1,size);
    fill_double_vector_with_zeros(C_eyh_2,size);

    // Electric field Ez:
    size = grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2];
    double *C_eze   = new double[size]();
    double *C_ezh_1 = new double[size]();
    double *C_ezh_2 = new double[size]();

    fill_double_vector_with_zeros(C_eze,size);
    fill_double_vector_with_zeros(C_ezh_1,size);
    fill_double_vector_with_zeros(C_ezh_2,size);


    /**
     * @brief Allocation of memory for ABC conditions:
     */

    // ABC Old Tangential Field Ey at the extrmities of x of the grid:
    size = (grid.size_Ey[1]-2)*(grid.size_Ey[2]-2);
    double *Eyx0    = NULL;
    double *Eyx1    = NULL;
    Eyx0 = new double[size]();   
    Eyx1 = new double[size]();
    //std::fill_n(Eyx0, size, 0);
    //std::fill_n(Eyx1, size, 0);

    // std::vector<double> Eyx0(size, 0.0);
    // std::vector<double> Eyx1(size, 0.0);

    // ABC Old Tangential Field Ez at the extrmities of x of the grid:
    size = (grid.size_Ez[1]-2)*(grid.size_Ez[2]-2);
    double *Ezx0    = NULL;
    double *Ezx1    = NULL;
    Ezx0 = new double[size]();  
    Ezx1 = new double[size]();
    //std::fill_n(Ezx0, size, 0);
    //std::fill_n(Ezx1, size, 0);

    // std::vector<double> Ezx0(size, 0.0);
    // std::vector<double> Ezx1(size, 0.0);

    // ABC Old Tangential Field Ex at the extrmities of y of the grid:
    size = (grid.size_Ex[0]-2)*(grid.size_Ex[2]-2);
    double *Exy0    = NULL;
    double *Exy1    = NULL;
    Exy0 = new double[size]();    
    Exy1 = new double[size]();
    //std::fill_n(Exy0, size, 0);
    //std::fill_n(Exy1, size, 0);

    // std::vector<double> Exy0(size, 0.0);
    // std::vector<double> Exy1(size, 0.0);

    // ABC Old Tangential Field Ez at the extrmities of y of the grid:
    size = (grid.size_Ez[0]-2)*(grid.size_Ez[2]-2);
    double *Ezy0    = NULL;    
    double *Ezy1    = NULL;
    Ezy0 = new double[size]();
    Ezy1 = new double[size]();
    //std::fill_n(Ezy0, size, 0);
    //std::fill_n(Ezy1, size, 0);
    // std::vector<double> Ezy0(size, 0.0);
    // std::vector<double> Ezy1(size, 0.0);

    // ABC Old Tangential Field Ex at the extrmities of z of the grid:
    size = (grid.size_Ex[0]-2)*(grid.size_Ex[1]-2);
    double *Exz0    = NULL;
    double *Exz1    = NULL;
    Exz0 = new double[size]();
    Exz1 = new double[size]();
    //std::fill_n(Exz0, size, 0);
    //std::fill_n(Exz1, size, 0);
    
    // std::vector<double> Exz0(size, 0.0);
    // std::vector<double> Exz1(size, 0.0);

    // ABC Old Tangential Field Ey at the extrmities of z of the grid:
    size = (grid.size_Ey[0]-2)*(grid.size_Ey[1]-2);
    double *Eyz0    = NULL;
    double *Eyz1    = NULL;
    Eyz0 = new double[size]();    
    Eyz1 = new double[size]();    
    //std::fill_n(Eyz0, size, 0);
    //std::fill_n(Eyz1, size, 0);
    // std::vector<double> Eyz0(size, 0.0);
    // std::vector<double> Eyz1(size, 0.0);







	if(this->VERBOSITY >= 1)
		printf("\t> [MPI %d] - Computing coefficients...\n",grid.MPI_communicator.getRank());
    /* COMPUTING COEFFICIENTS */
    #pragma omp parallel default(none)\
		firstprivate(C_exe,C_exh_1,C_exh_2)\
		firstprivate(C_eye,C_eyh_1,C_eyh_2)\
		firstprivate(C_eze,C_ezh_1,C_ezh_2)\
		firstprivate(C_hxh,C_hxe_1,C_hxe_2)\
		firstprivate(C_hyh,C_hye_1,C_hye_2)\
		firstprivate(C_hzh,C_hze_1,C_hze_2)\
		shared(grid)\
		firstprivate(dt)
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
    
    // Size 3 because 3 E components:
    std::vector<size_t>        *local_nodes_inside_source_NUMBER ;
    local_nodes_inside_source_NUMBER = new std::vector<size_t>[3];

    // Size 3 because 3 E components:
    std::vector<unsigned char> *ID_Source                         ;
    ID_Source = new std::vector<unsigned char>[3];

    std::vector<double>        local_nodes_inside_source_FREQ    ;

    std::vector<std::string> TYPE = {"Ex","Ey","Ez"};
	
	if(this->VERBOSITY >= 1){
		fflush_stdout();
		MPI_Barrier(MPI_COMM_WORLD);
		printf("\t> [MPI %d] - Computing nodes inside sources...\n",
				grid.MPI_communicator.getRank());
		fflush_stdout();
	}
	
    for(unsigned int i = 0 ; i < TYPE.size() ; i ++){
        grid.Compute_nodes_inside_sources(
            local_nodes_inside_source_NUMBER[i],
            ID_Source[i],
            local_nodes_inside_source_FREQ,
            TYPE[i]
        );
        if(local_nodes_inside_source_NUMBER[i].size() != ID_Source[i].size()){
            DISPLAY_ERROR_ABORT(
				"Sizes do not match (has %zu and %zu).",
				local_nodes_inside_source_NUMBER[i].size(),
				ID_Source[i].size()
			);
        }
    }
	/// Verify that there is at least one emitting element:
	int at_least_one_node = 0;
	for(size_t I = 0 ; I < grid.input_parser.source.get_number_of_sources() ; I++){
		if(grid.input_parser.source.there_is_at_least_one_element_non_zero_in_source[I]){
			at_least_one_node = 1;
			printf("[MPI %d over %d] - There is at least one node !\n",
				grid.MPI_communicator.getRank(),
				grid.MPI_communicator.getNumberOfMPIProcesses());
			break;
		}
	}
	/// Communicate between MPI processes to check if there is a source somewhere:
	int *checking_at_least_one_node 
		= (int*)malloc(sizeof(int)*grid.MPI_communicator.getNumberOfMPIProcesses());
	MPI_Gather(
		&at_least_one_node,
		1,
		MPI_INT,
		(int*)checking_at_least_one_node,
		1,//grid.MPI_communicator.getNumberOfMPIProcesses(),
		MPI_INT,
		grid.MPI_communicator.rootProcess,
		MPI_COMM_WORLD);
	int counter_false = 0;
	if(grid.MPI_communicator.getRank() == grid.MPI_communicator.rootProcess){
		for(int i = 0 ; i < grid.MPI_communicator.getNumberOfMPIProcesses() ; i++){
			if(checking_at_least_one_node[i] == 0){
				counter_false++;
			}
		}
		if(counter_false == grid.MPI_communicator.getNumberOfMPIProcesses()){
			DISPLAY_ERROR_ABORT(
				"There is no node emitting anything. Your solution will remain zero everywhere."
			);
		}else{
			DISPLAY_WARNING("Ok, there is at least on node inside your sources.");
		}
	}
	free(checking_at_least_one_node);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	if(this->VERBOSITY >= 1){
		fflush_stdout();
		MPI_Barrier(MPI_COMM_WORLD);
		printf("\t> [MPI %d] - Computing nodes inside sources...\n",
				grid.MPI_communicator.getRank());
		fflush_stdout();
	}

    /// Assign frequencies:
    local_nodes_inside_source_FREQ = grid.input_parser.source.frequency;
    
    /// Clean the output:
    if(grid.MPI_communicator.isRootProcess() == 0)
        {
			fflush(stdout);
            printf(">>> FDTD scheme started with time step of %.15lf seconds.\n",
                dt);
			fflush_stdout();
        }

    ///////////////////////////////////////////
    /// VARIABLES FOR STEADY-STATE CHECKING ///
    ///////////////////////////////////////////
    /* For the Ex field */
    size_t checkEvery = grid.input_parser.SteadyState_CheckEveryPoint;
    size_t nbr_pointsX = grid.size_Ex[0]/checkEvery;
    size_t nbr_pointsY = grid.size_Ex[1]/checkEvery;
    size_t nbr_pointsZ = grid.size_Ex[2]/checkEvery;
    size_t nbr_points  = nbr_pointsX * nbr_pointsY * nbr_pointsZ;
    std::vector<double> cumtrapz(nbr_points);
    std::vector<double> cummmean(nbr_points);
    std::vector<double> derivati(nbr_points);
    DISPLAY_WARNING("Steadiness : en cours de construction.");


    ////////////////////////////////////////////////
    // UPDATE WHILE LOOP - PARALLELIZED WITH      //
    // OPENMP THREADS    - MINIMUM 6 OPENMP       //
    // THREADS ARE REQUIRED FOR MPI COMMUNICATION //
    // TO WORK.                                   //
    //////////////////////////////////////////////// 

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
        firstprivate(local_nodes_inside_source_NUMBER)\
        firstprivate(local_nodes_inside_source_FREQ,ID_Source)\
        shared(interfaceParaview)\
        firstprivate(C_hxh,C_hxe_1,C_hxe_2)\
        firstprivate(C_hyh,C_hye_1,C_hye_2)\
        firstprivate(C_hzh,C_hze_1,C_hze_2)\
        firstprivate(C_exe,C_exh_1,C_exh_2)\
        firstprivate(C_eye,C_eyh_1,C_eyh_2)\
        firstprivate(C_eze,C_ezh_1,C_ezh_2)\
        shared(ompi_mpi_comm_world,ompi_mpi_int)\
        firstprivate(Electric_field_to_send,Electric_field_to_recv)\
        firstprivate(Magnetic_field_to_send,Magnetic_field_to_recv)\
        firstprivate(electric_field_sizes,magnetic_field_sizes,dt)\
        firstprivate(size_faces_electric,size_faces_magnetic)\
        firstprivate(Eyx0, Eyx1)\
        firstprivate(Ezx0, Ezx1)\
        firstprivate(Exy0, Exy1)\
        firstprivate(Ezy0, Ezy1)\
        firstprivate(Exz0, Exz1)\
        firstprivate(Eyz0, Eyz1)
    {
        /*#pragma omp master
        {
            fflush_stdout();
            MPI_Barrier(MPI_COMM_WORLD);
            printf("\t >>>> MPI %d enters the parallel region.\n",grid.MPI_communicator.getRank());
            fflush_stdout();
        }
        #pragma omp barrier*/

		double MIN_GAUSS_BEFORE_LET_BE = 1E-100;
        std::vector<bool> MODULATE_SOURCE(grid.input_parser.source.get_number_of_sources());
		for(size_t i = 0 ; i < MODULATE_SOURCE.size() ; i++){
			if(grid.input_parser.source_time[i] == "GAUSSIAN"){
				MODULATE_SOURCE[i]      = true;
				MIN_GAUSS_BEFORE_LET_BE = 1E-5;
			}else{
				MODULATE_SOURCE[i] = false;
			}
		}

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

        bool has_neighboor = false;
        for(unsigned int ii =  0; ii < NBR_FACES_CUBE ; ii ++){
            if(grid.MPI_communicator.RankNeighbour[ii] != -1){
                has_neighboor = true;
                break;
            }
        }

        // Some indexing variables:
        size_t I,J,K;

        /// Variables to monitor the time spent communicating:
        struct timeval start_mpi_comm;
        struct timeval end___mpi_comm;
        double         total_mpi_comm = 0.0;

        /// Variables to compute the time taken by each iteration:
        struct timeval start_while_iter;
        struct timeval end___while_iter;
        double         total_while_iter = 0.0;

        while(current_time < grid.input_parser.get_stopTime()
                && currentStep < grid.input_parser.maxStepsForOneCycleOfElectro){

            /*#pragma omp master
            {
                fflush_stdout();
                MPI_Barrier(MPI_COMM_WORLD);
                printf("\t >>> MPI %d enters the while loop.\n",grid.MPI_communicator.getRank());
                fflush_stdout();
            }
            #pragma omp barrier*/

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
                {
                    if(grid.MPI_communicator.getRank() == 0)
                    printf("%s>>> %s!!! WARNING !!!%s 'NDEBUG' is not defined. You are in debug mode."
                        " Be aware that the code is subsequently much slower. As an example,"
                        " a lot of asserts and printf's are performed.%s\n",
                            ANSI_COLOR_RED,
                            ANSI_COLOR_YELLOW,
                            ANSI_COLOR_GREEN,
                            ANSI_COLOR_RESET);
                }
            #endif

            #pragma omp for schedule(static) collapse(3) nowait
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
            /*#pragma omp master
            {
                fflush_stdout();
                MPI_Barrier(MPI_COMM_WORLD);
                printf("\t>[MPI %d] - Hx ok(EM) [step %zu].\n",
					grid.MPI_communicator.getRank(),currentStep);
            }
            #pragma omp barrier*/
			

            // Updating the magnetic field Hy.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Hy[0];
            size_y = grid.size_Hy[1];

            size_x_1 = grid.size_Ez[0];
            size_y_1 = grid.size_Ez[1];

            size_x_2 = grid.size_Ex[0];
            size_y_2 = grid.size_Ex[1];

            #pragma omp for schedule(static) collapse(3) nowait
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
			/*#pragma omp master
            {
                fflush_stdout();
                MPI_Barrier(MPI_COMM_WORLD);
                printf("\t>[MPI %d] - Hy ok(EM) [step %zu].\n",
					grid.MPI_communicator.getRank(),currentStep);
            }
            #pragma omp barrier*/

            // Updating the magnetic field Hz.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Hz[0];
            size_y = grid.size_Hz[1];

            size_x_1 = grid.size_Ex[0];
            size_y_1 = grid.size_Ex[1];

            size_x_2 = grid.size_Ey[0];
            size_y_2 = grid.size_Ey[1];

            #pragma omp for schedule(static) collapse(3) nowait
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
						
						if(false && grid.MPI_communicator.getRank() == 0){
							/*printf("Coucou Hz mpi 0 [%zu,%zu,%zu,%zu,%zu] / size[%zu,%zu,%zu]\n",
									index,
									index_1Plus,
									index_1Moins,
									index_2Plus,
									index_2Moins,
									grid.size_Hz[0]*grid.size_Hz[1]*grid.size_Hz[2],
									grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2],
									grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]);*/
                            printf("Access to H_z_tmp[%zu] %lf.\n",index,H_z_tmp[index]);
                            printf("Access to C_hzh[%zu]   %lf.\n",index,C_hzh[index]);
                            printf("Access to C_hze_1[%zu] %lf.\n",index,C_hze_1[index]);
                        }

                        H_z_tmp[index] = C_hzh[index] * H_z_tmp[index]
                                + C_hze_1[index] * (E_x_tmp[index_1Plus] - E_x_tmp[index_1Moins])
                                - C_hze_2[index] * (E_y_tmp[index_2Plus] - E_y_tmp[index_2Moins]);
                    }
                }
            }
			/*#pragma omp master
            {
                fflush_stdout();
                MPI_Barrier(MPI_COMM_WORLD);
                printf("\t>[MPI %d] - Hz ok(EM) [step %zu].\n",
					grid.MPI_communicator.getRank(),currentStep);
            }
            #pragma omp barrier*/
			

            /////////////////////////////////////////////////////
            /// OPENMP barrier because we must ensure all the ///
            /// magnetic fields have been updated.            ///
            /////////////////////////////////////////////////////
            #pragma omp barrier


            /////////////////////////
            /// MPI COMMUNICATION ///
            /////////////////////////
			/*if(this->VERBOSITY >= 3 && omp_get_thread_num() == 0 && grid.MPI_communicator.getRank() == 0){
				printf("[MPI %d] - While loop: communication of H.\n",
						grid.MPI_communicator.getRank());
			}	*/
            gettimeofday( &start_mpi_comm, NULL);
            /// Prepare the array to send:
            if(has_neighboor){
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
            }
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


            #pragma omp for schedule(static) collapse(3) nowait
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

            #pragma omp for schedule(static) collapse(3) nowait
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

            #pragma omp for schedule(static) collapse(3) nowait
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
            double frequency  = 0.0;
            double false_freq = 0.0;
            double gauss      = 0.0;
            double COEF_MEAN  = 1./3.;
            double COEF_STD   = 1./10.;

            int FIELD = 2;
            #pragma omp for schedule(static) nowait
            for(size_t it = 0 ; it < local_nodes_inside_source_NUMBER[FIELD].size() ; it ++){

                index = local_nodes_inside_source_NUMBER[FIELD][it];
                ASSERT(index,<,grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2]);

                if(ID_Source[FIELD][it] == UCHAR_MAX){
                    frequency = 0;
                    // False frequency is used to stop imposing zero.
                    false_freq = local_nodes_inside_source_FREQ[0];
                    if(MODULATE_SOURCE[0] == true){
                        double period    = 2*M_PI/false_freq;
                        double MEAN      = period*COEF_MEAN;
                        double STD       = period*COEF_STD;
                                        
                        double t = current_time;
                        
                        gauss = exp(-((t-MEAN)*(t-MEAN))/(2*STD*STD));
                        if(gauss < MIN_GAUSS_BEFORE_LET_BE){
                            // do nothing
                        }else{
                            E_z_tmp[index] = 0;//sin(2*M_PI*frequency*current_time);
                        }
                    }else{
                        E_z_tmp[index] = 0;//sin(2*M_PI*frequency*current_time);
                    }
                }else{

                    frequency = local_nodes_inside_source_FREQ[ID_Source[FIELD][it]];
                    if(ID_Source[FIELD][it] >= MODULATE_SOURCE.size()){
                        printf("Out of bound !\n%s:%d\n",__FILE__,__LINE__);
                        std::abort();
                    }
                    if(MODULATE_SOURCE[ID_Source[FIELD][it]] == true){
                        double period    = 2*M_PI/frequency;
                        double MEAN      = period*COEF_MEAN;
                        double STD       = period*COEF_STD;
                                        
                        double t = current_time;
                        
                        gauss = exp(-((t-MEAN)*(t-MEAN))/(2*STD*STD));
                        /*printf("GAUSS : %.10g | MEAN %.10g | "
                                "STD %.10g | period %.10g | frequency %.10g | "
                                "should be %.10g | coef_std %.10g.\n",
                            gauss,
                            MEAN,
                            STD,
                            period,
                            frequency,
                            period*COEF_STD,
                            COEF_STD);*/
                        if(gauss < MIN_GAUSS_BEFORE_LET_BE){
                            // do nothing
                        }else{
                            E_z_tmp[index] = gauss * sin(2*M_PI*frequency*current_time);
                            //printf("Applying.\n");
                        }
                    }else{
						E_z_tmp[index] = sin(2*M_PI*frequency*current_time);
                        //printf("Case pas modulated.\n");
					}

                }                
            }


            FIELD = 1;
            #pragma omp for schedule(static) nowait
            for(size_t it = 0 ; it < local_nodes_inside_source_NUMBER[FIELD].size() ; it ++){

                index = local_nodes_inside_source_NUMBER[FIELD][it];
                ASSERT(index,<,grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2]);

                if(ID_Source[FIELD][it] == UCHAR_MAX){
                    frequency = 0;
                    false_freq = local_nodes_inside_source_FREQ[0];
                    if(MODULATE_SOURCE[0] == true){
                        double period    = 2*M_PI/frequency;
                        double MEAN      = period*COEF_MEAN;
                        double STD       = period*COEF_STD;
                                        
                        double t = current_time;
                        
                        gauss = exp(-((t-MEAN)*(t-MEAN))/(2*STD*STD));
                        if(gauss < MIN_GAUSS_BEFORE_LET_BE){
                            // do nothing
                        }else{
                            E_y_tmp[index] = 0;//sin(2*M_PI*frequency*current_time);
                        }
                    }else{
						E_y_tmp[index] = 0;//sin(2*M_PI*frequency*current_time);
					}
                }else{

                    frequency = local_nodes_inside_source_FREQ[ID_Source[FIELD][it]];
                    if(ID_Source[FIELD][it] >= MODULATE_SOURCE.size()){
                        printf("Out of bound !\n%s:%d\n",__FILE__,__LINE__);
                        std::abort();
                    }
                    if(MODULATE_SOURCE[ID_Source[FIELD][it]] == true){
                        double period    = 2*M_PI/frequency;
                        double MEAN      = period*COEF_MEAN;
                        double STD       = period*COEF_STD;
                                        
                        double t = current_time;
                        
                        gauss = exp(-((t-MEAN)*(t-MEAN))/(2*STD*STD));
                        if(gauss < MIN_GAUSS_BEFORE_LET_BE){
                            // do nothing
                        }else{
                            E_y_tmp[index] = gauss * sin(2*M_PI*frequency*current_time);
                        }
                    }else{
						E_y_tmp[index] = sin(2*M_PI*frequency*current_time);
					}
                }
            }

            FIELD = 0;
            #pragma omp for schedule(static) nowait
            for(size_t it = 0 ; it < local_nodes_inside_source_NUMBER[FIELD].size() ; it ++){

                index = local_nodes_inside_source_NUMBER[FIELD][it];
                ASSERT(index,<,grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2]);

                if(ID_Source[FIELD][it] == UCHAR_MAX){
                    frequency = 0;
                    false_freq = local_nodes_inside_source_FREQ[0];
                    if(MODULATE_SOURCE[0] == true){
                        double period    = 2*M_PI/frequency;
                        double MEAN      = period*COEF_MEAN;
                        double STD       = period*COEF_STD;
                                        
                        double t = current_time;
                        
                        gauss = exp(-((t-MEAN)*(t-MEAN))/(2*STD*STD));

                        if(gauss < MIN_GAUSS_BEFORE_LET_BE){
                            // do nothing
                        }else{
                            E_x_tmp[index] = 0;//sin(2*M_PI*frequency*current_time);
                        }
                    }else{
						E_x_tmp[index] = 0;//sin(2*M_PI*frequency*current_time);
					}
                }else{
                    frequency = local_nodes_inside_source_FREQ[ID_Source[FIELD][it]];
                    if(ID_Source[FIELD][it] >= MODULATE_SOURCE.size()){
                        printf("Out of bound !\n%s:%d\n",__FILE__,__LINE__);
                        std::abort();
                    }
                    if(MODULATE_SOURCE[ID_Source[FIELD][it]] == true){
                        double period    = 2*M_PI/frequency;
                        double MEAN      = period*COEF_MEAN;
                        double STD       = period*COEF_STD;
                                        
                        double t = current_time;
                        
                        gauss = exp(-((t-MEAN)*(t-MEAN))/(2*STD*STD));

                        if(gauss < MIN_GAUSS_BEFORE_LET_BE){
                            // do nothing
                        }else{
                            E_x_tmp[index] = gauss * sin(2*M_PI*frequency*current_time);
                        }
                    }else{
						E_x_tmp[index] = sin(2*M_PI*frequency*current_time);
					}
                }
            }

            
            
            /////////////////////////
            /// MPI COMMUNICATION ///
            /////////////////////////

            /// Wait all OPENMP threads to be sure computations are done for this step:
            #pragma omp barrier
            gettimeofday( &start_mpi_comm , NULL);
            /// Prepare the array to send:
            if(has_neighboor){
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

            }

           
            gettimeofday( &end___mpi_comm , NULL);
            total_mpi_comm += end___mpi_comm.tv_sec  - start_mpi_comm.tv_sec + 
                                (end___mpi_comm.tv_usec - start_mpi_comm.tv_usec) / 1.e6;

            /////////////////////////
            ///      END OF       ///
            /// MPI COMMUNICATION ///
            /////////////////////////
            
            #pragma omp barrier
            #pragma omp master
            {
				/*MPI_Barrier(MPI_COMM_WORLD);
				fflush_stdout();
				MPI_Barrier(MPI_COMM_WORLD);
				printf("[MPI %d] - Starting ABC.\n",grid.MPI_communicator.getRank());*/
            this->abc(grid,
                E_x_tmp, E_y_tmp, E_z_tmp, 
                Eyx0, Ezx0, 
                Eyx1, Ezx1, 
                Exy0, Ezy0, 
                Exy1, Ezy1, 
                Exz0, Eyz0, 
                Exz1, Eyz1,
                dt
                );
				/*fflush_stdout();
				MPI_Barrier(MPI_COMM_WORLD);
				printf("[MPI %d] - Ending ABC.\n",grid.MPI_communicator.getRank());*/
            }
            #pragma omp barrier

            #pragma omp master


            gettimeofday( &end___while_iter , NULL);
            total_while_iter += end___while_iter.tv_sec  - start_while_iter.tv_sec + 
                                (end___while_iter.tv_usec - start_while_iter.tv_usec) / 1.e6;

            currentStep ++;

            

            #pragma omp master
            {
                /// PROBE POINTS IF NECESSARY
                if(!grid.input_parser.points_to_be_probed.empty()){
                    for(size_t curr_pt = 0 ; curr_pt < grid.input_parser.points_to_be_probed.size();
                                curr_pt ++)
                        {
                            std::string which_form_to_probe = "point";
                            probe_a_field(
                                //GridCreator_NEW &grid
                                grid,
                                //std::string &which_field
                                grid.input_parser.points_to_be_probed[curr_pt].type_field,
                                //std::string &filename
                                grid.input_parser.points_to_be_probed[curr_pt].filename,
                                //std::string &which_form_to_probe
                                which_form_to_probe,
                                //std::vector<double> &infoOnForm
                                grid.input_parser.points_to_be_probed[curr_pt].coordinates,
                                //std::vector<double> &electro_deltas
                                grid.delta_Electromagn,
                                //double current_time
                                current_time,
                                //double dt
                                dt
                            );
                        }
                }
                /// PROBE LINES IF NECESSARY
                if(!grid.input_parser.lines_to_be_probed.empty()){
                    for(size_t curr_pt = 0 ; curr_pt < grid.input_parser.lines_to_be_probed.size();
                                curr_pt ++)
                        {
                            std::string which_form_to_probe = "line";
                            probe_a_field(
                                //GridCreator_NEW &grid
                                grid,
                                //std::string &which_field
                                grid.input_parser.lines_to_be_probed[curr_pt].type_field,
                                //std::string &filename
                                grid.input_parser.lines_to_be_probed[curr_pt].filename,
                                //std::string &which_form_to_probe
                                which_form_to_probe,
                                //std::vector<double> &infoOnForm
                                grid.input_parser.lines_to_be_probed[curr_pt].coords,
                                //std::vector<double> &electro_deltas
                                grid.delta_Electromagn,
                                //double current_time
                                current_time,
                                //double dt
                                dt
                            );
                        }
                }

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

                if(    grid.MPI_communicator.isRootProcess() != INT_MIN 
                    /*&& currentStep == grid.input_parser.maxStepsForOneCycleOfElectro*/){
                        printf("%s[MPI %d - Electro - Update - step %zu]%s\n"
                               "\t> Current simulation time is   %.10g seconds (over %.10g) [dt = %.10g seconds].\n"
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
							   dt,
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

    /**
     * @brief Freeing memory of ABC conditions.
     * 
     * Note: The delete[] checks for NULLPTR so no need to do 
     *      if(ptr != NULL) {delete[] ptr;}
     * but simply do delete[] ptr;
     */
    delete[] Eyx0;
    delete[] Eyx1;
    delete[] Ezx0;
    delete[] Ezx1;
    delete[] Exy0;
    delete[] Exy1;
    delete[] Ezy0;
    delete[] Ezy1;
    delete[] Exz0;
    delete[] Exz1;
    delete[] Eyz0;
    delete[] Eyz1;

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
    //std::cout << "AlgoElectro_NEW_UPDATE => Time: " << delta << " s" << std::endl;
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




void AlgoElectro_NEW::abc(   GridCreator_NEW &grid, 
            double *Ex, double *Ey, double *Ez,  
            double *Eyx0, double *Ezx0, 
            double *Eyx1, double *Ezx1, 
            double *Exy0, double *Ezy0, 
            double *Exy1, double *Ezy1, 
            double *Exz0, double *Eyz0,
            double *Exz1, double *Eyz1,
            double dt
        )
{
    size_t i, j, k;


    std::vector<double> delta_Electromagn = grid.delta_Electromagn;
    size_t size_x =0 , size_y = 0, size_z = 0; 
    size_t index, index_1Plus, index_1Moins, indexTmp;
    double c, abccoef;

    printf("omp_get _num_threads() = %d ", omp_get_num_threads());


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! |
// !!!!!!!!!!!!!!!! A MODIFIER !!!!!!!!!!!!!!!!!!!!!!!!! |
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! |
   /* printf("Attention pas la bonne vitesse de la lumière (line %d)\n",__LINE__);*/
    c = 299792458;                                   //  |
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! |
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! |



    


    /*  Dans cette section, on va considérer que les dimensions dans chaque direction
        seront correctement donnee sans correction a appliquer:
        - Plus de "-1" 
        - Tous les indices ont la même expression 
        Par contre, une correction sera appliquee en prenant en compte le colonnes 
        untiles a la communication:
        - "-2" à size_x,y,z uniquement dans "IndexTmp"
        Il restera a discuter des equations BC sur les aretes d'intersections qui, 
        quoiqu'il advienne causeront des refections d'onde.
    */

    /*printf("delta t = %.15lf\n", dt);
    printf("E_y_eps = %f\n", delta_Electromagn[0]);
    printf("delta_Electromagn[1] = %f\n", delta_Electromagn[1]);
    printf("delta_Electromagn[2] = %f\n", delta_Electromagn[2]);
    printf("delta_Electromagn[0] = %f\n", delta_Electromagn[0]);*/



    /* ABC at "x0" */

    if(grid.MPI_communicator.RankNeighbour[1] == -1){
        i = 1;

        size_x = grid.size_Ey[0];
        size_y = grid.size_Ey[1];
        size_z = grid.size_Ey[2];


        for (j = 1; j < size_y - 1; j++) // Peut-etre inverser les boucles
            for (k = 1; k < size_z -1 ; k++) {
                index = i + size_x * ( j + size_y * k);
                index_1Plus = i+1 + size_x * (j + size_y *k);
                indexTmp = (j-1) * (size_z -2) + (k-1);
                // c = 1/(sqrt(grid.E_y_eps[index]*grid.H_y_mu[index]));
                abccoef = (c*dt/delta_Electromagn[1] -1) / (c*dt/delta_Electromagn[1] +1);
                Ey[index] = Eyx0[indexTmp] +
                    abccoef * (Ey[index_1Plus] - Ey[index]);
                Eyx0[indexTmp] = Ey[index_1Plus];
            }

    

        size_x = grid.size_Ez[0];
        size_y = grid.size_Ez[1];
        size_z = grid.size_Ez[2];

        for (j = 1; j < size_y - 1 ; j++)
            for (k = 1; k < size_z - 1; k++) {
                index = i + size_x * ( j + size_y * k);
                index_1Plus = i+1 + size_x * (j + size_y *k);
                indexTmp = (j-1) * (size_z -2) + (k-1);
                // c = 1/(sqrt(grid.E_z_eps[index]*grid.H_z_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ez[index] = Ezx0[indexTmp] +
                    abccoef * (Ez[index_1Plus] - Ez[index]);
                Ezx0[indexTmp] = Ez[index_1Plus];

            }

    } // End if

/* ABC at "x1" */   

    if(grid.MPI_communicator.RankNeighbour[0] == -1){
        
        size_x = grid.size_Ey[0];
        size_y = grid.size_Ey[1];
        size_z = grid.size_Ey[2];

        i = size_x - 2;

        for (j = 1; j < size_y - 1; j++) // -1 not to take the last column which is to send
            for (k = 1; k < size_z - 1; k++) { // -1 not to take the last column
                index = i + size_x * ( j + size_y * k); //
                index_1Moins = i-1 + size_x * (j + size_y *k);
                indexTmp = (j-1) * (size_z - 2) + (k-1) ;
                // c = 1.0/(sqrt(grid.E_y_eps[index]*grid.H_y_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ey[index] = Eyx1[indexTmp] +
                    abccoef * (Ey[index_1Moins] - Ey[index]);
                Eyx1[indexTmp] = Ey[index_1Moins];
            }

        

        size_x = grid.size_Ez[0];
        size_y = grid.size_Ez[1];
        size_z = grid.size_Ez[2];
        
        i = size_x - 2;

        for ( j = 1; j < size_y - 1; j++) 
            for (k = 1; k < size_z -  1; k++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Moins = i-1 + size_x * (j + size_y *k);
                indexTmp = (j-1) * (size_z - 2) + (k-1) ;
                // c = 1.0/(sqrt(grid.E_z_eps[index]*grid.H_z_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ez[index] = Ezx1[indexTmp] +
                    abccoef * (Ez[index_1Moins] - Ez[index]);
                Ezx1[indexTmp] = Ez[index_1Moins];
            }

    } // End if

    /* ABC at "y0" */

    if(grid.MPI_communicator.RankNeighbour[2] == -1){
        j = 1;

        size_x = grid.size_Ex[0];
        size_y = grid.size_Ex[1];
        size_z = grid.size_Ex[2];

        for (i = 1; i < size_x - 1; i++)
            for (k = 1; k < size_z - 1; k++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Plus = i + size_x * ((j+1) + size_y *k);
                indexTmp = (i-1) * (size_z - 2) + (k-1) ;
                // c = 1.0/(sqrt(grid.E_x_eps[index]*grid.H_x_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ex[index] = Exy0[indexTmp] +
                    abccoef * (Ex[index_1Plus] - Ex[index]);
                Exy0[indexTmp] = Ex[index_1Plus];
 
            }
        
        size_x = grid.size_Ez[0];
        size_y = grid.size_Ez[1];
        size_z = grid.size_Ez[2];

        for (i = 1; i < size_x - 1; i++)
            for (k = 1; k < size_z - 1; k++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Plus = i + size_x * ((j+1) + size_y *k);
                indexTmp = (i-1) * (size_z - 2) + (k-1) ;
                // c = 1.0/(sqrt(grid.E_z_eps[index]*grid.H_z_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ez[index] = Ezy0[indexTmp] +
                    abccoef * (Ez[index_1Plus] - Ez[index]);
                Ezy0[indexTmp] = Ez[index_1Plus];
            }
    } //End if


    
    /* ABC at "y1" */

    if(grid.MPI_communicator.RankNeighbour[3] == -1){

        size_x = grid.size_Ex[0];
        size_y = grid.size_Ex[1];
        size_z = grid.size_Ex[2];

        j = size_y - 2;

        for (i = 1; i < size_x - 1; i++)
            for (k = 1; k < size_z -1; k++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Moins = i + size_x * ((j-1) + size_y *k);
                indexTmp = (i-1) * (size_z - 2) + (k-1) ;
                // c = 1.0/(sqrt(grid.E_x_eps[index]*grid.H_x_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ex[index] = Exy1[indexTmp] +
                    abccoef * (Ex[index_1Moins] - Ex[index]);
                Exy1[indexTmp] = Ex[index_1Moins];

            }

        size_x = grid.size_Ez[0];
        size_y = grid.size_Ez[1];
        size_z = grid.size_Ez[2];

        j = size_y - 2;

        for (i = 1; i < size_x - 1; i++)
            for (k = 1; k < size_z - 1; k++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Moins = i + size_x * ((j-1) + size_y *k);
                indexTmp = (i-1) * (size_z - 2) + (k-1) ;
                // c = 1.0/(sqrt(grid.E_z_eps[index]*grid.H_z_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ez[index] = Ezy1[indexTmp] +
                abccoef * (Ez[index_1Moins] - Ez[index]);
                Ezy1[indexTmp] = Ez[index_1Moins];

            }
    } // End if



    /* ABC at "z0" (bottom) */

    if(grid.MPI_communicator.RankNeighbour[4] == -1){
        k = 1;


        size_x = grid.size_Ex[0];
        size_y = grid.size_Ex[1];
        size_z = grid.size_Ex[2];

        for (i = 1; i < size_x - 1; i++)
            for (j = 1; j < size_y -1 ; j++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Plus = i + size_x * (j + size_y *(k+1));
                indexTmp = (i-1) * (size_y - 2) + (j-1) ;
                // c = 1.0/(sqrt(grid.E_x_eps[index]*grid.H_x_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ex[index] = Exz0[indexTmp] +
                abccoef * (Ex[index_1Plus] - Ex[index]);
                Exz0[indexTmp] = Ex[index_1Plus];
 
            }

        size_x = grid.size_Ey[0];
        size_y = grid.size_Ey[1];
        size_z = grid.size_Ey[2];

        for (i = 1; i < size_x - 1; i++)
            for (j = 1; j < size_y - 1; j++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Plus = i + size_x * (j + size_y *(k+1));
                indexTmp = (i-1) * (size_y - 2) + (j-1) ;
                // c = 1.0/(sqrt(grid.E_y_eps[index]*grid.H_y_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ey[index] = Eyz0[indexTmp] +
                abccoef * (Ey[index_1Plus] - Ey[index]);
                Eyz0[indexTmp] = Ey[index_1Plus];

            }
    } //End if



    /* ABC at "z1" (top) */

    if(grid.MPI_communicator.RankNeighbour[5] == -1){

        size_x = grid.size_Ex[0];
        size_y = grid.size_Ex[1];
        size_z = grid.size_Ex[2];

        k = size_z - 2;

        for (i = 1; i < size_x - 1; i++)
            for (j = 1; j < size_y - 1; j++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Moins = i + size_x * (j + size_y *(k-1));
                indexTmp = (i-1) * (size_y - 2) + (j-1) ;
                // c = 1.0/(sqrt(grid.E_x_eps[index]*grid.H_x_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ex[index] = Exz1[indexTmp] +
                abccoef * (Ex[index_1Moins] - Ex[index]);
                Exz1[indexTmp] = Ex[index_1Moins];

            }

        size_x = grid.size_Ey[0];
        size_y = grid.size_Ey[1];
        size_z = grid.size_Ey[2];

        k = size_z - 2;

        for (i = 1; i < size_x - 1; i++)
            for (j = 1; j < size_y - 1; j++) {
                index = i + size_x * ( j + size_y * k); //
                index_1Moins = i + size_x * (j + size_y *(k-1));
                indexTmp = (i-1) * (size_y - 2) + (j-1) ;
                // c = 1.0/(sqrt(grid.E_y_eps[index]*grid.H_y_mu[index]));
                abccoef = (c*dt/delta_Electromagn[2] -1) / (c*dt/delta_Electromagn[2] +1);
                Ey[index] = Eyz1[indexTmp] +
                abccoef * (Ey[index_1Moins] - Ey[index]);
                Eyz1[indexTmp] = Ey[index_1Moins];

            }

    } // End if 


    return;
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

/**
 * @brief Function to probe the value of a field in time, and write it inside a file.
 */
void probe_a_field(
    GridCreator_NEW &grid,
    std::string &which_field,
    std::string &filename,
    std::string &which_form_to_probe,
    std::vector<double> &infoOnForm,
    std::vector<double> &electro_deltas,
    double current_time,
    double dt
)
{
    if(which_form_to_probe == "point"){
        /**
         * @brief Probing the value of the field at a given point.
         * The point is given in (x,y,z) coordinates or in global node number.
         */
        if(infoOnForm.size() != 3){
            fprintf(stderr,"In %s :: ERROR :: You asked for probing a point"
                            " but the 'infoOnForm' vector is not of size 3"
                            " (has %zu). Aborting.\n",
                            __FUNCTION__,infoOnForm.size());
            fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
            MPI_Abort(MPI_COMM_WORLD,-1);
        }

        /// Check if the point is given as (x,y,z) or node number:
        if(   round(infoOnForm[0]) == infoOnForm[0]
           && round(infoOnForm[1]) == infoOnForm[1]
           && round(infoOnForm[2]) == infoOnForm[2])
        {
            /// The info on the point is given as node number.
            DISPLAY_ERROR_ABORT("Not implemented.");
        }else{
            /// The info on the point is given as coordinates.
            /// Compute global node number and check if it is inside the current MPI process:
            size_t nbr_X_gl = infoOnForm[0] / electro_deltas[0];
            size_t nbr_Y_gl = infoOnForm[1] / electro_deltas[1];
            size_t nbr_Z_gl = infoOnForm[2] / electro_deltas[2];
            if(grid.is_global_inside_me(nbr_X_gl,nbr_Y_gl,nbr_Z_gl) == true){
                /// The node is inside this MPI process. Proceed.
                /// Open the file:
                int fd = -1;
                while (fd == -1) fd = tryGetLock(filename.c_str());

                FILE *file = NULL;
                if(NULL == (file = fopen(filename.c_str(),"a"))){
                    DISPLAY_ERROR_ABORT(
                        "Cannot open the file %s.",filename.c_str()
                    );
                }else{
                    /// Get the local node numbering:
                    size_t nbr_X_loc = 0;
                    size_t nbr_Y_loc = 0;
                    size_t nbr_Z_loc = 0;
                    bool is_ok = false;
                    grid.get_local_from_global_electro(
                        nbr_X_gl  ,nbr_Y_gl  ,nbr_Z_gl,
                        &nbr_X_loc,&nbr_Y_loc,&nbr_Z_loc,
                        &is_ok
                    );
                    if(!is_ok){
                        DISPLAY_ERROR_ABORT(
                            "There was an error inside get_local_from_global."
                        );
                    }
                    std::string size = "size_";
                    size.append(which_field);
                    std::vector<size_t> sizes = grid.get_fields_size(size);
                    size_t index = nbr_X_loc + sizes[0] * (nbr_Y_loc + sizes[1] * nbr_Z_loc);

                    if(index >= sizes[0]*sizes[1]*sizes[2]){
                        DISPLAY_ERROR_ABORT(
                            "index is out of bounds. Field is %s (size %zu) and index is %zu.",
                            which_field.c_str(),sizes[0]*sizes[1]*sizes[2],index
                        );
                    }

                    double value = grid.get_fields(which_field)[index];

                    fprintf(file,"(%.10g,%.10g,%.10g,%.10g) %s = %.10g [gl_node(%zu,%zu,%zu)| dt %.10g | step %zu]\n",
                                current_time,
                                infoOnForm[0],
                                infoOnForm[1],
                                infoOnForm[2],
                                which_field.c_str(),
                                value,
                                nbr_X_gl,nbr_Y_gl,nbr_Y_gl,
                                dt,
                                (size_t)(current_time/dt));
                    fclose(file);
                }
                releaseLock(fd);
            }else{
                /// The node is not inside this MPI process. Return.
                return;
            }
        }
	}else if(which_form_to_probe == "line"){
		/**
		 * @brief Export the value of a field over a specified line.
		 *	The informations must be given in terms of coordinates.
		 */
        if(infoOnForm.size() != 3){
            DISPLAY_ERROR_ABORT("You didn't provide 3 infos for a lines.");
        }
        
        /// In which case are we?
        unsigned int direction = 10;
        double longueur = -1;
        size_t global_nbr[3];
        size_t local[3];
        
        if(       std::isnan(infoOnForm[0])){
            direction = 0;
            longueur = grid.input_parser.lengthX_WholeDomain_Electro;
            global_nbr[1] = infoOnForm[1] / electro_deltas[1];
            global_nbr[2] = infoOnForm[2] / electro_deltas[2];
            printf("Case ALL along x. Direction %u. Length %lf. global_nbr{%zu,%zu,%zu}\n",
                direction,
                longueur,
                global_nbr[0],global_nbr[1],global_nbr[2]);
        }else if( std::isnan(infoOnForm[1])){
            direction = 1;
            longueur = grid.input_parser.lengthY_WholeDomain_Electro;
            global_nbr[0] = infoOnForm[0] / electro_deltas[0];
            global_nbr[2] = infoOnForm[2] / electro_deltas[2];
            printf("Case ALL along y. Direction %u. Length %lf. global_nbr{%zu,%zu,%zu}\n",
                direction,
                longueur,
                global_nbr[0],global_nbr[1],global_nbr[2]);
        }else if( std::isnan(infoOnForm[2])){
            direction = 2;
            longueur = grid.input_parser.lengthZ_WholeDomain_Electro;
            global_nbr[0] = infoOnForm[0] / electro_deltas[0];
            global_nbr[1] = infoOnForm[1] / electro_deltas[1];
            printf("Case ALL along z. Direction %u. Length %lf. global_nbr{%zu,%zu,%zu}\n",
                direction,
                longueur,
                global_nbr[0],global_nbr[1],global_nbr[2]);
        }else{
            printf("infoOnForm[%lf,%lf,%lf].\n",
                infoOnForm[0],
                infoOnForm[1],
                infoOnForm[2]);
            DISPLAY_ERROR_ABORT("I don't understand what you want.");
        }
        
        /// Retrieve size of field to probe:
        std::string size = "size_";
        size.append(which_field);
        std::vector<size_t> sizes = grid.get_fields_size(size);
        
        std::vector<double> data;
        std::vector<size_t> NBR;
        
        /// Loop in the direction needed:
        for(size_t I = 0 ; I < longueur/electro_deltas[direction] +1 ; I ++){
            global_nbr[direction] = I;
            /// Check if the global node is inside the MPI process:
            if(grid.is_global_inside_me(global_nbr[0],global_nbr[1],global_nbr[2]) == true){
                /// Get local node number:
                bool is_ok = false;
                grid.get_local_from_global_electro(
                    global_nbr[0],global_nbr[1],global_nbr[2],
                    &local[0],&local[1],&local[2],
                    &is_ok
                );
                if(!is_ok){
                    DISPLAY_ERROR_ABORT(
                        "There was an error inside get_local_from_global."
                    );
                }
                size_t index = local[0] + sizes[0] * (local[1] + sizes[1] * local[2]);
                data.push_back(grid.get_fields(which_field)[index]);
                NBR.push_back(I);
            }
        }
        
        /*for(size_t I = 0 ; I < data.size() ; I ++){
            printf("DATA[%zu] = %lf.\n",I,data[I]);
        }*/
        
        /// MPI safe file writing:
        int TAG = 15;
        if(grid.MPI_communicator.getRank() == grid.MPI_communicator.rootProcess){
            /// The root process starts
            /// Open the file:
            int fd = -1;
            while (fd == -1) fd = tryGetLock(filename.c_str());

            FILE *file = NULL;
            if(NULL == (file = fopen(filename.c_str(),"a"))){
                DISPLAY_ERROR_ABORT(
                    "Cannot open the file %s.",filename.c_str()
                );
            }
            /// Write:
            for(size_t K = 0 ; K < data.size() ; K ++){
                if(       direction == 0){
                    global_nbr[0]    = NBR[K];
                }else if( direction == 1){
                    global_nbr[1]    = NBR[K];
                }else if( direction == 2){
                    global_nbr[2]    = NBR[K];
                }
                fprintf(file,"%s at (%.10g,%.10g,%.10g) is %.10g at t = %.10g | step %zu | dt = %.10g.\n",
                    which_field.c_str(),
                    global_nbr[0]*electro_deltas[0],
                    global_nbr[1]*electro_deltas[1],
                    global_nbr[2]*electro_deltas[2],
                    data[K],
                    current_time,
                    (size_t)(current_time/dt),
                    dt);
            }
            /// Close file:
            fclose(file);
            releaseLock(fd);
            
            int OK = 1;
            if(grid.MPI_communicator.getNumberOfMPIProcesses() >= 2)
                MPI_Send(&OK,1,MPI_INT,
                         grid.MPI_communicator.getRank()+1,
                         TAG,MPI_COMM_WORLD);
            
        }else{
            /// Wait for previous MPI to finish:
            int OK = 1;
            MPI_Recv(&OK,1,MPI_INT,
                     grid.MPI_communicator.getRank()-1,
                     TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(grid.MPI_communicator.getRank() != grid.MPI_communicator.getNumberOfMPIProcesses()-1){
                MPI_Send(&OK,1,MPI_INT,
                     grid.MPI_communicator.getRank()+1,
                     TAG,MPI_COMM_WORLD);
            }
        }
            
        
		//DISPLAY_ERROR_ABORT("Not yet implemented.");
    }else{
		DISPLAY_ERROR_ABORT(
			"There is no option corresponding to %s.",
			which_form_to_probe.c_str()
		);
	}
}

/*! Try to get lock. Return its file descriptor or -1 if failed.
 *
 *  @param lockName Name of file used as lock (i.e. '/var/lock/myLock').
 *  @return File descriptor of lock file, or -1 if failed.
 */
int tryGetLock( char const *lockName )
{
    mode_t m = umask( 0 );
    int fd = open( lockName, O_RDWR|O_CREAT, 0666 );
    umask( m );
    if( fd >= 0 && flock( fd, LOCK_EX | LOCK_NB ) < 0 )
    {
        close( fd );
        fd = -1;
    }
    return fd;
}

/*! Release the lock obtained with tryGetLock( lockName ).
 *
 *  @param fd File descriptor of lock returned by tryGetLock( lockName ).
 */
void releaseLock( int fd)
{
    if( fd < 0 )
        return;
    close( fd );
}

inline double absolute_value_sinus_func(double x){
	return fabs(std::sin(x));
}

/**
 * @brief Steady-state analyser.
 */
bool AlgoElectro_NEW::SteadyStateAnalyser(void){

    /// Compute the integral of the absolute value of a sine wave:
    double res = GaussLobattoInt(&absolute_value_sinus_func,
					0, 10,
					1e-10,
					1E5);
    DISPLAY_ERROR_ABORT(
        "Not implemented yet."
    );
    return false;
}
