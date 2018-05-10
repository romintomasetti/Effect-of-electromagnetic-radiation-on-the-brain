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
    if(    mesh.input_parser.conditionsInsideSources[0] == "FACE_EX_Electric_along_Z"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_EX_Electric_along_Y"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_EY_Electric_along_Z"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_EY_Electric_along_X"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_EZ_Electric_along_Y"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_EZ_Electric_along_X"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_Minus_EX_Electric_along_Z"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_Minus_EX_Electric_along_Y"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_Minus_EY_Electric_along_Z"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_Minus_EY_Electric_along_X"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_Minus_EZ_Electric_along_X"
        || mesh.input_parser.conditionsInsideSources[0] == "FACE_Minus_EZ_Electric_along_Y")
        {
            double c = 1/(sqrt(VACUUM_PERMEABILITY*VACUUM_PERMITTIVITY));
            dt = dx / c;
            DISPLAY_WARNING(
                "Using the 1D case formula for computing dt !"
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
    if(grid.input_parser.apply_ABC_BCs == true){
        Eyx0 = new double[size]();   
        Eyx1 = new double[size]();
    }
    //std::fill_n(Eyx0, size, 0);
    //std::fill_n(Eyx1, size, 0);

    // std::vector<double> Eyx0(size, 0.0);
    // std::vector<double> Eyx1(size, 0.0);

    // ABC Old Tangential Field Ez at the extrmities of x of the grid:
    size = (grid.size_Ez[1]-2)*(grid.size_Ez[2]-2);
    double *Ezx0    = NULL;
    double *Ezx1    = NULL;
    if(grid.input_parser.apply_ABC_BCs == true){
        Ezx0 = new double[size]();  
        Ezx1 = new double[size]();
    }
    //std::fill_n(Ezx0, size, 0);
    //std::fill_n(Ezx1, size, 0);

    // std::vector<double> Ezx0(size, 0.0);
    // std::vector<double> Ezx1(size, 0.0);

    // ABC Old Tangential Field Ex at the extrmities of y of the grid:
    size = (grid.size_Ex[0]-2)*(grid.size_Ex[2]-2);
    double *Exy0    = NULL;
    double *Exy1    = NULL;
    if(grid.input_parser.apply_ABC_BCs == true){
        Exy0 = new double[size]();    
        Exy1 = new double[size]();
    }
    //std::fill_n(Exy0, size, 0);
    //std::fill_n(Exy1, size, 0);

    // std::vector<double> Exy0(size, 0.0);
    // std::vector<double> Exy1(size, 0.0);

    // ABC Old Tangential Field Ez at the extrmities of y of the grid:
    size = (grid.size_Ez[0]-2)*(grid.size_Ez[2]-2);
    double *Ezy0    = NULL;    
    double *Ezy1    = NULL;
    if(grid.input_parser.apply_ABC_BCs == true){
        Ezy0 = new double[size]();
        Ezy1 = new double[size]();
    }
    //std::fill_n(Ezy0, size, 0);
    //std::fill_n(Ezy1, size, 0);
    // std::vector<double> Ezy0(size, 0.0);
    // std::vector<double> Ezy1(size, 0.0);

    // ABC Old Tangential Field Ex at the extrmities of z of the grid:
    size = (grid.size_Ex[0]-2)*(grid.size_Ex[1]-2);
    double *Exz0    = NULL;
    double *Exz1    = NULL;
    if(grid.input_parser.apply_ABC_BCs == true){
        Exz0 = new double[size]();
        Exz1 = new double[size]();
    }
    //std::fill_n(Exz0, size, 0);
    //std::fill_n(Exz1, size, 0);
    
    // std::vector<double> Exz0(size, 0.0);
    // std::vector<double> Exz1(size, 0.0);

    // ABC Old Tangential Field Ey at the extrmities of z of the grid:
    size = (grid.size_Ey[0]-2)*(grid.size_Ey[1]-2);
    double *Eyz0    = NULL;
    double *Eyz1    = NULL;
    if(grid.input_parser.apply_ABC_BCs == true){
        Eyz0 = new double[size]();    
        Eyz1 = new double[size]();    
    }
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
    bool DO_3D_SOURCE = true;
    if(    grid.input_parser.conditionsInsideSources[0] == "FACE_EX_Electric_along_Z"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_EX_Electric_along_Y"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_EY_Electric_along_Z"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_EY_Electric_along_X"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_EZ_Electric_along_Y"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_EZ_Electric_along_X"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EX_Electric_along_Z"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EX_Electric_along_Y"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EY_Electric_along_Z"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EY_Electric_along_X"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EZ_Electric_along_X"
        || grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EZ_Electric_along_Y")
        {
            DO_3D_SOURCE = false;
        }
    
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
	
    if(DO_3D_SOURCE == true){
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
    }
	
	
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

    /* For the Ez field */
    size_t nbr_points_Ez = grid.size_Ez[0] * grid.size_Ez[1] * grid.size_Ez[2];
    std::vector<double> Ez_trapz_absolute(nbr_points_Ez,0);

    /* For the Ey field */
    size_t nbr_points_Ey = grid.size_Ey[0] * grid.size_Ey[1] * grid.size_Ey[2];
    std::vector<double> Ey_trapz_absolute(nbr_points_Ey,0);

    /* For the Ey field */
    size_t nbr_points_Ex = grid.size_Ex[0] * grid.size_Ex[1] * grid.size_Ex[2];
    std::vector<double> Ex_trapz_absolute(nbr_points_Ex,0);

    /* Look over a given period of time, derived from the smallest source frequency **/
    double min_frequency = DBL_MAX;
    for(size_t I = 0 ; I < grid.input_parser.source.frequency.size() ; I ++){
        if(grid.input_parser.source.frequency[I] < min_frequency)
            min_frequency = grid.input_parser.source.frequency[I];
    }
        // Look over 35 periods of the signal:
    double look_over_period_SS = 100 * 1 / min_frequency;
    size_t look_over_steps__SS = std::floor(look_over_period_SS/dt);
    size_t counter_step_____SS = 0;
        // Number of nodes at steady state:
    size_t number_of_points_EZ_at_speady_state = 0;
        // Boolean value:
    bool is_steady_state_for_all_mpi  = false;
    bool is_steady_state_for_this_MPI = false;
    size_t counter_steady_state_segments = 0;
        // Begin and end time of the time period:
    double time_beg = 0.0, time_end = 0.0;
        // Wait that the wave has time to travel through the whole domain at least twice
        // before checking for steadiness:
    double max_length = 0.0;
    double LX = grid.input_parser.lengthX_WholeDomain_Electro;
    double LY = grid.input_parser.lengthY_WholeDomain_Electro;
    double LZ = grid.input_parser.lengthZ_WholeDomain_Electro;

    max_length = LX;
    if(LX < LY || LX < LZ){
        if(LY < LZ)
            max_length = LZ;
        else
            max_length = LY;
    }
    double min_time_before_checking_steadiness
        = 2* max_length / (1 / std::sqrt( VACUUM_PERMEABILITY * VACUUM_PERMITTIVITY));
        // at least 2 segments must be steady:
    bool is_previous_segment_steady = false;
    size_t numero_du_segment_precedent = 1E9;

    ///////////////////////////////////////////
    /// END VARIABLES FOR STEADY-STATE CHECKING ///
    ///////////////////////////////////////////

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
        firstprivate(Eyz0, Eyz1)\
        firstprivate(Ez_trapz_absolute)\
        firstprivate(Ey_trapz_absolute)\
        firstprivate(Ex_trapz_absolute)\
        shared(number_of_points_EZ_at_speady_state)\
        shared(counter_step_____SS,look_over_steps__SS)\
        shared(is_steady_state_for_all_mpi)\
        shared(is_steady_state_for_this_MPI)\
        shared(counter_steady_state_segments)\
        shared(time_beg,time_end)\
        shared(min_time_before_checking_steadiness)\
        firstprivate(is_previous_segment_steady)\
        firstprivate(numero_du_segment_precedent)
    {
        /**
         * @brief Determine if the simulation must be 1D.
         * 
         * There are 6 types of 1D simulations:
         *      1) Face with normal minus e_y is imposed. Propagation in Y direction.
         *          Parameter: FACE_Minus_EY.
         *      2) Face with normal plus e_y is imposed. Propagation in -Y direction.
         *          Parameter: FACE_EY.
         *      3) Face with normal minus e_x is imposed. Propagation in X direction.
         *          Parameter: FACE_Minus_EX.
         *      4) Face with normal plus e_x is imposed. Propagation in -X direction.
         *          Parameter: FACE_EX.
         *      5) Face with normal minus e_z is imposed. Propagation in Z direction.
         *          Parameter: FACE_Minus_EZ.
         *      6) Face with normal plus e_z is imposed. Propagation in -Z direction.
         *          Parameter: FACE_EZ.
         * 
         * Note: the convention for normal is unit outward normal.
         */
        bool IS_1D_FACE_EX_Electric_along_Z       = false;
        bool IS_1D_FACE_EX_Electric_along_Y       = false;
        bool IS_1D_FACE_EY_Electric_along_Z       = false;
        bool IS_1D_FACE_EY_Electric_along_X       = false;
        bool IS_1D_FACE_EZ_Electric_along_Y       = false;
        bool IS_1D_FACE_EZ_Electric_along_X       = false;
        bool IS_1D_FACE_Minus_EX_Electric_along_Z = false;
        bool IS_1D_FACE_Minus_EX_Electric_along_Y = false;
        bool IS_1D_FACE_Minus_EY_Electric_along_Z = false;
        bool IS_1D_FACE_Minus_EY_Electric_along_X = false;
        bool IS_1D_FACE_Minus_EZ_Electric_along_X = false;
        bool IS_1D_FACE_Minus_EZ_Electric_along_Y = false;

        bool IS_1D_PROPAGATION_IN_X = false;
        bool IS_1D_PROPAGATION_IN_Y = false;
        bool IS_1D_PROPAGATION_IN_Z = false;

        bool IS_3D_CASE = true;

        if(grid.input_parser.conditionsInsideSources[0] == "FACE_EX_Electric_along_Z"){
            printf(">>> [MPI %d] - Using 1D with FACE_EX_Electric_along_Z.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_EX_Electric_along_Z = true;
            IS_3D_CASE                     = false;
            IS_1D_PROPAGATION_IN_X         = true;
        
        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_EX_Electric_along_Y"){
            printf(">>> [MPI %d] - Using 1D with FACE_EX_Electric_along_Y.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_EX_Electric_along_Y = true;
            IS_3D_CASE                     = false;
            IS_1D_PROPAGATION_IN_X         = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EX_Electric_along_Z"){
            printf(">>> [MPI %d] - Using 1D with FACE_Minus_EX_Electric_along_Z.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_Minus_EX_Electric_along_Z = true;
            IS_3D_CASE                           = false;
            IS_1D_PROPAGATION_IN_X               = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EX_Electric_along_Y"){
            printf(">>> [MPI %d] - Using 1D with FACE_Minus_EX_Electric_along_Y.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_Minus_EX_Electric_along_Y = true;
            IS_3D_CASE                           = false;
            IS_1D_PROPAGATION_IN_X               = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_EY_Electric_along_Z"){
            printf(">>> [MPI %d] - Using 1D with FACE_EY_Electric_along_Z.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_EY_Electric_along_Z = true;
            IS_3D_CASE                     = false;
            IS_1D_PROPAGATION_IN_Y         = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_EY_Electric_along_X"){
            printf(">>> [MPI %d] - Using 1D with FACE_EY_Electric_along_X.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_EY_Electric_along_X = true;
            IS_3D_CASE                     = false;
            IS_1D_PROPAGATION_IN_Y         = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EY_Electric_along_Z"){
            printf(">>> [MPI %d] - Using 1D with FACE_Minus_EY_Electric_along_Z.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_Minus_EY_Electric_along_Z = true;
            IS_3D_CASE                           = false;
            IS_1D_PROPAGATION_IN_Y               = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EY_Electric_along_X"){
            printf(">>> [MPI %d] - Using 1D with FACE_Minus_EY_Electric_along_X.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_Minus_EY_Electric_along_X = true;
            IS_3D_CASE                           = false;
            IS_1D_PROPAGATION_IN_Y               = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_EZ_Electric_along_Y"){
            printf(">>> [MPI %d] - Using 1D with FACE_EZ_Electric_along_Y.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_EZ_Electric_along_Y = true;
            IS_3D_CASE                     = false;
            IS_1D_PROPAGATION_IN_Z         = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_EZ_Electric_along_X"){
            printf(">>> [MPI %d] - Using 1D with FACE_EZ_Electric_along_X.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_EZ_Electric_along_X = true;
            IS_3D_CASE                     = false;
            IS_1D_PROPAGATION_IN_Z         = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EZ_Electric_along_Y"){
            printf(">>> [MPI %d] - Using 1D with FACE_Minus_EZ_Electric_along_ZY.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_Minus_EZ_Electric_along_Y = true;
            IS_3D_CASE                           = false;
            IS_1D_PROPAGATION_IN_Z               = true;

        }else if(grid.input_parser.conditionsInsideSources[0] == "FACE_Minus_EZ_Electric_along_X"){
            printf(">>> [MPI %d] - Using 1D with FACE_Minus_EZ_Electric_along_X.\n",
                grid.MPI_communicator.getRank());
            IS_1D_FACE_Minus_EZ_Electric_along_X = true;
            IS_3D_CASE                           = false;
            IS_1D_PROPAGATION_IN_Z               = true;
        }
        
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
            if( IS_1D_FACE_EY_Electric_along_Z
             || IS_1D_FACE_EZ_Electric_along_Y
             || IS_1D_FACE_Minus_EZ_Electric_along_Y
             || IS_1D_FACE_Minus_EY_Electric_along_Z){
                // Do nothing.
            }else{
                IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX = 1;
            }
        }

        if(grid.MPI_communicator.MPI_POSITION[1] == grid.MPI_communicator.MPI_MAX_POSI[1]){
            if(   IS_1D_FACE_EX_Electric_along_Z 
               || IS_1D_FACE_Minus_EX_Electric_along_Z
               || IS_1D_FACE_EZ_Electric_along_X
               || IS_1D_FACE_Minus_EZ_Electric_along_X){
                // Do nothing.
            }else{
                IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY = 1;
            }
        }

        if(grid.MPI_communicator.MPI_POSITION[2] == grid.MPI_communicator.MPI_MAX_POSI[2]){
            if( IS_1D_FACE_EX_Electric_along_Y
             || IS_1D_FACE_Minus_EX_Electric_along_Y
             || IS_1D_FACE_EY_Electric_along_X
             || IS_1D_FACE_Minus_EY_Electric_along_X){
                // Do nothing.
            }else{
                IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ = 1;
            }
        }

        if(grid.MPI_communicator.MPI_POSITION[0] == 0){
            if( IS_1D_FACE_EY_Electric_along_Z
             || IS_1D_FACE_EZ_Electric_along_Y
             || IS_1D_FACE_Minus_EZ_Electric_along_Y
             || IS_1D_FACE_Minus_EY_Electric_along_Z){
                // Do nothing.
            }else{
                IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX = 1;
            }
        }

        if(grid.MPI_communicator.MPI_POSITION[1] == 0){
            if(   IS_1D_FACE_EX_Electric_along_Z
               || IS_1D_FACE_Minus_EX_Electric_along_Z
               || IS_1D_FACE_EZ_Electric_along_X
               || IS_1D_FACE_Minus_EZ_Electric_along_X){
                // Do nothing.
            }else{
                IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY = 1;
            }
        }

        if(grid.MPI_communicator.MPI_POSITION[2] == 0){
            if( IS_1D_FACE_EX_Electric_along_Y
             || IS_1D_FACE_Minus_EX_Electric_along_Y
             || IS_1D_FACE_EY_Electric_along_X
             || IS_1D_FACE_Minus_EY_Electric_along_X){
                // Do nothing.
            }else{
                IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ = 1;
            }
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
                && currentStep < grid.input_parser.maxStepsForOneCycleOfElectro
                && !is_steady_state_for_all_mpi){
            #pragma omp barrier
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

            size_x = grid.size_Hx[0];
            size_y = grid.size_Hx[1];

            size_x_1 = grid.size_Ey[0];
            size_y_1 = grid.size_Ey[1];

            size_x_2 = grid.size_Ez[0];
            size_y_2 = grid.size_Ez[1];

            #pragma omp for //schedule(static) collapse(3) nowait
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
            #pragma omp barrier
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

            #pragma omp for //schedule(static) collapse(3) nowait
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
            #pragma omp barrier
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

            #pragma omp for //schedule(static) collapse(3) nowait
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
            #pragma omp barrier
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


            ///////////////////////////////////////
            // 1D CASE - ONLY WITH 1 MPI PROCESS //
            ///////////////////////////////////////
            #pragma omp master
            {
                if(IS_3D_CASE == false){
                    this->apply_1D_case_on_magnetic_field(
                        IS_1D_FACE_EX_Electric_along_Z,      
                        IS_1D_FACE_EX_Electric_along_Y ,     
                        IS_1D_FACE_EY_Electric_along_Z  ,     
                        IS_1D_FACE_EY_Electric_along_X   ,    
                        IS_1D_FACE_EZ_Electric_along_Y    ,   
                        IS_1D_FACE_EZ_Electric_along_X     , 
                        IS_1D_FACE_Minus_EX_Electric_along_Z, 
                        IS_1D_FACE_Minus_EX_Electric_along_Y ,
                        IS_1D_FACE_Minus_EY_Electric_along_Z ,
                        IS_1D_FACE_Minus_EY_Electric_along_X ,
                        IS_1D_FACE_Minus_EZ_Electric_along_X ,
                        IS_1D_FACE_Minus_EZ_Electric_along_Y,
                        grid,
                        H_x_tmp,
                        H_y_tmp,
                        H_z_tmp
                    );
                }
            }
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


            #pragma omp for //schedule(static) collapse(3) nowait
            for(K = 1 + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ ;
                    K < grid.size_Ex[2]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ ; 
                    K ++){
                for(J = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY ; 
                        J < grid.size_Ex[1]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY ; 
                        J ++){
                    for(I = 1  /*+ IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX*/; 
                            I < grid.size_Ex[0]-1 /*- IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX*/; 
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

                        double E_x_prev = E_x_tmp[index];

                        E_x_tmp[index] = C_exe[index] * E_x_tmp[index]
                                + C_exh_1[index] * (H_z_tmp[index_1Plus] - H_z_tmp[index_1Moins])
                                - C_exh_2[index] * (H_y_tmp[index_2Plus] - H_y_tmp[index_2Moins]);

                        /* COMPUTE ABSOLUTE TRAPEZOIDAL INTEGRATION */
                            trapz_without_dt(
                                std::abs(E_x_prev),
                                std::abs(E_x_tmp[index]),
                                &Ex_trapz_absolute[index]
                            );

                        /* END OF COMPUTE ABSOLUTE TRAPEZOIDAL INTEGRATION */

                    }
                }
            }
            #pragma omp barrier

            // Updating the electric field Ey.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Ey[0];
            size_y = grid.size_Ey[1];

            size_x_1 = grid.size_Hx[0];
            size_y_1 = grid.size_Hx[1];

            size_x_2 = grid.size_Hz[0];
            size_y_2 = grid.size_Hz[1];

            #pragma omp for //schedule(static) collapse(3) nowait
            for(K = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ ; 
                    K < grid.size_Ey[2]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ ; 
                    K ++){
                for(J = 1  /*+ IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY*/; 
                        J < grid.size_Ey[1]-1 /*- IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY*/ ; 
                        J ++){
                    for(I = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX ; 
                            I < grid.size_Ey[0]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX ; 
                            I ++){

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

                        double E_y_prev = E_y_tmp[index];

                        E_y_tmp[index] = C_eye[index] * E_y_tmp[index]
                                + C_eyh_1[index] * (H_x_tmp[index_1Plus] - H_x_tmp[index_1Moins])
                                - C_eyh_2[index] * (H_z_tmp[index_2Plus] - H_z_tmp[index_2Moins]);

                        /* COMPUTE ABSOLUTE TRAPEZOIDAL INTEGRATION */
                            trapz_without_dt(
                                std::abs(E_y_prev),
                                std::abs(E_y_tmp[index]),
                                &Ey_trapz_absolute[index]
                            );

                        /* END OF COMPUTE ABSOLUTE TRAPEZOIDAL INTEGRATION */

                    }
                }
            }
            #pragma omp barrier

            // Updating the electric field Ez.
            // Don't update neighboors ! Start at 1. Go to size-1.

            size_x = grid.size_Ez[0];
            size_y = grid.size_Ez[1];

            size_x_1 = grid.size_Hy[0];
            size_y_1 = grid.size_Hy[1];

            size_x_2 = grid.size_Hx[0];
            size_y_2 = grid.size_Hx[1];

            #pragma omp for //schedule(static) collapse(3) nowait
            for(K = 1  /*+ IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ*/;
                     K < grid.size_Ez[2]-1 /*- IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ*/ ; 
                     K ++){
                for(J = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY ; 
                        J < grid.size_Ez[1]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY ; 
                        J ++){
                    for(I = 1  + IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX ; 
                            I < grid.size_Ez[0]-1 - IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX ; 
                            I ++){

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

                        double Ez_prev = E_z_tmp[index];

                        E_z_tmp[index] = C_eze[index] * E_z_tmp[index]
                                + C_ezh_1[index] * (H_y_tmp[index_1Plus] - H_y_tmp[index_1Moins])
                                - C_ezh_2[index] * (H_x_tmp[index_2Plus] - H_x_tmp[index_2Moins]);

                        /** COMPUTE ABSOLUTE TRAPEZOIDAL INTEGRATION OF EZ **/
                            
                            double trapz_prev = Ez_trapz_absolute[index];
                            trapz_without_dt(
                                std::abs(Ez_prev),
                                std::abs(E_z_tmp[index]),
                                &Ez_trapz_absolute[index]
                            );

                            if(counter_step_____SS > 10){
                                if(std::abs( 
                                    Ez_trapz_absolute[index]/(double)counter_step_____SS
                                        - trapz_prev/(double)(counter_step_____SS-1))
                                    < 1E-3*Ez_trapz_absolute[index]/(double)counter_step_____SS)
                                    {
                                        #pragma omp critical
                                        {
                                            number_of_points_EZ_at_speady_state++;
                                        }
                                    }
                            }

                        /** END OF COMPUTE ABSOLUTE TRAPEZOIDAL INTEGRATION OF EZ **/
                    }
                }
            }
            #pragma omp barrier

            #pragma omp master
            {
                // Total size without neighboors:
                size_t tmp = (grid.size_Ez[0]-4) * (grid.size_Ez[1]-2) * (grid.size_Ez[2]-2);
                printf("\t>>> Number of Ez nodes that are in the steady state is %zu over %zu.\n",
                            number_of_points_EZ_at_speady_state,
                            tmp);
                if(     number_of_points_EZ_at_speady_state >= std::floor(0.99 * tmp)-1
                    &&  counter_step_____SS > 100){

                            printf("\t>>> [MPI %d] - Steady state is reached!\n",
                            grid.MPI_communicator.getRank());

                            if(is_previous_segment_steady == true
                                && numero_du_segment_precedent != counter_steady_state_segments){
                                is_steady_state_for_this_MPI = true;
                            }else{
                                is_previous_segment_steady  = true;
                                numero_du_segment_precedent = counter_steady_state_segments;
                            }
                    
                }else{
                    number_of_points_EZ_at_speady_state = 0;
                }

                if(counter_step_____SS >= look_over_steps__SS){
                    for(size_t I = 0 ; I < Ez_trapz_absolute.size() ; I ++)
                        Ez_trapz_absolute[I] = 0.0;
                    counter_step_____SS = 1;
                    printf("look_over_steps__SS %zu\n",look_over_steps__SS);
                    counter_steady_state_segments++;
                    time_beg = current_time;

                }else{
                    counter_step_____SS ++;
                }

                /* Communicate with other MPI's to know is everybody is in steady state: */
                time_end = current_time;
                /// Attention ! Check only after a given number of iterations !
                if( current_time > min_time_before_checking_steadiness ){
                    this->SteadyStateAnalyser(
                        is_steady_state_for_this_MPI,
                        &is_steady_state_for_all_mpi,
                        grid,
                        time_beg,
                        time_end,
                        Ez_trapz_absolute,
                        Ey_trapz_absolute,
                        Ex_trapz_absolute,
                        dt
                    );
                }
            }
            #pragma omp barrier

            ///////////////////////////////////////
            // 1D CASE - ONLY WITH 1 MPI PROCESS //
            ///////////////////////////////////////
            #pragma omp master
            {
                if(IS_3D_CASE == false){
                    this->apply_1D_case_on_electric_field(
                        IS_1D_FACE_EX_Electric_along_Z,      
                        IS_1D_FACE_EX_Electric_along_Y ,     
                        IS_1D_FACE_EY_Electric_along_Z  ,     
                        IS_1D_FACE_EY_Electric_along_X   ,    
                        IS_1D_FACE_EZ_Electric_along_Y    ,   
                        IS_1D_FACE_EZ_Electric_along_X     , 
                        IS_1D_FACE_Minus_EX_Electric_along_Z, 
                        IS_1D_FACE_Minus_EX_Electric_along_Y ,
                        IS_1D_FACE_Minus_EY_Electric_along_Z ,
                        IS_1D_FACE_Minus_EY_Electric_along_X ,
                        IS_1D_FACE_Minus_EZ_Electric_along_X ,
                        IS_1D_FACE_Minus_EZ_Electric_along_Y,
                        current_time,
                        currentStep,
                        grid,
                        E_x_tmp,
                        E_y_tmp,
                        E_z_tmp
                    );
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

            // Do sources only if 3D case.
            if(IS_3D_CASE == true ){
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
                            double period    = 1/false_freq;
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
                            double period    = 1/frequency;
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
                            double period    = 1/frequency;
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
                            double period    = 1 / frequency;
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
                            double period    = 1 / frequency;
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
                            double period    = 1 / frequency;
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
                if(grid.input_parser.apply_ABC_BCs == true){
                    this->abc(grid,
                        E_x_tmp, E_y_tmp, E_z_tmp, 
                        Eyx0, Ezx0, 
                        Eyx1, Ezx1, 
                        Exy0, Ezy0, 
                        Exy1, Ezy1, 
                        Exz0, Eyz0, 
                        Exz1, Eyz1,
                        dt,
                        IS_1D_FACE_EX_Electric_along_Z  ,    
                        IS_1D_FACE_EX_Electric_along_Y,      
                        IS_1D_FACE_EY_Electric_along_Z ,      
                        IS_1D_FACE_EY_Electric_along_X  ,     
                        IS_1D_FACE_EZ_Electric_along_Y   ,    
                        IS_1D_FACE_EZ_Electric_along_X    ,  
                        IS_1D_FACE_Minus_EX_Electric_along_Z, 
                        IS_1D_FACE_Minus_EX_Electric_along_Y ,
                        IS_1D_FACE_Minus_EY_Electric_along_Z ,
                        IS_1D_FACE_Minus_EY_Electric_along_X ,
                        IS_1D_FACE_Minus_EZ_Electric_along_X ,
                        IS_1D_FACE_Minus_EZ_Electric_along_Y
                    );
                }
				/*fflush_stdout();
				MPI_Barrier(MPI_COMM_WORLD);
				printf("[MPI %d] - Ending ABC.\n",grid.MPI_communicator.getRank());*/
            }
            #pragma omp barrier


            gettimeofday( &end___while_iter , NULL);
            total_while_iter += end___while_iter.tv_sec  - start_while_iter.tv_sec + 
                                (end___while_iter.tv_usec - start_while_iter.tv_usec) / 1.e6;

            currentStep ++;

            
            #pragma omp barrier
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
                        std::string applyABC_str;
                        if( grid.input_parser.apply_ABC_BCs == true ){
                            applyABC_str = "true";
                        }else{
                            applyABC_str = "false";
                        }
                        std::string applyPML_str;
                        if( grid.input_parser.apply_PML_BCs == true){
                            applyPML_str = "true";
                        }else{
                            applyPML_str = "false";
                        }
                        
                        printf("%s[MPI %d - Electro - Update - step %zu]%s\n"
                               "\t> Current simulation time is   %.10g seconds (over %.10g) [dt = %.10g seconds].\n"
                               "\t> Current step is              %zu over %zu.\n"
                               "\t> Time elapsed inside while is %.6lf seconds (TOTAL).\n"
                               "\t> Time elaspsed per iter. is   %.6lf seconds (ITER).\n"
                               "\t> Time spent in MPI comm. is   %.6lf seconds (TOTAL).\n"
                               "\t> Time spent in MPI comm. is   %.6lf seconds (ITER).\n"
                               "\t> Time spent in writing is     %.6lf seconds (TOTAL).\n"
                               "\t> Using %d MPI process(es) and %d OMP thread(s).\n"
                               "\t> Apply ABC boundary conditions is set to %s.\n"
                               "\t> Apply PML boundary conditions is set to %s.\n\n",
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
                               omp_get_num_threads(),
                               applyABC_str.c_str(),
                               applyPML_str.c_str()
                               );
                    }

                if(!IS_3D_CASE){
                    if(IS_1D_PROPAGATION_IN_X){
                        printf("\t> Your 1D wave propagates in the X direction.\n");
                    }else if(IS_1D_PROPAGATION_IN_Y){
                        printf("\t> Your 1D wave propagates in the Y direction.\n");
                    }else if(IS_1D_PROPAGATION_IN_Z){
                        printf("\t> Your 1D wave propagates in the Z direction.\n");
                    }else{
                        printf(
                            "The simulation is not 3D but there is no propagation f a 1D wave."
                        );
                        std::abort();
                    }
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




void AlgoElectro_NEW::abc(   GridCreator_NEW &grid, 
            double *Ex, double *Ey, double *Ez,  
            double *Eyx0, double *Ezx0, 
            double *Eyx1, double *Ezx1, 
            double *Exy0, double *Ezy0, 
            double *Exy1, double *Ezy1, 
            double *Exz0, double *Eyz0,
            double *Exz1, double *Eyz1,
            double dt,
            bool IS_1D_FACE_EX_Electric_along_Z  ,
            bool IS_1D_FACE_EX_Electric_along_Y ,
            bool IS_1D_FACE_EY_Electric_along_Z ,
            bool IS_1D_FACE_EY_Electric_along_X  ,
            bool IS_1D_FACE_EZ_Electric_along_Y , 
            bool IS_1D_FACE_EZ_Electric_along_X  ,
            bool IS_1D_FACE_Minus_EX_Electric_along_Z ,
            bool IS_1D_FACE_Minus_EX_Electric_along_Y,
            bool IS_1D_FACE_Minus_EY_Electric_along_Z ,
            bool IS_1D_FACE_Minus_EY_Electric_along_X ,
            bool IS_1D_FACE_Minus_EZ_Electric_along_X ,
            bool IS_1D_FACE_Minus_EZ_Electric_along_Y
        )
{
    size_t i, j, k;


    std::vector<double> delta_Electromagn = grid.delta_Electromagn;
    size_t size_x =0 , size_y = 0, size_z = 0; 
    size_t index, index_1Plus, index_1Moins, indexTmp;
    double c,abccoef;

    //printf("omp_get _num_threads() = %d ", omp_get_num_threads());


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! |
// !!!!!!!!!!!!!!!! A MODIFIER !!!!!!!!!!!!!!!!!!!!!!!!! |
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! |
   /* printf("Attention pas la bonne vitesse de la lumière (line %d)\n",__LINE__);*/

    c = 1/sqrt(VACUUM_PERMEABILITY*VACUUM_PERMITTIVITY);                                //  |
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! |
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! |

    bool apply_x0_ABC = true;
    bool apply_x1_ABC = true;
    bool apply_y0_ABC = true;
    bool apply_y1_ABC = true;
    bool apply_z0_ABC = true;
    bool apply_z1_ABC = true;

    // Convention: unit OUTWARD normal.

    if( IS_1D_FACE_Minus_EY_Electric_along_Z || IS_1D_FACE_Minus_EY_Electric_along_X){
        // Propagation in +Y direction. Apply ABC on y1.
        apply_x0_ABC = false;
        apply_x1_ABC = false;
        apply_y0_ABC = false;
        apply_y1_ABC = true;
        apply_z0_ABC = false;
        apply_z1_ABC = false;
        printf("\t> In ABC : applying only y1.\n");

    }else if( IS_1D_FACE_EY_Electric_along_Z || IS_1D_FACE_EY_Electric_along_X){
        // Propagation in -Y direction. Apply ABC on y0.
        apply_x0_ABC = false;
        apply_x1_ABC = false;
        apply_y0_ABC = true;
        apply_y1_ABC = false;
        apply_z0_ABC = false;
        apply_z1_ABC = false;
        printf("\t> In ABC : applying only y0.\n");
    
    }else if( IS_1D_FACE_Minus_EX_Electric_along_Z || IS_1D_FACE_Minus_EX_Electric_along_Y){
        // Propagation in +X direction. Apply ABC on x1.
        apply_x0_ABC = false;
        apply_x1_ABC = true;
        apply_y0_ABC = false;
        apply_y1_ABC = false;
        apply_z0_ABC = false;
        apply_z1_ABC = false;
        printf("\t> In ABC : applying only x1.\n");

    }else if( IS_1D_FACE_EX_Electric_along_Z || IS_1D_FACE_EX_Electric_along_Y){
        // Propagation in -X direction. Apply ABC on x0.
        apply_x0_ABC = true;
        apply_x1_ABC = false;
        apply_y0_ABC = false;
        apply_y1_ABC = false;
        apply_z0_ABC = false;
        apply_z1_ABC = false;
        printf("\t> In ABC : applying only x0.\n");

    }else if( IS_1D_FACE_Minus_EZ_Electric_along_Y || IS_1D_FACE_Minus_EZ_Electric_along_X){
        // Propagation in +Z direction. Apply ABC on z1.
        apply_x0_ABC = false;
        apply_x1_ABC = false;
        apply_y0_ABC = false;
        apply_y1_ABC = false;
        apply_z0_ABC = false;
        apply_z1_ABC = true;
        printf("\t> In ABC : applying only z1.\n");

    }else if( IS_1D_FACE_EZ_Electric_along_Y || IS_1D_FACE_EZ_Electric_along_X){
        // Propagation in -Z direction. Apply ABC on z0.
        apply_x0_ABC = false;
        apply_x1_ABC = false;
        apply_y0_ABC = false;
        apply_y1_ABC = false;
        apply_z0_ABC = true;
        apply_z1_ABC = false;
        printf("\t> In ABC : applying only z0.\n");
    }



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

    if(grid.MPI_communicator.RankNeighbour[1] == -1 && apply_x0_ABC){
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

    if(grid.MPI_communicator.RankNeighbour[0] == -1 && apply_x1_ABC){
        
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

    if(grid.MPI_communicator.RankNeighbour[2] == -1 && apply_y0_ABC){
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

    if(grid.MPI_communicator.RankNeighbour[3] == -1 && apply_y1_ABC){

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

    if(grid.MPI_communicator.RankNeighbour[4] == -1 && apply_z0_ABC){
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

    if(grid.MPI_communicator.RankNeighbour[5] == -1 && apply_z1_ABC){

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
	return std::abs(std::sin(x));
}

/**
 * @brief This function communicates with other MPI's to know is every MPI is in steady state.
 *      Than, the amplitude of the electric field at each node is computed.
 */
bool AlgoElectro_NEW::SteadyStateAnalyser(
    const bool is_steady_state_for_this_MPI,
    bool *is_steady_state_for_all_mpi,
    GridCreator_NEW &grid,
    const double time_beg,
    const double time_end,
    const std::vector<double> &Ez_trapz_absolute,
    const std::vector<double> &Ey_trapz_absolute,
    const std::vector<double> &Ex_trapz_absolute,
    double dt
){

    /* Gather all the is_steady_state_for_this_MPI booleans */
        // MPI cannot communicate boolean, use unsigned int instead:
    std::vector<unsigned int> steady_state_per_MPI(
        grid.MPI_communicator.getNumberOfMPIProcesses(),
        0
    );
        // Transform the boolean value in an unsigned int:
    unsigned int my_steady_state_boolean = 0;

    if(is_steady_state_for_this_MPI){
        my_steady_state_boolean = 1;
    }
    MPI_Allgather(
        &my_steady_state_boolean,
        1,
        MPI_UNSIGNED,
        &steady_state_per_MPI[0],
        1,
        MPI_UNSIGNED,
        MPI_COMM_WORLD
    );
    unsigned int sum = 0;
    for(size_t i = 0 ; i < steady_state_per_MPI.size() ; i ++)
        sum += steady_state_per_MPI[i];

    if(sum != steady_state_per_MPI.size()){
        // At least one MPI is not in the steady state.
        *is_steady_state_for_all_mpi = false;
        return false;
    }

    /** Compute the amplitude of the signal at each node **/
    // Integral(numerical) = Amplitude * integral(abs(sin(2*pi*f*t)) dt) 
    //                     = Amplitude * Integral(analytical)
    // The analytical integral is computed with GaussLobatto method.
    
    
    /// Compute the integral of the absolute value of a sine wave:
    double frequency = grid.input_parser.source.frequency[0];
    printf("\n>>> time_beg %.20g\n>>> time_end %.20g\n2*M_PI*freq %.20g\n",
        time_beg,time_end,2*M_PI*frequency);
    double res = GaussLobattoInt(&absolute_value_sinus_func,
					time_beg,time_end,
                    2*M_PI*frequency,
					1e-10,
					1E5);

    size_t index;

    double average_amplitude = 0.0;

    size_t counter = 0;

    ////////////////////////////////////////////////////////////////
    // WE DIRECTLY PUT THE AMPLITUDE IN THE ELECTRIC FIELD VECTOR //
    ////////////////////////////////////////////////////////////////

    // Electric field along Z:
    for(size_t K = 1 ; K < grid.size_Ez[2]-1 ; K ++){
        for(size_t J = 1 ; J < grid.size_Ez[1]-1 ; J ++){
            for(size_t I = 1 ; I < grid.size_Ez[0]-1 ; I ++){

                index = I + grid.size_Ez[0] * ( J + grid.size_Ez[1] * K);

                if(index >= grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2]){
                    DISPLAY_ERROR_ABORT_CLASS(
                        "Index out of bound..."
                    );
                }

                grid.E_z[index] = Ez_trapz_absolute[index] / res * dt;

                if(I == 2 && J == 2 && K == 2){
                    printf( "\n>>> res     %.20g"
                            "\n>>> Ez_trap %.20g"
                            "\n>>> Amplit  %.20g\n",
                            res,
                            Ez_trapz_absolute[index]*dt,
                            grid.E_z[index]);
                }

                if(counter == 0){
                    average_amplitude = grid.E_z[index];
                }else{
                    average_amplitude = average_amplitude * (counter-1) + grid.E_z[index];
                    average_amplitude /= counter;
                }
                counter++;

            }
        }
    }

    printf("\n\t>>> Average ampltide of Ez is %.15g.\n\n",average_amplitude);

    average_amplitude = 0.0;

    counter = 0;

    // Electric field along X:
    for(size_t K = 1 ; K < grid.size_Ex[2]-1 ; K ++){
        for(size_t J = 1 ; J < grid.size_Ex[1]-1 ; J ++){
            for(size_t I = 1 ; I < grid.size_Ex[0]-1 ; I ++){

                index = I + grid.size_Ex[0] * ( J + grid.size_Ex[1] * K);

                if(index >= grid.size_Ex[2]*grid.size_Ex[2]*grid.size_Ex[2]){
                    DISPLAY_ERROR_ABORT_CLASS(
                        "Index out of bound..."
                    );
                }

                grid.E_x[index] = Ex_trapz_absolute[index] / res * dt;

                if(I == 2 && J == 2 && K == 2){
                    printf( "\n>>> res     %.20g"
                            "\n>>> Ex_trap %.20g"
                            "\n>>> Amplit  %.20g\n",
                            res,
                            Ex_trapz_absolute[index]*dt,
                            grid.E_x[index]);
                }

                if(counter == 0){
                    average_amplitude = grid.E_x[index];
                }else{
                    average_amplitude = average_amplitude * (counter-1) + grid.E_x[index];
                    average_amplitude /= counter;
                }
                counter++;

            }
        }
    }

    printf("\n\t>>> Average ampltide of Ex is %.15g.\n\n",average_amplitude);

    average_amplitude = 0.0;

    counter = 0;

    // Electric field along X:
    for(size_t K = 1 ; K < grid.size_Ey[2]-1 ; K ++){
        for(size_t J = 1 ; J < grid.size_Ey[1]-1 ; J ++){
            for(size_t I = 1 ; I < grid.size_Ey[0]-1 ; I ++){

                index = I + grid.size_Ey[0] * ( J + grid.size_Ey[1] * K);

                if(index >= grid.size_Ey[2]*grid.size_Ey[2]*grid.size_Ey[2]){
                    DISPLAY_ERROR_ABORT_CLASS(
                        "Index out of bound..."
                    );
                }

                grid.E_y[index] = Ey_trapz_absolute[index] / res * dt;

                if(I == 2 && J == 2 && K == 2){
                    printf( "\n>>> res     %.20g"
                            "\n>>> Ey_trap %.20g"
                            "\n>>> Amplit  %.20g\n",
                            res,
                            Ey_trapz_absolute[index]*dt,
                            grid.E_y[index]);
                }

                if(counter == 0){
                    average_amplitude = grid.E_y[index];
                }else{
                    average_amplitude = average_amplitude * (counter-1) + grid.E_y[index];
                    average_amplitude /= counter;
                }
                counter++;

            }
        }
    }

    printf("\n\t>>> Average ampltide of Ey is %.15g.\n\n",average_amplitude);


    *is_steady_state_for_all_mpi = true;
    return true;
}





double AlgoElectro_NEW::interpolationX(size_t x1, size_t x2, size_t x3,
                                      size_t x4, size_t x5, size_t x6,
                                      size_t x7, size_t x8,  GridCreator_NEW &grid)
{
    
    double valueInterpolationX = 0.0;

    size_t max_size_Ex = grid.size_Ex[0]*grid.size_Ex[1]*grid.size_Ex[2];

    if(x1 >= max_size_Ex){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ex avec x1."
        );
    }
    if(x2 >= max_size_Ex){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ex avec x2."
        );
    }
    if(x3 >= max_size_Ex){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ex avec x3."
        );
    }
    if(x4 >= max_size_Ex){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ex avec x4."
        );
    }
    if(x5 >= max_size_Ex){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ex avec x5."
        );
    }
    if(x6 >= max_size_Ex){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ex avec x6."
        );
    }
    if(x7 >= max_size_Ex){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ex avec x7."
        );
    }
    if(x8 >= max_size_Ex){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ex avec x8."
        );
    }

    double c_000_Ex = grid.E_x[x1];
    double c_001_Ex = grid.E_x[x2];
    double c_010_Ex = grid.E_x[x3];
    double c_011_Ex = grid.E_x[x4];
    double c_100_Ex = grid.E_x[x5];
    double c_101_Ex = grid.E_x[x6];
    double c_110_Ex = grid.E_x[x7];
    double c_111_Ex = grid.E_x[x8];

    double Mean15_x = (c_000_Ex + c_100_Ex) / 2;
    double Mean26_x = (c_001_Ex + c_101_Ex) / 2;
    double Mean37_x = (c_010_Ex + c_110_Ex) / 2;
    double Mean48_x = (c_011_Ex + c_111_Ex) / 2;

    double Mean1537_x = (Mean15_x + Mean37_x) / 2;
    double Mean2648_x = (Mean26_x + Mean48_x) / 2;

    valueInterpolationX = (Mean1537_x + Mean2648_x) / 2;

    return valueInterpolationX;
}


double AlgoElectro_NEW::interpolationY(size_t y1, size_t y2, size_t y3,
                                      size_t y4, size_t y5, size_t y6,
                                      size_t y7, size_t y8,  GridCreator_NEW &grid)
{
    
    double valueInterpolationY = 0.0;
    
    size_t max_size_Ey = grid.size_Ey[0]*grid.size_Ey[1]*grid.size_Ey[2];

    if(y1 >= max_size_Ey){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ey avec y1."
        );
    }
    if(y2 >= max_size_Ey){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ey avec y2."
        );
    }
    if(y3 >= max_size_Ey){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ey avec y3."
        );
    }
    if(y4 >= max_size_Ey){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ey avec y4."
        );
    }
    if(y5 >= max_size_Ey){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ey avec y5."
        );
    }
    if(y6 >= max_size_Ey){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ey avec y6."
        );
    }
    if(y7 >= max_size_Ey){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ey avec y7."
        );
    }
    if(y8 >= max_size_Ey){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ey avec y8."
        );
    }
    
    
    // Will contain the value of Ex along the 3 direction of space that will be used to compute the modulus
    // double *EyForModulus = (double *)calloc(3, sizeof(double));
    
    // if(EyForModulus == NULL)
    // {
    //     printf("The pointer EyForModulus could not be allocated.\n");
    //     printf("This error comes from line %d from file %s\n", __LINE__, __FILE__);
    //     printf("Aborting...");
    //     exit(EXIT_FAILURE);
    // }

    double c_000_Ey = grid.E_y[y1];
    double c_001_Ey = grid.E_y[y2];
    double c_010_Ey = grid.E_y[y3];
    double c_011_Ey = grid.E_y[y4];
    double c_100_Ey = grid.E_y[y5];
    double c_101_Ey = grid.E_y[y6];
    double c_110_Ey = grid.E_y[y7];
    double c_111_Ey = grid.E_y[y8];

    double Mean15_y = (c_000_Ey + c_100_Ey) / 2;
    double Mean26_y = (c_001_Ey + c_101_Ey) / 2;
    double Mean37_y = (c_010_Ey + c_110_Ey) / 2;
    double Mean48_y = (c_011_Ey + c_111_Ey) / 2;

    double Mean1537_y = (Mean15_y + Mean37_y) / 2;
    double Mean2648_y = (Mean26_y + Mean48_y) / 2;

    valueInterpolationY = (Mean1537_y + Mean2648_y) / 2;

    return valueInterpolationY;
}


double AlgoElectro_NEW::interpolationZ(size_t z1, size_t z2, size_t z3,
                                      size_t z4, size_t z5, size_t z6,
                                      size_t z7, size_t z8,  GridCreator_NEW &grid)
{
    
    double valueInterpolationZ = 0.0;

    size_t max_size_Ez = grid.size_Ez[0]*grid.size_Ez[1]*grid.size_Ez[2];

    if(z1 >= max_size_Ez){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ez avec z1."
        );
    }
    if(z2 >= max_size_Ez){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ez avec z2."
        );
    }
    if(z3 >= max_size_Ez){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ez avec z3."
        );
    }
    if(z4 >= max_size_Ez){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ez avec z4."
        );
    }
    if(z5 >= max_size_Ez){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ez avec z5."
        );
    }
    if(z6 >= max_size_Ez){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ez avec z6."
        );
    }
    if(z7 >= max_size_Ez){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ez avec z7."
        );
    }
    if(z8 >= max_size_Ez){
        DISPLAY_ERROR_ABORT_CLASS(
            "Dépassement de tabeau pour Ez avec z8."
        );
    }
    
    double c_000_Ez = grid.E_z[z1];
    double c_001_Ez = grid.E_z[z2];
    double c_010_Ez = grid.E_z[z3];
    double c_011_Ez = grid.E_z[z4];
    double c_100_Ez = grid.E_z[z5];
    double c_101_Ez = grid.E_z[z6];
    double c_110_Ez = grid.E_z[z7];
    double c_111_Ez = grid.E_z[z8];

    double Mean15_z = (c_000_Ez + c_100_Ez) / 2;
    double Mean26_z = (c_001_Ez + c_101_Ez) / 2;
    double Mean37_z = (c_010_Ez + c_110_Ez) / 2;
    double Mean48_z = (c_011_Ez + c_111_Ez) / 2;

    double Mean1537_z = (Mean15_z + Mean37_z) / 2;
    double Mean2648_z = (Mean26_z + Mean48_z) / 2;

    valueInterpolationZ = (Mean1537_z + Mean2648_z) / 2;

    return valueInterpolationZ;
}





std::vector<double> AlgoElectro_NEW::ComputeNormE2square(std::vector<double> Ex,
                                                   std::vector<double> Ey, std::vector<double> Ez)
{
    // Those indices will serve to go through all the nodes of the domain
    size_t i = 0;
    size_t j = 0;
    size_t k = 0;

    // This vector will contain the norm of the electric field at each node
    std::vector<double> ModulusE;

    // Number of centers for Ex along the 3 directions of space
    size_t NbCentersExx = 1;
    size_t NbCentersExy = 1;
    size_t NbCentersExz = 1;

    // Number of centers for Ey along the 3 directions of space
    size_t NbCentersEyx = 1;
    size_t NbCentersEyy = 1;
    size_t NbCentersEyz = 1;

    // Number of centers for Ez along the 3 directions of space
    size_t NbCentersEzx = 1;
    size_t NbCentersEzy = 1;
    size_t NbCentersEzz = 1;

    size_t x1 = 0; // Correspond to point (i+1, j-1, k-1)
    size_t x2 = 0; // Correspond to point (i-1, j-1, k-1)
    size_t x3 = 0; // Correspond to point (i+1, j-1, k+1)
    size_t x4 = 0; // Correspond to point (i-1, j-1, k+1)
    size_t x5 = 0; // Correspond to point (i+1, j+1, k-1)
    size_t x6 = 0; // Correspond to point (i-1, j+1, k-1)
    size_t x7 = 0; // Correspond to point (i+1, j+1, k-1)
    size_t x8 = 0; // Correspond to point (i+1, j+1, k+1)

    size_t y1 = 0; // Correspond to point (i+1, j-1, k-1)
    size_t y2 = 0; // Correspond to point (i-1, j-1, k-1)
    size_t y3 = 0; // Correspond to point (i+1, j-1, k+1)
    size_t y4 = 0; // Correspond to point (i-1, j-1, k+1)
    size_t y5 = 0; // Correspond to point (i+1, j+1, k-1)
    size_t y6 = 0; // Correspond to point (i-1, j+1, k-1)
    size_t y7 = 0; // Correspond to point (i+1, j+1, k-1)
    size_t y8 = 0; // Correspond to point (i+1, j+1, k+1)

    size_t z1 = 0; // Correspond to point (i+1, j-1, k-1)
    size_t z2 = 0; // Correspond to point (i-1, j-1, k-1)
    size_t z3 = 0; // Correspond to point (i+1, j-1, k+1)
    size_t z4 = 0; // Correspond to point (i-1, j-1, k+1)
    size_t z5 = 0; // Correspond to point (i+1, j+1, k-1)
    size_t z6 = 0; // Correspond to point (i-1, j+1, k-1)
    size_t z7 = 0; // Correspond to point (i+1, j+1, k-1)
    size_t z8 = 0; // Correspond to point (i+1, j+1, k+1)

    // Will contain the interpolation for a given cube
    double interpolationEx = 0.0;
    double interpolationEy = 0.0;
    double interpolationEz = 0.0;

    // Will contain the interpolation for each cube
    std::vector<double> allInterpolationEx;
    std::vector<double> allInterpolationEy;
    std::vector<double> allInterpolationEz;

    printf("Ex.size = %zu\n", Ex.size());
    //size_t Nx = 1;
    //size_t Ny = 1;
    //size_t Nz = 1;
    // Computations for Ex
    for(i=0; i<NbCentersExx; i++)
    {
        for(j=0; j<NbCentersExy; j++)
        {
            for(k=0; k<NbCentersExz; k++)
            {
                // x1 = (i+1) + Nx * (j + k * Ny);
                // x2 = i + Nx * (j + k * Ny);
                // x3 = (i+1) + Nx * (j + (k+1) * Ny);
                // x4 = i + Nx * (j + (k+1) * Ny);
                // x5 = i+1 + j + 2*Nx + k;
                // x6 = i + Nx * ((j+1) + k * Ny);
                // x7 = (i+1) + Nx * ((j+1) + k * Ny);
                // x8 = (i+1) + Nx * ((j+1) + (k+1) * Ny);
                x1 = 0;
                x2 = 1;
                x4 = 5;
                x3 = 4;
                x6 = 3;
                x5 = 2;
                x8 = 7;
                x7 = 6;

                interpolationEx = interpolationX2(x1, x2, x3, x4, x5, x6, x7, x8, Ex);

                printf("Valuex = %lf\n", interpolationEx);

                allInterpolationEx.push_back(interpolationEx);
            }
        }
    }

    // Computations for Ey
    for(i=0; i<NbCentersEyx; i++)
    {
        for(j=0; j<NbCentersEyy; j++)
        {
            for(k=0; k<NbCentersEyz; k++)
            {
                // y1 = (i+1) + Ey.size() * (j + k * Ey.size());
                // y2 = i + Ey.size() * (j + k * Ey.size());
                // y3 = (i+1) + Ey.size() * (j + (k+1) * Ey.size());
                // y4 = i + Ey.size() * (j + (k+1) * Ey.size());
                // y5 = (i+1) + Ey.size() * ((j+1) + k * Ey.size());
                // y6 = i + Ey.size() * ((j+1) + k * Ey.size());
                // y7 = (i+1) + Ey.size() * ((j+1) + k * Ey.size());
                // y8 = (i+1) + Ey.size() * ((j+1) + (k+1) * Ey.size());
                y1 = 1;
                y2 = 0;
                y3 = 5;
                y4 = 4;
                y5 = 3;
                y6 = 2;
                y7 = 7;
                y8 = 6;

                interpolationEy = interpolationY2(y1, y2, y3, y4, y5, y6, y7, y8, Ey);

                printf("Valuey = %lf\n", interpolationEy);

                allInterpolationEy.push_back(interpolationEy);
            }
        }
    }


    // Computations for Ez
    for(i=0; i<NbCentersEzx; i++)
    {
        for(j=0; j<NbCentersEzy; j++)
        {
            for(k=0; k<NbCentersEzz; k++)
            {
                // z1 = (i+1) + Ez.size() * (j + k * Ez.size());
                // z2 = i + Ez.size() * (j + k * Ez.size());
                // z3 = (i+1) + Ez.size() * (j + (k+1) * Ez.size());
                // z4 = i + Ez.size() * (j + (k+1) * Ez.size());
                // z5 = (i+1) + Ez.size() * ((j+1) + k * Ez.size());
                // z6 = i + Ez.size() * ((j+1) + k * Ez.size());
                // z7 = (i+1) + Ez.size() * ((j+1) + k * Ez.size());
                // z8 = (i+1) + Ez.size() * ((j+1) + (k+1) * Ez.size());
                z1 = 1;
                z2 = 0;
                z3 = 5;
                z4 = 4;
                z5 = 3;
                z6 = 2;
                z7 = 7;
                z8 = 6;

                interpolationEz = interpolationZ2(z1, z2, z3, z4, z5, z6, z7, z8, Ez);

                allInterpolationEz.push_back(interpolationEz);
                printf("Valuez = %lf\n", interpolationEz);
            }
        }
    }

    size_t lengthInterpX = 0;
    size_t lengthInterpY = 0;
    size_t lengthInterpZ = 0;

    for(lengthInterpX=0; lengthInterpX<allInterpolationEx.size(); lengthInterpX++)
    {
        for(lengthInterpY=0; lengthInterpY<allInterpolationEy.size(); lengthInterpY++)
        {
            for(lengthInterpZ=0; lengthInterpZ<allInterpolationEz.size(); lengthInterpZ++)
            {
                // Recall : modulus = sqrt(x1^2 + x2^2 + x3^2)
                ModulusE.push_back(allInterpolationEx[lengthInterpX]*allInterpolationEx[lengthInterpX]
                                   + allInterpolationEy[lengthInterpY]*allInterpolationEy[lengthInterpY]
                                   + allInterpolationEz[lengthInterpZ]*allInterpolationEz[lengthInterpZ]);
            }
        }
    }

    int iiiiii=0;
    for(iiiiii=0; iiiiii<1; iiiiii++)
        printf("Modulus = %lf\n", ModulusE[iiiiii]);
    
    
    return ModulusE;
}


double AlgoElectro_NEW::interpolationX2(size_t x1, size_t x2, size_t x3,
                        size_t x4, size_t x5, size_t x6,
                        size_t x7, size_t x8,  std::vector<double> Ex)
{
    double valueInterpolationX = 0.0;

    // printf("In interpolationX2 \n");
    // printf("x1 = %zu\n", x1);
    // printf("x2 = %zu\n", x2);
    // printf("x3 = %zu\n", x3);
    // printf("x4 = %zu\n", x4);
    // printf("x5 = %zu\n", x5);
    // printf("x6 = %zu\n", x6);
    // printf("x7 = %zu\n", x7);
    // printf("x8 = %zu\n", x8);

    double c_000_Ex = Ex[x1];
    double c_001_Ex = Ex[x2];
    double c_010_Ex = Ex[x3];
    double c_011_Ex = Ex[x4];
    double c_100_Ex = Ex[x5];
    double c_101_Ex = Ex[x6];
    double c_110_Ex = Ex[x7];
    double c_111_Ex = Ex[x8];
    printf("x1 = %lf\n", Ex[x1]);
    printf("x2 = %lf\n", Ex[x2]);
    printf("x3 = %lf\n", Ex[x3]);
    printf("x4 = %lf\n", Ex[x4]);
    printf("x5 = %lf\n", Ex[x5]);
    printf("x6 = %lf\n", Ex[x6]);
    printf("x7 = %lf\n", Ex[x7]);
    printf("x8 = %lf\n", Ex[x8]);

    double Mean15_x = (c_000_Ex + c_100_Ex) / 2;
    double Mean26_x = (c_001_Ex + c_101_Ex) / 2;
    double Mean37_x = (c_010_Ex + c_110_Ex) / 2;
    double Mean48_x = (c_011_Ex + c_111_Ex) / 2;

    double Mean1537_x = (Mean15_x + Mean37_x) / 2;
    double Mean2648_x = (Mean26_x + Mean48_x) / 2;

    valueInterpolationX = (Mean1537_x + Mean2648_x) / 2;

    printf("Mean15x = %lf\n", Mean15_x);
    printf("Mean26x = %lf\n", Mean26_x);
    printf("Mean37x = %lf\n", Mean37_x);
    printf("Mean48x = %lf\n", Mean48_x);
    printf("Mean1537x = %lf\n", Mean1537_x);
    printf("Mean2648x = %lf\n", Mean2648_x);
    printf("valueInterpolationX = %lf\n", valueInterpolationX);
    
    return valueInterpolationX;
}

double AlgoElectro_NEW::interpolationY2(size_t y1, size_t y2, size_t y3,
                        size_t y4, size_t y5, size_t y6,
                        size_t y7, size_t y8,  std::vector<double> Ey)
{
     double valueInterpolationY = 0.0;
    
    double c_000_Ey = Ey[y1];
    double c_001_Ey = Ey[y2];
    double c_010_Ey = Ey[y3];
    double c_011_Ey = Ey[y4];
    double c_100_Ey = Ey[y5];
    double c_101_Ey = Ey[y6];
    double c_110_Ey = Ey[y7];
    double c_111_Ey = Ey[y8];

    double Mean15_y = (c_000_Ey + c_100_Ey) / 2;
    double Mean26_y = (c_001_Ey + c_101_Ey) / 2;
    double Mean37_y = (c_010_Ey + c_110_Ey) / 2;
    double Mean48_y = (c_011_Ey + c_111_Ey) / 2;

    double Mean1537_y = (Mean15_y + Mean37_y) / 2;
    double Mean2648_y = (Mean26_y + Mean48_y) / 2;

    valueInterpolationY = (Mean1537_y + Mean2648_y) / 2;

    printf("Mean15y = %lf\n", Mean15_y);
    printf("Mean26y = %lf\n", Mean26_y);
    printf("Mean37y = %lf\n", Mean37_y);
    printf("Mean48y = %lf\n", Mean48_y);
    printf("Mean1537y = %lf\n", Mean1537_y);
    printf("Mean2648y = %lf\n", Mean2648_y);
    printf("valueInterpolationY = %lf\n", valueInterpolationY);

    return valueInterpolationY;
}

double AlgoElectro_NEW::interpolationZ2(size_t z1, size_t z2, size_t z3,
                        size_t z4, size_t z5, size_t z6,
                        size_t z7, size_t z8,  std::vector<double> Ez)
{
    double valueInterpolationZ = 0.0;
    
    double c_000_Ez = Ez[z1];
    double c_001_Ez = Ez[z2];
    double c_010_Ez = Ez[z3];
    double c_011_Ez = Ez[z4];
    double c_100_Ez = Ez[z5];
    double c_101_Ez = Ez[z6];
    double c_110_Ez = Ez[z7];
    double c_111_Ez = Ez[z8];

    double Mean15_z = (c_000_Ez + c_100_Ez) / 2;
    double Mean26_z = (c_001_Ez + c_101_Ez) / 2;
    double Mean37_z = (c_010_Ez + c_110_Ez) / 2;
    double Mean48_z = (c_011_Ez + c_111_Ez) / 2;

    double Mean1537_z = (Mean15_z + Mean37_z) / 2;
    double Mean2648_z = (Mean26_z + Mean48_z) / 2;

    valueInterpolationZ = (Mean1537_z + Mean2648_z) / 2;

    printf("Mean15z = %lf\n", Mean15_z);
    printf("Mean26z = %lf\n", Mean26_z);
    printf("Mean37z = %lf\n", Mean37_z);
    printf("Mean48z = %lf\n", Mean48_z);
    printf("Mean1537z = %lf\n", Mean1537_z);
    printf("Mean2648z = %lf\n", Mean2648_z);
    printf("valueInterpolationZ = %lf\n", valueInterpolationZ);

    return valueInterpolationZ;
}




std::vector<double> AlgoElectro_NEW::ComputeNormEsquareBIS(GridCreator_NEW &grid)
{
    if(grid.MPI_communicator.getNumberOfMPIProcesses() != 1)
    {
        DISPLAY_ERROR_ABORT("More than one process is used. The functino is intended for only one");
    }

    size_t M = grid.sizes_EH[0];
    size_t N = grid.sizes_EH[1];
    size_t P = grid.sizes_EH[2];

    printf("M = %zu, N = %zu, P = %zu\n", M, N, P);

    // Will serve to go through the grid
    size_t i,j,k;

    // This vector will contain the norm of the electric field at each node
    std::vector<double> SquareModulusE;

    // Will contain the interpolation for a given cube
    double interpolationEx = 0.0;
    double interpolationEy = 0.0;
    double interpolationEz = 0.0;

    // Will contain the interpolation for each cube
    std::vector<double> allInterpolationEx;
    std::vector<double> allInterpolationEy;
    std::vector<double> allInterpolationEz;

    // Will serve to have the 8 corners of the cube
    size_t x1 = 0;
    size_t x2 = 0;
    size_t x3 = 0;
    size_t x4 = 0;
    size_t x5 = 0;
    size_t x6 = 0;
    size_t x7 = 0;
    size_t x8 = 0;

    // Will serve to have the 8 corners of the cube
    size_t y1 = 0;
    size_t y2 = 0;
    size_t y3 = 0;
    size_t y4 = 0;
    size_t y5 = 0;
    size_t y6 = 0;
    size_t y7 = 0;
    size_t y8 = 0;

    // Will serve to have the 8 corners of the cube
    size_t z1 = 0;
    size_t z2 = 0;
    size_t z3 = 0;
    size_t z4 = 0;
    size_t z5 = 0;
    size_t z6 = 0;
    size_t z7 = 0;
    size_t z8 = 0;

    size_t count_ex = 0;
    size_t count_ey = 0;
    size_t count_ez = 0;

    
    // Computation for Ex
    for(i=0; i<M-2; i++)
    {
        for(j=0; j<N-2; j++)
        {
            for(k=0; k<P-2; k++)
            {
                x1 = (i+1) + grid.size_Ex[0] * (j + k * grid.size_Ex[1]);           // Correspond to point (i+1, j, k)
                x2 = i     + grid.size_Ex[0] * (j + k * grid.size_Ex[1]);           // Correspond to point (i, j, k)
                x3 = (i+1) + grid.size_Ex[0] * (j + (k+1) * grid.size_Ex[1]);       // Correspond to point (i+1, j, k+1)
                x4 = i     + grid.size_Ex[0] * ((j+1) + k * grid.size_Ex[1]);       // Correspond to point (i, j+1, k)
                x5 = (i+1) + grid.size_Ex[0] * ((j+1) + k * grid.size_Ex[1]);       // Correspond to point (i+1, j+1, k)
                x6 = i     + grid.size_Ex[0] * ((j+1) + k * grid.size_Ex[1]);           // Correspond to point (i, j+1, k)
                x7 = (i+1) + grid.size_Ex[0] * ((j+1) + (k+1) * grid.size_Ex[1]);   // Correspond to point (i+1, j+1, k+1)
                x8 = i     + grid.size_Ex[0] * ((j+1) + (k+1) * grid.size_Ex[1]);       // Correspond to point (i, j+1, k+1)

                interpolationEx = interpolationX(x1, x2, x3, x4, x5, x6, x7, x8, grid);

                allInterpolationEx.push_back(interpolationEx);
                count_ex++;
            }
        }
    }

    // Computation for Ey
    for(i=0; i<M-2; i++)
    {
        for(j=0; j<N-2; j++)
        {
            for(k=0; k<P-2; k++)
            {
                y1 = (i+1) + grid.size_Ey[0] * (j + k * grid.size_Ey[1]);           // Correspond to point (i+1, j, k)
                y2 = i     + grid.size_Ey[0] * (j + k * grid.size_Ey[1]);               // Correspond to point (i, j, k)
                y3 = (i+1) + grid.size_Ey[0] * (j + (k+1) * grid.size_Ey[1]);       // Correspond to point (i+1, j, k+1)
                y4 = i     + grid.size_Ey[0] * ((j+1) + k * grid.size_Ey[1]);           // Correspond to point (i, j+1, k)
                y5 = (i+1) + grid.size_Ey[0] * ((j+1) + k * grid.size_Ey[1]);       // Correspond to point (i+1, j+1, k)
                y6 = i     + grid.size_Ey[0] * ((j+1) + k * grid.size_Ey[1]);           // Correspond to point (i, j+1, k)
                y7 = (i+1) + grid.size_Ey[0] * ((j+1) + (k+1) * grid.size_Ey[1]);   // Correspond to point (i+1, j+1, k+1)
                y8 = i     + grid.size_Ey[0] * ((j+1) + (k+1) * grid.size_Ey[1]);       // Correspond to point (i, j+1, k+1)

                interpolationEy = interpolationY(y1, y2, y3, y4, y5, y6, y7, y8, grid);

                allInterpolationEy.push_back(interpolationEy);
                count_ey++;
            }
        }
    }

    // Computation for Ez
    for(i=0; i<M-2; i++)
    {
        for(j=0; j<N-2; j++)
        {
            for(k=0; k<P-2; k++)
            {
                z1 = (i+1) + grid.size_Ez[0] * (j + k * grid.size_Ez[1]);           // Correspond to point (i+1, j, k)
                z2 = i     + grid.size_Ez[0] * (j + k * grid.size_Ez[1]);               // Correspond to point (i, j, k)
                z3 = (i+1) + grid.size_Ez[0] * (j + (k+1) * grid.size_Ez[1]);       // Correspond to point (i+1, j, k+1)
                z4 = i     + grid.size_Ez[0] * ((j+1) + k * grid.size_Ez[1]);           // Correspond to point (i, j+1, k)
                z5 = (i+1) + grid.size_Ez[0] * ((j+1) + k * grid.size_Ez[1]);       // Correspond to point (i+1, j+1, k)
                z6 = i     + grid.size_Ez[0] * ((j+1) + k * grid.size_Ez[1]);           // Correspond to point (i, j+1, k)
                z7 = (i+1) + grid.size_Ez[0] * ((j+1) + (k+1) * grid.size_Ez[1]);   // Correspond to point (i+1, j+1, k+1)
                z8 = i     + grid.size_Ez[0] * ((j+1) + (k+1) * grid.size_Ez[1]);       // Correspond to point (i, j+1, k+1)

                interpolationEz = interpolationZ(z1, z2, z3, z4, z5, z6, z7, z8, grid);

                allInterpolationEz.push_back(interpolationEz);
                count_ez++;
            }
        }
    }

    printf("There were %zu loop iterations for the Ex interpolations\n", count_ex);
    printf("There were %zu loop iterations for the Ey interpolations\n", count_ey);
    printf("There were %zu loop iterations for the Ez interpolations\n", count_ez);

    size_t Centers = 0;

    for(Centers = 0; Centers < (M-2)*(N-2)*(P-2); Centers++)
    {
        if(Centers >= allInterpolationEx.size()){
            DISPLAY_ERROR_ABORT_CLASS(
                "Dépassement de tableau sur allInterpolationEx."
            );
        }
        if(Centers >= allInterpolationEy.size()){
            DISPLAY_ERROR_ABORT_CLASS(
                "Dépassement de tableau sur allInterpolationEy."
            );
        }
        if(Centers >= allInterpolationEz.size()){
            DISPLAY_ERROR_ABORT_CLASS(
                "Dépassement de tableau sur allInterpolationEz."
            );
        }
        SquareModulusE.push_back(allInterpolationEx[Centers]*allInterpolationEx[Centers]
                                + allInterpolationEy[Centers]*allInterpolationEy[Centers]
                                + allInterpolationEz[Centers]*allInterpolationEz[Centers]);
    }

    printf("Size of the vector with the Ex interpolations = %zu\n", allInterpolationEx.size());
    printf("Size of the vector with the Ex interpolations = %zu\n", allInterpolationEy.size());
    printf("Size of the vector with the Ex interpolations = %zu\n", allInterpolationEz.size());
    printf("Size of the vector with all the moduli = %zu\n", SquareModulusE.size());
    // abort();

    
    return SquareModulusE;
}


/**
 * @brief Apply 1D conditions to the magnetic field.
 */
void AlgoElectro_NEW::apply_1D_case_on_magnetic_field(
    bool IS_1D_FACE_EX_Electric_along_Z  ,
    bool IS_1D_FACE_EX_Electric_along_Y ,
    bool IS_1D_FACE_EY_Electric_along_Z ,
    bool IS_1D_FACE_EY_Electric_along_X  ,
    bool IS_1D_FACE_EZ_Electric_along_Y , 
    bool IS_1D_FACE_EZ_Electric_along_X  ,
    bool IS_1D_FACE_Minus_EX_Electric_along_Z ,
    bool IS_1D_FACE_Minus_EX_Electric_along_Y,
    bool IS_1D_FACE_Minus_EY_Electric_along_Z ,
    bool IS_1D_FACE_Minus_EY_Electric_along_X ,
    bool IS_1D_FACE_Minus_EZ_Electric_along_X ,
    bool IS_1D_FACE_Minus_EZ_Electric_along_Y ,
    GridCreator_NEW &grid,
    double *H_x_tmp,
    double *H_y_tmp,
    double *H_z_tmp
)
{
    // Check the number of MPI processes:
    if(grid.MPI_communicator.getNumberOfMPIProcesses() != 1){
        printf("For the 1D case, use only 1 MPI process.\n");
        DISPLAY_ERROR_ABORT_CLASS(
            "Only one MPI process is allowed for the 1D cases."
        );
    }

    // Check that only one boolean is true:
    size_t counterBooleanTrue = 0;
    if(       IS_1D_FACE_EX_Electric_along_Y == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EX_Electric_along_Z == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EX_Electric_along_Z == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EX_Electric_along_Y == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EY_Electric_along_X == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EY_Electric_along_Z == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EY_Electric_along_X == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EY_Electric_along_Z == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EZ_Electric_along_X == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EZ_Electric_along_Y == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EZ_Electric_along_X == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EZ_Electric_along_Y == true){
        counterBooleanTrue ++;
    }
    if(counterBooleanTrue != 1){
        DISPLAY_ERROR_ABORT_CLASS(
            "Should have only one true boolean but has %zu.",
            counterBooleanTrue
        );
    }

    if(IS_1D_FACE_Minus_EY_Electric_along_Z){       
        /**
         * What will be imposed:
         *  1) On faces with normal e_x or -e_x:
         *          - Hy is zero.
         */
        size_t i, index;
        printf("\t> > 1D case is IS_1D_FACE_Minus_EY.\n");
        printf("\t\t> Imposing HY on faces with normal +e_x and -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Hy[1] - 1 ; j ++ ){
            for(size_t k = 1 ; k < grid.size_Hy[2] - 1 ; k ++){

                i = 1;
                index = i + grid.size_Hy[0] * ( j + grid.size_Hy[1] * k);
                H_y_tmp[index] = 0;

                i = grid.size_Hy[0] - 2;
                index = i + grid.size_Hy[0] * ( j + grid.size_Hy[1] * k);
                H_y_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EY_Electric_along_Z){
        /**
         * What will be imposed:
         *  1) On faces with normal e_x or -e_x:
         *          - Hy is zero.
         */
        size_t i, index;
        printf("\t> > 1D case is IS_1D_FACE_EY_Electric_along_Z.\n");
        printf("\t\t> Imposing HY on faces with normal +e_x and -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Hy[1] - 1 ; j ++ ){
            for(size_t k = 1 ; k < grid.size_Hy[2] - 1 ; k ++){

                i = 1;
                index = i + grid.size_Hy[0] * ( j + grid.size_Hy[1] * k);
                H_y_tmp[index] = 0;

                i = grid.size_Hy[0] - 2;
                index = i + grid.size_Hy[0] * ( j + grid.size_Hy[1] * k);
                H_y_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EY_Electric_along_X){
        /**
         * What will be imposed:
         *  1) On faces with normal e_x or -e_x:
         *          - Hy is zero.
         */
        size_t k, index;
        printf("\t> > 1D case is IS_1D_FACE_EY.\n");
        printf("\t\t> Imposing HY on faces with normal +e_x and -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Hy[1] - 1 ; j ++ ){
            for(size_t i = 1 ; i < grid.size_Hy[0] - 1 ; i ++){

                k = 1;
                index = i + grid.size_Hy[0] * ( j + grid.size_Hy[1] * k);
                H_y_tmp[index] = 0;

                k = grid.size_Hy[2] - 2;
                index = i + grid.size_Hy[0] * ( j + grid.size_Hy[1] * k);
                H_y_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_Minus_EX_Electric_along_Z){
        /**
         * What wil be imposed:
         *  1) On faces with normal +e_y or -e_y:
         *          - Hx is zero.
         */
        size_t j, index;
        printf("\t> > 1D case is IS_1D_FACE_Minus_EX_Electric_along_Z.\n");
        printf("\t\t> Imposing HX on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Hx[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Hx[2] - 1 ; k ++){

                j = 1;
                index = i + grid.size_Hx[0] * ( j + grid.size_Hx[1] * k);
                H_x_tmp[index] = 0;

                j = grid.size_Hx[1] - 2 ;
                index = i + grid.size_Hx[0] * ( j + grid.size_Hx[1] * k);
                H_x_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_Minus_EX_Electric_along_Y){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y or -e_y:
         *          - Hx is zero.
         */
        size_t j, index;
        printf("\t> > 1D case is IS_1D_FACE_Minus_EX_Electric_along_Y.\n");
        printf("\t\t> Imposing HX on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Hx[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Hx[2] - 1 ; k ++){

                j = 1;
                index = i + grid.size_Hx[0] * ( j + grid.size_Hx[1] * k);
                //H_x_tmp[index] = 0;

                j = grid.size_Hx[1] - 2 ;
                index = i + grid.size_Hx[0] * ( j + grid.size_Hx[1] * k);
                //H_x_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EX_Electric_along_Z){
        /**
         * What wil be imposed:
         *  1) On faces with normal +e_y or -e_y:
         *          - Hx is zero.
         */
        size_t j, index;
        printf("\t> > 1D case is IS_1D_FACE_EX_Electric_along_Z.\n");
        printf("\t\t> Imposing HX on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Hx[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Hx[2] - 1 ; k ++){

                j = 1;
                index = i + grid.size_Hx[0] * ( j + grid.size_Hx[1] * k);
                H_x_tmp[index] = 0;

                j = grid.size_Hx[1] - 2 ;
                index = i + grid.size_Hx[0] * ( j + grid.size_Hx[1] * k);
                H_x_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EX_Electric_along_Y){
        /**
         * What wil be imposed:
         *  1) On faces with normal +e_z or -e_z:
         *          - Hx is zero.
         */
        size_t k, index;
        printf("\t> > 1D case is IS_1D_FACE_EX_Electric_along_Y.\n");
        printf("\t\t> Imposing HX on faces with normal +e_z and -e_z...\n");
        for(size_t i = 1 ; i < grid.size_Hx[0] - 1 ; i ++){
            for(size_t j = 1 ; j < grid.size_Hx[1] - 1 ; j ++){

                k = 1;
                index = i + grid.size_Hx[0] * ( j + grid.size_Hx[1] * k);
                H_x_tmp[index] = 0;

                k = grid.size_Hx[2] - 2 ;
                index = i + grid.size_Hx[0] * ( j + grid.size_Hx[1] * k);
                H_x_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_Minus_EZ_Electric_along_Y){
        /**
         * What wil be imposed:
         *  1) On faces with normal +e_x or -e_x:
         *          - Hz is zero.
         */
        size_t i, index;
        printf("\t> > 1D case is IS_1D_FACE_Minus_EZ.\n");
        printf("\t\t> Imposing Hz on faces with normal +e_x and -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Hz[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Hz[2] - 1 ; k ++){

                i = 1;
                index = i + grid.size_Hz[0] * ( j + grid.size_Hz[1] * k);
                H_z_tmp[index] = 0;

                i = grid.size_Hz[0] - 2 ;
                index = i + grid.size_Hz[0] * ( j + grid.size_Hz[1] * k);
                H_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EZ_Electric_along_Y){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_x or -e_x:
         *          - Hz is zero.
         */
        size_t i, index;
        printf("\t> > 1D case is IS_1D_FACE_EZ.\n");
        printf("\t\t> Imposing Hz on faces with normal +e_x and -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Hz[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Hz[2] - 1 ; k ++){

                i = 1;
                index = i + grid.size_Hz[0] * ( j + grid.size_Hz[1] * k);
                H_z_tmp[index] = 0;

                i = grid.size_Hz[0] - 2 ;
                index = i + grid.size_Hz[0] * ( j + grid.size_Hz[1] * k);
                H_z_tmp[index] = 0;
            }
        }
    
    }else if(IS_1D_FACE_EZ_Electric_along_X){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_x or -e_x:
         *          - Hz is zero.
         */
        size_t j, index;
        printf("\t> > 1D case is IS_1D_FACE_EZ_Electric_along_X.\n");
        printf("\t\t> Imposing Hz on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Hz[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Hz[2] - 1 ; k ++){

                j = 1;
                index = i + grid.size_Hz[0] * ( j + grid.size_Hz[1] * k);
                H_z_tmp[index] = 0;

                j = grid.size_Hz[0] - 2 ;
                index = i + grid.size_Hz[0] * ( j + grid.size_Hz[1] * k);
                H_z_tmp[index] = 0;
            }
        }
    
    }else if(IS_1D_FACE_Minus_EZ_Electric_along_X){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_x or -e_x:
         *          - Hz is zero.
         */
        size_t j, index;
        printf("\t> > 1D case is IS_1D_FACE_Minus_EZ_Electric_along_X.\n");
        printf("\t\t> Imposing Hz on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Hz[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Hz[2] - 1 ; k ++){

                j = 1;
                index = i + grid.size_Hz[0] * ( j + grid.size_Hz[1] * k);
                H_z_tmp[index] = 0;

                j = grid.size_Hz[0] - 2 ;
                index = i + grid.size_Hz[0] * ( j + grid.size_Hz[1] * k);
                H_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_Minus_EY_Electric_along_X){
        /**
         * Propagation of electric field Ex in the +Y direction.
         * Propagation of magnetic field Hz in the +Y direction.
         * What will be applied:
         *  1) On faces with normal +e_z and -e_z:
         *          - Hy is imposedto zero.
         */
        size_t k, index;
        printf("\t> > 1D case is IS_1D_FACE_Minus_EY_Electric_along_X.\n");
        printf("\t\t> Imposing Hy on faces with normal +e_z and -e_z...\n");
        for(size_t i = 1 ; i < grid.size_Hy[0] - 1 ; i ++){
            for(size_t j = 1 ; j < grid.size_Hy[1] - 1 ; j ++){
                k = 1;
                index = i + grid.size_Hy[0] * ( j + k * grid.size_Hy[1]);
                H_y_tmp[index] = 0;

                k = grid.size_Hy[2] - 2;
                index = i + grid.size_Hy[0] * ( j + k * grid.size_Hy[1]);
                H_y_tmp[index] = 0;
            }
        }
    
    }else{
        DISPLAY_ERROR_ABORT_CLASS(
            "None of the boolean values is true."
        );
    }
}

/**
 * @brief Apply 1D conditions to the electric field.
 */
void AlgoElectro_NEW::apply_1D_case_on_electric_field(
    bool IS_1D_FACE_EX_Electric_along_Z  ,
    bool IS_1D_FACE_EX_Electric_along_Y ,
    bool IS_1D_FACE_EY_Electric_along_Z ,
    bool IS_1D_FACE_EY_Electric_along_X  ,
    bool IS_1D_FACE_EZ_Electric_along_Y , 
    bool IS_1D_FACE_EZ_Electric_along_X  ,
    bool IS_1D_FACE_Minus_EX_Electric_along_Z ,
    bool IS_1D_FACE_Minus_EX_Electric_along_Y,
    bool IS_1D_FACE_Minus_EY_Electric_along_Z ,
    bool IS_1D_FACE_Minus_EY_Electric_along_X ,
    bool IS_1D_FACE_Minus_EZ_Electric_along_X ,
    bool IS_1D_FACE_Minus_EZ_Electric_along_Y ,
    double current_time,
    size_t currentStep,
    GridCreator_NEW &grid,
    double *E_x_tmp,
    double *E_y_tmp,
    double *E_z_tmp
)
{
    // Check that we only have 1 MPI process:
    if(grid.MPI_communicator.getNumberOfMPIProcesses() != 1){
        DISPLAY_ERROR_ABORT_CLASS(
            "Only one MPI process is supported for 1D cases."
        );
    }
    // Check that only one boolean is true:
    size_t counterBooleanTrue = 0;
    if(       IS_1D_FACE_EX_Electric_along_Y == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EX_Electric_along_Z == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EX_Electric_along_Z == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EX_Electric_along_Y == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EY_Electric_along_X == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EY_Electric_along_Z == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EY_Electric_along_X == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EY_Electric_along_Z == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EZ_Electric_along_X == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_EZ_Electric_along_Y == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EZ_Electric_along_X == true){
        counterBooleanTrue ++;
    }else if( IS_1D_FACE_Minus_EZ_Electric_along_Y == true){
        counterBooleanTrue ++;
    }

    if(counterBooleanTrue != 1){
        DISPLAY_ERROR_ABORT_CLASS(
            "Should have only one true boolean but has %zu.",
            counterBooleanTrue
        );
    }

    /////////////////////////////////
    /// IMPOSED SOURCE PARAMETERS ///
    /////////////////////////////////
    double frequency = grid.input_parser.source.frequency[0];
    double period    = 1 / frequency * 10;
    double MEAN      = 2*period;
    double STD       = period / 4.;
    double gauss     
            = exp( - (current_time-MEAN)*(current_time-MEAN) / (2*STD*STD));
    printf("\t> Using frequency %.9g for 1D source.\n",frequency);

    std::vector<double> test_steady_frequencies = {
        1.0,
        0.7394,
        1.9822,
        0.9882,
        1.9785,
        0.0390,
        2.1542,
        1.1734,
        0.1005
    };
    for(size_t i = 0 ; i < test_steady_frequencies.size() ; i ++){
        test_steady_frequencies[i] = test_steady_frequencies[i]*frequency;
    }
    std::vector<double> test_steady_ampl = {
        2.0,
        3.2481,
        5.7305,
        7.3707,
        7.8721,
        7.8674,
        7.1704,
        6.9256,
        6.4077
    };

    double source_value = 0.0;
    if(grid.input_parser.source_time[0] == "GAUSSIAN"){
        source_value = gauss * sin(2*M_PI*frequency*current_time);
    }else if(grid.input_parser.source_time[0] == "SINE"){
        source_value = sin(2*M_PI*frequency*current_time);
    }else if(grid.input_parser.source_time[0] == "TEST_STEADY_STATE_1D"){
        for(size_t I = 0 ; I < test_steady_frequencies.size() ; I ++){
            source_value += test_steady_ampl[I]*sin(
                                2*M_PI*current_time*test_steady_frequencies[I])
                            * exp(-1E17*current_time*current_time*I);
        }
        printf(">>> test steady state source value is %lf.\n",source_value);
    }else{
        DISPLAY_ERROR_ABORT_CLASS(
            "Unknown source time type %s.",
            grid.input_parser.source_time[0].c_str()
        );
    }

    if(      IS_1D_FACE_Minus_EY_Electric_along_Z){          
        /**
         * What will be imposed:
         *  1) On faces with normal e_x:
         *          - Ex is zero.
         *  3) On face with normal -e_y:
         *          - gaussian modulated sinusoidal source for Ez.
         *          - Ex and Ey are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_Minus_EY.\n");
        size_t i, j, index;
        /** IMPOSING EX ON FACE WITH NORMAL E_X **/
        printf("\t\t> Imposing EX on faces with normal +e_x and -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++ ){
            for(size_t k = 1  ; k < grid.size_Ex[2] - 1; k ++){

                i = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

                i = grid.size_Ex[0] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;
            }
        }

        /** Imposing Ex to zero on face with normal -e_y. **/
        printf("\t\t> Imposing EX on face with normal -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ex[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ex[2] - 1 ; k ++){

                j = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

            }
        }

        /** Imposing Ey to zero on face with normal -e_y. **/
        printf("\t\t> Imposing Ey on face with normal -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ey[2] - 1 ; k ++){

                j = 1;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source **/
        /** We impose on face with normal -e_y. **/
        printf("\t\t> Imposing the source EZ on face with normal -e_y...\n");
        j = 1;
        for(size_t i = 1 ; i < grid.size_Ez[0] -1; i ++){
            for(size_t k = 1 ; k < grid.size_Ez[2] - 1 ; k ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);

                E_z_tmp[index] = source_value;
            }
        }

    }else if(IS_1D_FACE_EY_Electric_along_Z){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_x and -e_x:
         *          - Ex is zero.
         *  3) On face with normal +e_y:
         *          - gaussian modulated sinusoidal source for Ez.
         *          - Ex and Ey are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_EY_Electric_along_Z.\n");
        size_t i, j, index;
        /** IMPOSING EX ON FACE WITH NORMAL E_X **/
        printf("\t\t> Imposing EX on faces with normal +e_x and -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++ ){
            for(size_t k = 1  ; k < grid.size_Ex[2] - 1; k ++){

                i = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

                i = grid.size_Ex[0] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;
            }
        }

        /** Imposing Ex to zero on face with normal +e_y. **/
        printf("\t\t> Imposing EX on face with normal +e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ex[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ex[2] - 1 ; k ++){

                j = grid.size_Ex[1] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

            }
        }

        /** Imposing Ey to zero on face with normal +e_y. **/
        printf("\t\t> Imposing Ey on face with normal +e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ey[2] - 1 ; k ++){

                j = grid.size_Ey[1] - 2;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source **/
        /** We impose on face with normal +e_y. **/
        printf("\t\t> Imposing the source EZ on face with normal +e_y...\n");
        j = grid.size_Ez[1] - 2;
        for(size_t i = 1 ; i < grid.size_Ez[0] -1; i ++){
            for(size_t k = 1 ; k < grid.size_Ez[2] - 1 ; k ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);

                E_z_tmp[index] = source_value;
            }
        }

    }else if(IS_1D_FACE_EY_Electric_along_X){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_x and -e_x:
         *          - Ex is zero.
         *  3) On face with normal +e_y:
         *          - gaussian modulated sinusoidal source for Ez.
         *          - Ex and Ey are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_EY_Electric_along_X.\n");
        size_t i, j, index;
        /** IMPOSING EX ON FACE WITH NORMAL E_X **/
        //printf("\t\t> Imposing EX on faces with normal +e_x and -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++ ){
            for(size_t k = 1  ; k < grid.size_Ex[2] - 1; k ++){

                i = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                //E_x_tmp[index] = 0;

                i = grid.size_Ex[0] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                //E_x_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source **/
        /** We impose on face with normal +e_y. **/
        printf("\t\t> Imposing the source EX on face with normal +e_y...\n");
        j = grid.size_Ex[1] - 2;
        for(size_t i = 1 ; i < grid.size_Ex[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ex[2] - 1 ; k ++){

                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = source_value;

            }
        }

        /** Imposing Ey to zero on face with normal +e_y. **/
        printf("\t\t> Imposing Ey on face with normal +e_y...\n");
        j = grid.size_Ey[1] - 2;
        for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ey[2] - 1 ; k ++){

                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing Ey to zero on face with normal +e_y. **/
        printf("\t\t> Imposing EZ on face with normal +e_y...\n");
        j = grid.size_Ez[1] - 2;
        for(size_t i = 1 ; i < grid.size_Ez[0] -1; i ++){
            for(size_t k = 1 ; k < grid.size_Ez[2] - 1 ; k ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);
                E_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_Minus_EX_Electric_along_Z){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y and -e_y:
         *          - Ey is zero.
         *  3) On face with normal -e_x:
         *          - gaussian modulated sinusoidal source for Ez.
         *          - Ex and Ey are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_Minus_EX.\n");
        size_t i, j, index;
        /** IMPOSING EY ON FACE WITH NORMAL E_Y **/
        printf("\t\t> Imposing EY on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++ ){
            for(size_t k = 1  ; k < grid.size_Ey[2] - 1; k ++){

                j = 1;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;

                j = grid.size_Ey[1] - 2;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing Ex to zero on face with normal -e_x. **/
        printf("\t\t> Imposing EX on face with normal -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Ex[2] - 1 ; k ++){

                i = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

            }
        }

        /** Imposing Ey to zero on face with normal -e_x. **/
        printf("\t\t> Imposing Ey on face with normal -e_x...\n");
        for(size_t j = 1 ; j < grid.size_Ey[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Ey[2] - 1 ; k ++){

                i = 1;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source **/
        /** We impose on face with normal -e_x. **/
        printf("\t\t> Imposing the source EZ on face with normal -e_x...\n");
        i = 1;
        for(size_t j = 1 ; j < grid.size_Ez[1] -1; j ++){
            for(size_t k = 1 ; k < grid.size_Ez[2] - 1 ; k ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);

                E_z_tmp[index] = source_value;
            }
        }

    }else if(IS_1D_FACE_Minus_EX_Electric_along_Y){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y and -e_y:
         *          - Ey is zero.
         *  3) On face with normal -e_x:
         *          - gaussian modulated sinusoidal source for Ey.
         *          - Ex and Ez are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_Minus_EX_Electric_along_Y.\n");
        size_t i, j, index;
        /** IMPOSING EY ON FACE WITH NORMAL E_Y **/
        printf("\t\t> Imposing EY on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++ ){
            for(size_t k = 1  ; k < grid.size_Ey[2] - 1; k ++){

                j = 1;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                //E_y_tmp[index] = 0;

                j = grid.size_Ey[1] - 2;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                //E_y_tmp[index] = 0;
            }
        }

        /** Imposing Ex to zero on face with normal -e_x. **/
        printf("\t\t> Imposing EX on face with normal -e_x...\n");
        i = 1;
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Ex[2] - 1 ; k ++){

                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

            }
        }

        /** Imposing gaussian modulated source **/
        /** We impose on face with normal -e_x. **/
        printf("\t\t> Imposing the source EY on face with normal -e_x...\n");
        i = 1;
        for(size_t j = 1 ; j < grid.size_Ey[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Ey[2] - 1 ; k ++){

                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = source_value;
            }
        }
        printf(">>> Imposed is %.9g.\n",
            gauss * sin(2*M_PI*frequency*current_time));

        /** Imposing Ey to zero on face with normal -e_x. **/
        printf("\t\t> Imposing Ez on face with normal -e_x...\n");
        i = 1;
        for(size_t j = 1 ; j < grid.size_Ez[1] -1; j ++){
            for(size_t k = 1 ; k < grid.size_Ez[2] - 1 ; k ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);
                E_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EX_Electric_along_Z){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y and -e_y:
         *          - Ey is zero.
         *  3) On face with normal +e_x:
         *          - gaussian modulated sinusoidal source for Ez.
         *          - Ex and Ey are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_EX_Electric_along_Z.\n");
        size_t i, j, index;
        /** IMPOSING EY ON FACE WITH NORMAL E_Y **/
        printf("\t\t> Imposing EY on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++ ){
            for(size_t k = 1  ; k < grid.size_Ey[2] - 1; k ++){

                j = 1;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;

                j = grid.size_Ey[1] - 2;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing Ex to zero on face with normal e_x. **/
        printf("\t\t> Imposing EX on face with normal e_x...\n");
        i = grid.size_Ex[0] - 2;
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Ex[2] - 1 ; k ++){

                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

            }
        }

        /** Imposing Ey to zero on face with normal -e_x. **/
        printf("\t\t> Imposing Ey on face with normal -e_x...\n");
        i = grid.size_Ey[0] - 2;
        for(size_t j = 1 ; j < grid.size_Ey[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Ey[2] - 1 ; k ++){

                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source **/
        /** We impose on face with normal e_x. **/
        printf("\t\t> Imposing the source EZ on face with normal e_x...\n");
        i = grid.size_Ez[0]-2;
        for(size_t j = 1 ; j < grid.size_Ez[1] -1; j ++){
            for(size_t k = 1 ; k < grid.size_Ez[2] - 1 ; k ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);
                E_z_tmp[index] = source_value;
            }
        }

    }else if(IS_1D_FACE_Minus_EZ_Electric_along_Y){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y and -e_y:
         *          - EX is zero.
         *  3) On face with normal -e_z:
         *          - gaussian modulated sinusoidal source for Ey.
         *          - Ex and Ez are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_Minus_EZ.\n");
        size_t i, k, index;
        /** IMPOSING EX ON FACE WITH NORMAL E_Y **/
        printf("\t\t> Imposing EX on faces with normal +e_y and -e_y...\n");
        for(size_t j = 1 ; j < grid.size_Ex[0] - 1 ; j ++ ){
            for(size_t k = 1  ; k < grid.size_Ex[2] - 1; k ++){

                i = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

                i = grid.size_Ex[0] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;
            }
        }

        /** Imposing Ex to zero on face with normal -e_z. **/
        printf("\t\t> Imposing EX on face with normal -e_z...\n");
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++){
            for(size_t i = 1 ; i < grid.size_Ex[0] - 1 ; i ++){

                k = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source on Ey **/
        /** We impose on face with normal -e_z. **/
        printf("\t\t> Imposing the source EY on face with normal -e_z...\n");
        for(size_t j = 1 ; j < grid.size_Ey[1] - 1 ; j ++){
            for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++){

                k = 1;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = source_value;
            }
        }

        /** Imposing EZ to zero on face with normal -e_z. **/
        printf("\t\t> Imposing Ez on face with normal -e_z...\n");
        k = 1;
        for(size_t j = 1 ; j < grid.size_Ez[1] -1; j ++){
            for(size_t i = 1 ; i < grid.size_Ez[0] - 1 ; i ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);

                E_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EZ_Electric_along_Y){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y and -e_y:
         *          - EX is zero.
         *  3) On face with normal +e_z:
         *          - gaussian modulated sinusoidal source for Ey.
         *          - Ex and Ez are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_EZ.\n");
        size_t i, k, index;
        /** IMPOSING EX ON FACE WITH NORMAL E_Y **/
        printf("\t\t> Imposing EX on faces with normal +e_y and -e_y...\n");
        for(size_t j = 1 ; j < grid.size_Ex[0] - 1 ; j ++ ){
            for(size_t k = 1  ; k < grid.size_Ex[2] - 1; k ++){

                i = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

                i = grid.size_Ex[0] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;
            }
        }

        /** Imposing Ex to zero on face with normal +e_z. **/
        printf("\t\t> Imposing EX on face with normal +e_z...\n");
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++){
            for(size_t i = 1 ; i < grid.size_Ex[0] - 1 ; i ++){

                k = grid.size_Ex[2] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source on Ey**/
        /** We impose on face with normal +e_z. **/
        printf("\t\t> Imposing the source EY on face with normal +e_z...\n");
        k = grid.size_Ey[2] - 2;
        for(size_t j = 1 ; j < grid.size_Ey[1] - 1 ; j ++){
            for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++){

                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = source_value;
            }
        }

        /** Imposing Ey to zero on face with normal +e_x. **/
        printf("\t\t> Imposing Ez to zero on face with normal +e_z...\n");
        k = grid.size_Ez[2] - 2;
        for(size_t j = 1 ; j < grid.size_Ez[1] -1; j ++){
            for(size_t i = 1 ; i < grid.size_Ez[0] - 1 ; i ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);

                E_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EZ_Electric_along_X){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y and -e_y:
         *          - EX is zero.
         *  3) On face with normal +e_z:
         *          - gaussian modulated sinusoidal source for Ey.
         *          - Ex and Ez are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_EZ_Electric_along_X.\n");
        size_t i, k, index;
        /** IMPOSING EX ON FACE WITH NORMAL E_Y **/
        printf("\t\t> Imposing EX on faces with normal +e_y and -e_y...\n");
        for(size_t j = 1 ; j < grid.size_Ex[0] - 1 ; j ++ ){
            for(size_t k = 1  ; k < grid.size_Ex[2] - 1; k ++){

                i = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                //E_x_tmp[index] = 0;

                i = grid.size_Ex[0] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                //E_x_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source on Ey**/
        /** We impose on face with normal +e_z. **/
        printf("\t\t> Imposing the source EX on face with normal +e_z...\n");
        k = grid.size_Ex[2] - 2;
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++){
            for(size_t i = 1 ; i < grid.size_Ex[0] - 1 ; i ++){

                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = source_value;
            }
        }

        /** Imposing Ey to zero on face with normal +e_z. **/
        printf("\t\t> Imposing EY on face with normal +e_z...\n");
        k = grid.size_Ey[2] - 2;
        for(size_t j = 1 ; j < grid.size_Ey[1] - 1 ; j ++){
            for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++){

                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing Ey to zero on face with normal +e_x. **/
        printf("\t\t> Imposing Ez to zero on face with normal +e_z...\n");
        k = grid.size_Ez[2] - 2;
        for(size_t j = 1 ; j < grid.size_Ez[1] -1; j ++){
            for(size_t i = 1 ; i < grid.size_Ez[0] - 1 ; i ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);

                E_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_Minus_EZ_Electric_along_X){
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y and -e_y:
         *          - EX is zero.
         *  3) On face with normal +e_z:
         *          - gaussian modulated sinusoidal source for Ey.
         *          - Ex and Ez are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_Minus_EZ_Electric_along_X.\n");
        size_t i, k, index;
        /** IMPOSING EX ON FACE WITH NORMAL E_Y **/
        printf("\t\t> Imposing EX on faces with normal +e_y and -e_y...\n");
        for(size_t j = 1 ; j < grid.size_Ex[0] - 1 ; j ++ ){
            for(size_t k = 1  ; k < grid.size_Ex[2] - 1; k ++){

                i = 1;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                //E_x_tmp[index] = 0;

                i = grid.size_Ex[0] - 2;
                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                //E_x_tmp[index] = 0;
            }
        }

        /** Imposing gaussian modulated source on Ex**/
        /** We impose on face with normal -e_z. **/
        printf("\t\t> Imposing the source EX on face with normal -e_z...\n");
        k = 1;
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++){
            for(size_t i = 1 ; i < grid.size_Ex[0] - 1 ; i ++){

                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = source_value;
            }
        }

        /** Imposing Ey to zero on face with normal -e_z. **/
        printf("\t\t> Imposing EY on face with normal -e_z...\n");
        k = 1;
        for(size_t j = 1 ; j < grid.size_Ey[1] - 1 ; j ++){
            for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++){

                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing Ez to zero on face with normal -e_z. **/
        printf("\t\t> Imposing Ez to zero on face with normal -e_z...\n");
        k = 1;
        for(size_t j = 1 ; j < grid.size_Ez[1] -1; j ++){
            for(size_t i = 1 ; i < grid.size_Ez[0] - 1 ; i ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);

                E_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_Minus_EY_Electric_along_X){
        /**
         * Propagation of Ex along +y direction.
         * Propagation of Hz along +y direction.
         * What will be imposed:
         *      1) On face with normal -e_y:
         *              - Ex is a gaussian modulated sine wave.
         *              - Ey and Ez are imposed to zero.
         *      2) On faces with normal +e_z and -e_z:
         *              - Ez is imposed to zero.
         */
        size_t j, k, index;
        printf("\t> > 1D case is IS_1D_FACE_Minus_EY_Electric_along_X.\n");

        /** Imposing Ez on faces with normal +e_z and -e_z. **/
        printf("\t\t> Imposing EZ to zero on faces with normal +e_z and -e_z...\n");
        for(size_t i = 1 ; i < grid.size_Ez[0] - 1 ; i ++){
            for(size_t j = 1 ; j < grid.size_Ez[1] - 1 ; j ++){

                k = 1;
                index = i + grid.size_Ez[0] * ( j + k * grid.size_Ez[1]);
                E_z_tmp[index] = 0;

                k = grid.size_Ez[2] - 2;
                index = i + grid.size_Ez[0] * ( j + k * grid.size_Ez[1]);
                E_z_tmp[index] = 0;
            }
        }

        /** Imposing Ex to a gaussian modulated sine wave on face with normal -e_y. **/
        j = 1;
        printf("\t\t> Imposing the source EX on face with normal -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ex[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ex[2] - 1 ; k ++){

                index = i + grid.size_Ex[0] * ( j + k * grid.size_Ex[1]);
                E_x_tmp[index] = source_value;
            }
        }

        /** Imposing Ey to zero on face with normal -e_y. **/
        printf("\t\t> Imposing EY to zero on face with normal -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ey[2] - 1 ; k ++){

                index = i + grid.size_Ey[0] * ( j + k * grid.size_Ey[1]);
                E_y_tmp[index] = 0;
            }
        }

        /** Imposing Ez to zero on face with normal -e_y. **/
        printf("\t\t> Imposing EZ to zero on face with normal -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ez[0] - 1 ; i ++){
            for(size_t k = 1 ; k < grid.size_Ez[2] - 1 ; k ++){

                index = i + grid.size_Ez[0] * ( j + k * grid.size_Ez[1]);
                E_z_tmp[index] = 0;
            }
        }

    }else if(IS_1D_FACE_EX_Electric_along_Y){
        // Propagation of electric field Ey along -x direction.
        // Propagation of magnetic field Hz along -x direction.
        /**
         * What will be imposed:
         *  1) On faces with normal +e_y and -e_y:
         *          - Ey is zero.
         *  3) On face with normal +e_x:
         *          - gaussian modulated sinusoidal source for Ey.
         *          - Ex and Ez are imposed to zero.
         */
        printf("\t> > 1D case is IS_1D_FACE_EX_Electric_along_Y.\n");
        size_t i, j, index;
        /** IMPOSING EY ON FACE WITH NORMAL E_Y **/
        printf("\t\t> Imposing EY on faces with normal +e_y and -e_y...\n");
        for(size_t i = 1 ; i < grid.size_Ey[0] - 1 ; i ++ ){
            for(size_t k = 1  ; k < grid.size_Ey[2] - 1; k ++){

                j = 1;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                //E_z_tmp[index] = 0;

                j = grid.size_Ey[1] - 2;
                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                //E_z_tmp[index] = 0;
            }
        }

        /** Imposing Ex to zero on face with normal e_x. **/
        printf("\t\t> Imposing EX on face with normal e_x...\n");
        i = grid.size_Ex[0] - 2;
        for(size_t j = 1 ; j < grid.size_Ex[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Ex[2] - 1 ; k ++){

                index = i + grid.size_Ex[0] * ( j + grid.size_Ex[1] * k);
                E_x_tmp[index] = 0;

            }
        }

        /** Imposing gaussian modulated source **/
        /** We impose on face with normal e_x. **/
        printf("\t\t> Imposing the source EY on face with normal e_x...\n");
        i = grid.size_Ey[0] - 2;
        for(size_t j = 1 ; j < grid.size_Ey[1] - 1 ; j ++){
            for(size_t k = 1 ; k < grid.size_Ey[2] - 1 ; k ++){

                index = i + grid.size_Ey[0] * ( j + grid.size_Ey[1] * k);
                E_y_tmp[index] = source_value;
            }
        }

        /** Imposing Ez to zero on face with normal -e_x. **/
        printf("\t\t> Imposing Ez on face with normal -e_x...\n");
        i = grid.size_Ez[0] - 2;
        for(size_t j = 1 ; j < grid.size_Ez[1] -1; j ++){
            for(size_t k = 1 ; k < grid.size_Ez[2] - 1 ; k ++){

                index = i + grid.size_Ez[0] * ( j + grid.size_Ez[1] * k);

                E_z_tmp[index] = 0;
            }
        }

    }else{
        DISPLAY_ERROR_ABORT_CLASS(
            "None of the boolean values is true."
        );
    }
}
