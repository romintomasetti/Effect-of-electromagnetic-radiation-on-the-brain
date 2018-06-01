#ifndef ALGOELECTRO_NEW_HPP
#define ALGOELECTRO_NEW_HPP

#include "GridCreator_NEW.h"

#include <float.h>
#include <limits.h>
#include <cmath>


#include "InterfaceToParaviewer.h"

class AlgoElectro_NEW{
    public:
        /* MEMBERS */
		unsigned int VERBOSITY = 0;
        
        /* FUNCTIONS */

        // Compute the time step required to fulfill the theoretical stability condition:
        double Compute_dt(GridCreator_NEW & /*grid*/);


        /* CONSTRUCTOR */
        AlgoElectro_NEW(unsigned int VERBOSITY){
			this->VERBOSITY = VERBOSITY;
		}

        /* DESTRUCTOR  */
        ~AlgoElectro_NEW(void){}

        /* Update function */
        void update(GridCreator_NEW &,InterfaceToParaviewer &);


        // Check that OMP_DYNAMIC is set to false:
        void check_OMP_DYNAMIC_envVar(void);


        /* Update the points with boundary conditions  */
        void abc( GridCreator_NEW &grid, 
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
        );
        

        void determine_OMP_thread_role_in_MPI_communication(
            int omp_thread_id /* omp_get_thread_num()*/,
            bool *OMP_thread_has_neighboor,
            char *direction,
            GridCreator_NEW &grid,
            std::vector<size_t> &,
            double **,
            double **,
            size_t * /*size_of_sent_vector*/,
            size_t * /*size_of_recv_vector*/
        );

        bool SteadyStateAnalyser(
            const bool is_steady_state_for_this_MPI,
            bool *is_steady_state_for_all_mpi,
            GridCreator_NEW &grid,
            const double time_beg,
            const double time_end,
            const std::vector<double> &Ez_trapz_absolute,
            const std::vector<double> &Ey_trapz_absolute,
            const std::vector<double> &Ex_trapz_absolute,
            double dt
        );

        /**
         * @brief Apply 1D conditions on the magnetic field.
         */
        void apply_1D_case_on_magnetic_field(
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
        );

        /**
         * @brief Apply 1D conditions on the electric field.
         */
        void apply_1D_case_on_electric_field(
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
        );
        
        void pmlE( GridCreator_NEW &grid,
            double *Ex, double *Ey, double *Ez,
            double *Ex_pml_x0, double *Ex_pml_x1,
            double *Ex_pml_y0, double *Ex_pml_y1, 
            double *Ex_pml_z0, double *Ex_pml_z1,
            double *Ey_pml_x0, double *Ey_pml_x1,
            double *Ey_pml_y0, double *Ey_pml_y1,
            double *Ey_pml_z0, double *Ey_pml_z1,
            double *Ez_pml_x0, double *Ez_pml_x1, 
            double *Ez_pml_y0, double *Ez_pml_y1,
            double *Ez_pml_z0, double *Ez_pml_z1,

            double *Hx, double *Hy, double *Hz,
            size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ, size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ,
            size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY, size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY,
            size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX, size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX,
            
            double *C_exe, double *C_exh_1, double *C_exh_2, double *C_exe2,
            double *C_eye, double *C_eyh_1, double *C_eyh_2, double *C_eye2,
            double *C_eze, double *C_ezh_1, double *C_ezh_2, double *C_eze2,
            unsigned int rhoX0, unsigned int rhoX1,
            unsigned int rhoY0, unsigned int rhoY1,
            unsigned int rhoZ0, unsigned int rhoZ1
            );
            
        void pmlH(  GridCreator_NEW &grid,
            double *Hx, double *Hy, double *Hz,
            double *Hx_pml_x0, double *Hx_pml_x1,
            double *Hx_pml_y0, double *Hx_pml_y1, 
            double *Hx_pml_z0, double *Hx_pml_z1,
            double *Hy_pml_x0, double *Hy_pml_x1,
            double *Hy_pml_y0, double *Hy_pml_y1,
            double *Hy_pml_z0, double *Hy_pml_z1,
            double *Hz_pml_x0, double *Hz_pml_x1,
            double *Hz_pml_y0, double *Hz_pml_y1,
            double *Hz_pml_z0, double *Hz_pml_z1,

            double *Ex, double *Ey, double *Ez,

            double *C_hxh, double *C_hxe_1, double *C_hxe_2, double *C_hxh2,
            double *C_hyh, double *C_hye_1, double *C_hye_2, double *C_hyh2,
            double *C_hzh, double *C_hze_1, double *C_hze_2, double *C_hzh2,
            unsigned int rhoX0, unsigned int rhoX1,
            unsigned int rhoY0, unsigned int rhoY1,
            unsigned int rhoZ0, unsigned int rhoZ1
            );

    void WriteData(int MPI_my_rank, GridCreator_NEW &grid);

    void ComputePowerInBrain(unsigned int *GlobalIndexVECTOR,
                            unsigned int *LocalIndexVector, 
                            unsigned int counterInBrain,
                            GridCreator_NEW &grid,
                            double *PowerInBrain);

    void WritePowerInVacuum(double current_time,
                            GridCreator_NEW &grid,
                            size_t rhoX0, size_t rhoX1,
                            size_t rhoY0, size_t rhoY1,
                            size_t rhoZ0, size_t rhoZ1,
                            size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX,
                            size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY,
                            size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ,
                            size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX,
                            size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY,
                            size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ);

    double ComputePowerInVacuum( GridCreator_NEW &grid,
                                size_t rhoX0, size_t rhoX1,
                                size_t rhoY0, size_t rhoY1,
                                size_t rhoZ0, size_t rhoZ1,
                                size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDX,
                                size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDY,
                                size_t IS_THE_FIRST_MPI_FOR_ELECRIC_FIELDZ,
                                size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDX,
                                size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDY,
                                size_t IS_THE_LAST_MPI_FOR_ELECTRIC_FIELDZ);
};

#endif
