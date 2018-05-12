#ifndef ALGOELECTRO_NEW_HPP
#define ALGOELECTRO_NEW_HPP

#include "GridCreator_NEW.h"

#include <float.h>
#include <limits.h>


#include "InterfaceToParaviewer.h"

class AlgoElectro_NEW{
    private:
        /* MEMBERS */
		unsigned int VERBOSITY = 0;
        
        /* FUNCTIONS */

        // Compute the time step required to fulfill the theoretical stability condition:
        double Compute_dt(GridCreator_NEW & /*grid*/);

    public:

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

        std::vector<double> ComputeNormE2square(std::vector<double> Ex,
                                          std::vector<double> Ey, std::vector<double> Ez);


        double interpolationX(size_t x1, size_t x2, size_t x3,
                              size_t x4, size_t x5, size_t x6,
                              size_t x7, size_t x8,  GridCreator_NEW &grid);

        double interpolationX2(size_t x1, size_t x2, size_t x3,
                              size_t x4, size_t x5, size_t x6,
                              size_t x7, size_t x8,  std::vector<double> Ex);

        
        double interpolationY(size_t y1, size_t y2, size_t y3,
                              size_t y4, size_t y5, size_t y6,
                              size_t y7, size_t y8,  GridCreator_NEW &grid);

        double interpolationY2(size_t y1, size_t y2, size_t y3,
                              size_t y4, size_t y5, size_t y6,
                              size_t y7, size_t y8,  std::vector<double> Ey);


        double interpolationZ(size_t z1, size_t z2, size_t z3,
                              size_t z4, size_t z5, size_t z6,
                              size_t z7, size_t z8,  GridCreator_NEW &grid);

        double interpolationZ2(size_t z1, size_t z2, size_t z3,
                              size_t z4, size_t z5, size_t z6,
                              size_t z7, size_t z8,  std::vector<double> Ez);


        std::vector<double> ComputeNormEsquareBIS(GridCreator_NEW &grid);

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

};

#endif
