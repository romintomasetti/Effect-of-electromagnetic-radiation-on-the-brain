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

        // std::vector<double> ComputeNormE2square(std::vector<double> Ex,
        //                                   std::vector<double> Ey, std::vector<double> Ez, GridCreator_NEW &grid);
        std::vector<double> ComputeNormE2square(GridCreator_NEW &grid);
        

        double interpolationX(size_t x1, size_t x2, size_t x3,
                              size_t x4, size_t x5, size_t x6,
                              size_t x7, size_t x8,  GridCreator_NEW &grid);

        // double interpolationX2(size_t x1, size_t x2, size_t x3,
        //                       size_t x4, size_t x5, size_t x6,
        //                       size_t x7, size_t x8,  std::vector<double> Ex);
        double interpolationX2(size_t x1, size_t x2, size_t x3,
                              size_t x4, size_t x5, size_t x6,
                              size_t x7, size_t x8,  double *Ex);

        
        double interpolationY(size_t y1, size_t y2, size_t y3,
                              size_t y4, size_t y5, size_t y6,
                              size_t y7, size_t y8,  GridCreator_NEW &grid);

        // double interpolationY2(size_t y1, size_t y2, size_t y3,
        //                       size_t y4, size_t y5, size_t y6,
        //                       size_t y7, size_t y8,  std::vector<double> Ey);
        double interpolationY2(size_t y1, size_t y2, size_t y3,
                              size_t y4, size_t y5, size_t y6,
                              size_t y7, size_t y8,  double* Ey);


        double interpolationZ(size_t z1, size_t z2, size_t z3,
                              size_t z4, size_t z5, size_t z6,
                              size_t z7, size_t z8,  GridCreator_NEW &grid);

        // double interpolationZ2(size_t z1, size_t z2, size_t z3,
        //                       size_t z4, size_t z5, size_t z6,
        //                       size_t z7, size_t z8,  std::vector<double> Ez);
        double interpolationZ2(size_t z1, size_t z2, size_t z3,
                              size_t z4, size_t z5, size_t z6,
                              size_t z7, size_t z8,  double* Ez);


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


    size_t findMin(size_t a,
                    size_t b,
                    size_t c);

    size_t findMinInVec(std::vector<size_t> vec);

    void WriteData(int MPI_my_rank, GridCreator_NEW &grid);

    size_t findMaxVectorX(std::vector<size_t> GlobalIndexVECTOR);

    size_t findMinVectorX(std::vector<size_t> GlobalIndexVECTOR);

    size_t findMaxVectorY(std::vector<size_t> GlobalIndexVECTOR);

    size_t findMinVectorY(std::vector<size_t> GlobalIndexVECTOR);

    size_t findMaxVectorZ(std::vector<size_t> GlobalIndexVECTOR);

    size_t findMinVectorZ(std::vector<size_t> GlobalIndexVECTOR);

    void ComputePowerInBrain(unsigned int *GlobalIndexVECTOR,
                            unsigned int *LocalIndexVector, 
                            unsigned int counterInBrain,
                            GridCreator_NEW &grid,
                            double *PowerInBrain);

};

#endif
