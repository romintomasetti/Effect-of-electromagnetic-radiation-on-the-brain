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
                bool   IS_1D_FACE_EX
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

        bool SteadyStateAnalyser(void);

};

#endif
