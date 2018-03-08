#ifndef ALGOELECTRO_NEW_HPP
#define ALGOELECTRO_NEW_HPP

#include "GridCreator_NEW.h"
class AlgoElectro_NEW{
    private:
        /* MEMBERS */
        
        /* FUNCTIONS */

        // Compute the time step required to fulfill the theoretical stability condition:
        double Compute_dt(GridCreator_NEW & /*grid*/);

    public:

        /* CONSTRUCTOR */
        AlgoElectro_NEW(void){}

        /* DESTRUCTOR  */
        ~AlgoElectro_NEW(void){}

        /* Update function */
        void update(GridCreator_NEW &);

        // Check that OMP_DYNAMIC is set to false:
        void check_OMP_DYNAMIC_envVar(void);
};

#endif
