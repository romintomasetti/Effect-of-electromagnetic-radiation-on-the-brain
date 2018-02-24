#ifndef INTERFACETOPARAVIEWER_H
#define INTERFACETOPARAVIEWER_H

#include "GridCreator.h"
#include "vtl/vtl.h"
#include "Array_3D_Template.h"

class InterfaceToParaviewer{
    private:
        // Reference to a GridCreator object:
        GridCreator     &grid_Creator;
        // Reference to a MPI_Initializer object:
        MPI_Initializer &MPI_communicator;
        
        // List of subgrids to pass at pvti file encoder:
        std::vector<vtl::SPoints> sgrids;

        // Whole grid:
        vtl::SPoints grid;

        // My grid:
        vtl::SPoints mygrid;

    public:
        // Default constructor:
        InterfaceToParaviewer(GridCreator   &grid_Creator,
                            MPI_Initializer &MPI_communicator):
                            grid_Creator(grid_Creator),
                            MPI_communicator(MPI_communicator){
                                this->initializeAll();
                                printf("InterfaceToParaviewer::constructor::OUT\n");
                            };
        // Default destructor:
        ~InterfaceToParaviewer(void){printf("InterfaceToParaviewer::destructor::out\n");};

        // Initilize everything:
        void initializeAll(void);

        // Convert and write output:
        void convertAndWriteData(unsigned long currentStep);
};

#endif