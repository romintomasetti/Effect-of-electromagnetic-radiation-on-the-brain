#ifndef INTERFACETOPARAVIEWER_H
#define INTERFACETOPARAVIEWER_H

#include "GridCreator.h"

#include "GridCreator_NEW.h"

#include "vtl.h"
#include "Array_3D_Template.h"

#include <cstring>

class InterfaceToParaviewer{
    private:
        // Reference to a GridCreator object:
        GridCreator     &grid_Creator;
        // Reference to a GridCreator_NEW object:
        GridCreator_NEW &grid_Creator_NEW;
        // Reference to a MPI_Initializer object:
        MPI_Initializer &MPI_communicator;
        
        // List of subgrids to pass at pvti file encoder, for electromagnetic grid:
        std::vector<vtl::SPoints> sgrids_Electro;

        // List of subgrids to pass at pvti file encoder, for thermal grid:
        std::vector<vtl::SPoints> sgrids_Thermal;

        // Whole electromagnetic grid:
        vtl::SPoints grid_Electro;

        // Whole thermal grid:
        vtl::SPoints grid_Thermal;

        // My electromagnetic grid:
        vtl::SPoints mygrid_Electro;

        // My thermal grid:
        vtl::SPoints mygrid_Thermal;

        // Decide if it is GidCreator or GridCreator_NEW:
        bool is_grid_creator_new;

    public:
        // Default constructor:
        InterfaceToParaviewer(GridCreator   &grid_Creator,
                            MPI_Initializer &MPI_communicator,
                            GridCreator_NEW &grid_Creator_NEW,
                            bool is_grid_creator_new):
                            grid_Creator(grid_Creator),
                            MPI_communicator(MPI_communicator),
                            grid_Creator_NEW(grid_Creator_NEW){
                                this->initializeAll();
                                printf("InterfaceToParaviewer::constructor::OUT\n");
                                this->is_grid_creator_new = is_grid_creator_new;
                            };
        // Default destructor:
        ~InterfaceToParaviewer(void){printf("InterfaceToParaviewer::destructor::out\n");};

        // Initilize everything:
        void initializeAll(void);

        // Convert and write output:
        void convertAndWriteData(unsigned long currentStep,
                std::string type /*"thermal" or "electro", case sensitive*/);
};

#endif