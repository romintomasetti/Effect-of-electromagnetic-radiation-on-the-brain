#ifndef INTERFACETOPARAVIEWER_H
#define INTERFACETOPARAVIEWER_H


#include "GridCreator_NEW.h"

#include "vtl_romin.h"
#include "Array_3D_Template.h"

#include <cstring>

class InterfaceToParaviewer{
    private:

        // Reference to a GridCreator_NEW object:
        GridCreator_NEW &grid_Creator_NEW;
        // Reference to a MPI_Initializer object:
        MPI_Initializer &MPI_communicator;
        
        // List of subgrids to pass at pvti file encoder, for electromagnetic grid:
        std::vector<vtl_romin::SPoints> sgrids_Electro;

        // List of subgrids to pass at pvti file encoder, for thermal grid:
        std::vector<vtl_romin::SPoints> sgrids_Thermal;

        // Whole electromagnetic grid:
        vtl_romin::SPoints grid_Electro;

        // Whole thermal grid:
        vtl_romin::SPoints grid_Thermal;

        // My electromagnetic grid:
        vtl_romin::SPoints mygrid_Electro;

        // My thermal grid:
        vtl_romin::SPoints mygrid_Thermal;

    public:
        // Default constructor:
        InterfaceToParaviewer(MPI_Initializer &MPI_communicator,
                            GridCreator_NEW &grid_Creator_NEW):
                            grid_Creator_NEW(grid_Creator_NEW),
                            MPI_communicator(MPI_communicator)
                            {
                                this->initializeAll();
                            };
        // Default destructor:
        ~InterfaceToParaviewer(void){};

        // Initilize everything:
        void initializeAll(void);

        // Convert and write output:
        void convertAndWriteData(unsigned long currentStep,
                std::string type /*"thermal" or "electro", case sensitive*/);

        std::string create_folder_and_go_in(std::string folderName);
};

#endif