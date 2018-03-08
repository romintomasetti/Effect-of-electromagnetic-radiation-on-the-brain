#include "InterfaceToParaviewer.h"
#include "vtlSPoints.h"

#include "unistd.h"

#include <string>
#include <string.h>
#include <cstring>
#include <map>

#define TAG_NP1_ELECTRO 20
#define TAG_NP2_ELECTRO 21
#define TAG_NP1_THERMAL 22
#define TAG_NP2_THERMAL 23

/**
 * 
 * @brief Initialize all the grid fields via MPI communication.
 * 
 * Detailed explanation goes here.
 */
void InterfaceToParaviewer::initializeAll(void){
    /* MAIN GOAL : create a vector of subgrids */

    // Retrieve the number of MPI processes:
    int nb_MPI = this->MPI_communicator.getNumberOfMPIProcesses();

    double dx_Electro = 0.0;
    double dy_Electro = 0.0;
    double dz_Electro = 0.0;
    double delta_thermal = 0.0;

    if(this->is_grid_creator_new == false){
        dx_Electro = this->grid_Creator.input_parser.deltaX_Electro;
        dy_Electro = this->grid_Creator.input_parser.deltaY_Electro;
        dz_Electro = this->grid_Creator.input_parser.deltaZ_Electro;

        delta_thermal = this->grid_Creator.input_parser.delta_Thermal;
    }else{
        dx_Electro = this->grid_Creator_NEW.input_parser.deltaX_Electro;
        dy_Electro = this->grid_Creator_NEW.input_parser.deltaY_Electro;
        dz_Electro = this->grid_Creator_NEW.input_parser.deltaZ_Electro;

        delta_thermal = this->grid_Creator_NEW.input_parser.delta_Thermal;
    }

    /* CHECK DELTAS */
    if(dx_Electro <= 0.0 || dy_Electro <= 0.0 || dz_Electro <= 0.0 || delta_thermal <= 0.0){
        fprintf(stderr,"InterfaceToParaviewer::initializeAll::ERROR :: One of the spatial step for the electromagnetic/thermal mesh is negative or zeros.\n");
        fprintf(stderr,"InterfaceToParaviewer::initializeAll::ABORTING (line %d, file %s).\n",__LINE__,__FILE__);
        abort();
    }

    /* INITIALIZE 'grid_thermal'/'grid_Electro' for each MPI process (even for non-root ones) */
    if(this->is_grid_creator_new == false){
        this->grid_Electro.o = this->grid_Creator.originOfWholeSimulation;
        this->grid_Thermal.o = this->grid_Creator.originOfWholeSimulation;
    }else{
        this->grid_Electro.o = this->grid_Creator_NEW.originOfWholeSimulation_Electro;
        this->grid_Thermal.o = this->grid_Creator_NEW.originOfWholeSimulation_Thermal;
    }

    this->grid_Electro.dx = vtl::Vec3d(dx_Electro,dy_Electro,dz_Electro);
    this->grid_Thermal.dx = vtl::Vec3d(delta_thermal,delta_thermal,delta_thermal);
    
    this->grid_Electro.np1 = vtl::Vec3i(0,0,0); 
    this->grid_Thermal.np1 = vtl::Vec3i(0,0,0);

    double length_X_whole_dom_electro = 0.0;
    double length_Y_whole_dom_electro = 0.0;
    double length_Z_whole_dom_electro = 0.0;
    double length_x_whole_dom_thermal = 0.0;
    double length_y_whole_dom_thermal = 0.0;
    double length_z_whole_dom_thermal = 0.0;

    if(this->is_grid_creator_new == false){
        length_X_whole_dom_electro = this->grid_Creator.input_parser.lengthX_WholeDomain_Electro;
        length_Y_whole_dom_electro = this->grid_Creator.input_parser.lengthY_WholeDomain_Electro;
        length_Z_whole_dom_electro = this->grid_Creator.input_parser.lengthZ_WholeDomain_Electro;
        length_x_whole_dom_thermal = this->grid_Creator.input_parser.lengthX_WholeDomain_Thermal;
        length_y_whole_dom_thermal = this->grid_Creator.input_parser.lengthY_WholeDomain_Thermal;
        length_z_whole_dom_thermal = this->grid_Creator.input_parser.lengthZ_WholeDomain_Thermal;
    }else{
        length_X_whole_dom_electro = this->grid_Creator_NEW.input_parser.lengthX_WholeDomain_Electro;
        length_Y_whole_dom_electro = this->grid_Creator_NEW.input_parser.lengthY_WholeDomain_Electro;
        length_Z_whole_dom_electro = this->grid_Creator_NEW.input_parser.lengthZ_WholeDomain_Electro;
        length_x_whole_dom_thermal = this->grid_Creator_NEW.input_parser.lengthX_WholeDomain_Thermal;
        length_y_whole_dom_thermal = this->grid_Creator_NEW.input_parser.lengthY_WholeDomain_Thermal;
        length_z_whole_dom_thermal = this->grid_Creator_NEW.input_parser.lengthZ_WholeDomain_Thermal;
    }

    size_t nodesWholeDom_X_Electro = (size_t) length_X_whole_dom_electro / dx_Electro +1;
    size_t nodesWholeDom_Y_Electro = (size_t) length_Y_whole_dom_electro / dy_Electro +1;
    size_t nodesWholeDom_Z_Electro = (size_t) length_Z_whole_dom_electro / dz_Electro +1;

    this->grid_Electro.np2 = vtl::Vec3i(nodesWholeDom_X_Electro,nodesWholeDom_Y_Electro,nodesWholeDom_Z_Electro);

    size_t nodesWholeDom_X_Thermal = (size_t) length_x_whole_dom_thermal / delta_thermal +1;
    size_t nodesWholeDom_Y_Thermal = (size_t) length_y_whole_dom_thermal / delta_thermal +1;
    size_t nodesWholeDom_Z_Thermal = (size_t) length_z_whole_dom_thermal / delta_thermal +1;
    
    this->grid_Thermal.np2 = vtl::Vec3i(nodesWholeDom_X_Thermal,nodesWholeDom_Y_Thermal,nodesWholeDom_Z_Thermal);

    // Initialize the subgrid if I am the root process:
    if(this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess){
        printf("InterfaceToParaviewer::initializeAll\n");
        printf("\t> From MPI %d, I am the root !\n",this->MPI_communicator.getRank());
        printf("\t> Initializing the subgrids (EM and TH) of class-type SPoints with %d processes.\n",
                    nb_MPI);

        // Loop over each grid and get the starting/ending indices:
        for(unsigned char I = 0 ; I < nb_MPI ; I ++){

            printf("\t> MPI %d :: Initializing sgrid_electro[%d] and sgrid_thermal[%d]...\n",
                    this->MPI_communicator.getRank(),I,I);

            // Giving the MPI ID to the element sgrids[I]:
            this->sgrids_Electro[I].id = I;
            this->sgrids_Thermal[I].id = I;

            // Giving the spacing to the element sgrids[I]:
            this->sgrids_Electro[I].dx = vtl::Vec3d(dx_Electro,dy_Electro,dz_Electro);
            this->sgrids_Thermal[I].dx = vtl::Vec3d(delta_thermal,delta_thermal,delta_thermal);

            // Setting the origin:
            if(this->is_grid_creator_new == false){
                this->sgrids_Electro[I].o = this->grid_Creator.originOfWholeSimulation;
                this->sgrids_Thermal[I].o = this->grid_Creator.originOfWholeSimulation;
            }else{
                this->sgrids_Electro[I].o = this->grid_Creator_NEW.originOfWholeSimulation_Electro;
                this->sgrids_Thermal[I].o = this->grid_Creator_NEW.originOfWholeSimulation_Thermal;
            }

            // Setting the origin and end indices of each subgrid:
            if(I == this->MPI_communicator.rootProcess){
                this->mygrid_Electro.id = this->MPI_communicator.getRank();
                this->mygrid_Thermal.id = this->MPI_communicator.getRank();
                // Don't communicate:
                for(int k = 0 ; k < 3 ; k ++){
                    if(this->is_grid_creator_new == false){
                        // Grid creator:
                        this->mygrid_Electro.np1[k] = this->mygrid_Thermal.np1[k] =
                                        this->grid_Creator.originIndices[k];

                        this->mygrid_Electro.np2[k] = this->mygrid_Thermal.np2[k] =
                                        this->mygrid_Electro.np1[k] + this->grid_Creator.numberOfNodesInEachDir[k];

                        this->sgrids_Electro[I].np1[k] = this->sgrids_Thermal[I].np1[k] =
                                        this->grid_Creator.originIndices[k];

                        this->sgrids_Electro[I].np2[k] = this->sgrids_Thermal[I].np1[k] = 
                                        this->sgrids_Electro[I].np1[k] + this->grid_Creator.numberOfNodesInEachDir[k];
                    }else{
                        // Grid creator new:
                        this->mygrid_Electro.np1[k] = this->grid_Creator_NEW.originIndices_Electro[k];

                        this->mygrid_Thermal.np1[k] = this->grid_Creator_NEW.originIndices_Thermal[k];

                        this->mygrid_Electro.np2[k] = this->grid_Creator_NEW.sizes_EH[k];

                        this->mygrid_Thermal.np2[k] = this->grid_Creator_NEW.size_Thermal[k];

                        this->sgrids_Electro[I].np1[k] = this->mygrid_Electro.np1[k];

                        this->sgrids_Thermal[I].np1[k] = this->mygrid_Thermal.np1[k];

                        this->sgrids_Electro[I].np2[k] = this->mygrid_Electro.np1[k]
                                + this->mygrid_Electro.np2[k];

                        this->sgrids_Thermal[I].np2[k] = this->mygrid_Thermal.np1[k]
                                + this->mygrid_Electro.np2[k];
                    }
                }
            }else{
                unsigned long np1_Electro[3];
                unsigned long np1_Thermal[3];
                unsigned long np2_Electro[3];
                unsigned long np2_Thermal[3];

                // Receive np1_Electro:
                MPI_Recv(&np1_Electro,3,MPI_UNSIGNED_LONG,I,TAG_NP1_ELECTRO,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                printf("\t> From root, received np1_electro from %d.\n",I);

                // Receive np2_Electro:
                MPI_Recv(&np2_Electro,3,MPI_UNSIGNED_LONG,I,TAG_NP2_ELECTRO,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                printf("\t> From root, received np2_electro from %d.\n",I);

                // Receive np1_Thermal:
                MPI_Recv(&np1_Thermal,3,MPI_UNSIGNED_LONG,I,TAG_NP1_THERMAL,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                printf("\t From root, received np1_thermal from %d.\n",I);

                // Receive np2_Thermal:
                MPI_Recv(&np2_Thermal,3,MPI_UNSIGNED_LONG,I,TAG_NP2_THERMAL,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                printf("\t From root, received np2_thermal from %d.\n",I);

                for(int k = 0 ; k < 3 ; k ++){
                    this->sgrids_Electro[I].np1[k] = np1_Electro[k];
                    this->sgrids_Electro[I].np2[k] = np2_Electro[k];

                    this->sgrids_Thermal[I].np1[k] = np1_Thermal[k];
                    this->sgrids_Thermal[I].np2[k] = np2_Thermal[k];
                }
            }

            cout << this->sgrids_Electro[I] << endl;
            cout << this->sgrids_Thermal[I] << endl;
        }

    }else if(this->MPI_communicator.isRootProcess() == INT_MIN){

        printf("InterfaceToParaviewer::initializeAll\n");
        printf("\t> From MPI process %d, I am not the root.\n",this->MPI_communicator.getRank());

        unsigned long np1_Electro[3];
        unsigned long np2_Electro[3];
        unsigned long np1_Thermal[3];
        unsigned long np2_Thermal[3];

        for(int k = 0 ; k  <  3 ; k ++ ){
            if(this->is_grid_creator_new == false){
                np1_Electro[k] = this->grid_Creator.originIndices[k];
                np2_Electro[k] = np1_Electro[k] + this->grid_Creator.numberOfNodesInEachDir[k];
                np1_Thermal[k] = np1_Electro[k];
                np2_Thermal[k] = np2_Electro[k];
            }else{
                
                np1_Electro[k] = this->grid_Creator_NEW.originIndices_Electro[k];
                np2_Electro[k] = np1_Electro[k] + this->grid_Creator_NEW.sizes_EH[k];

                np1_Thermal[k] = this->grid_Creator_NEW.originIndices_Thermal[k];
                np2_Thermal[k] = np1_Thermal[k] + this->grid_Creator_NEW.size_Thermal[k];
            }
        }   

        // Send np1_Electro:
        MPI_Send(&np1_Electro,3,MPI_UNSIGNED_LONG,this->MPI_communicator.rootProcess,TAG_NP1_ELECTRO,MPI_COMM_WORLD);
        printf("\t> From MPI process %d, np1 sent to root.\n",this->MPI_communicator.getRank());

        // Send np2_Electro:
        MPI_Send(&np2_Electro,3,MPI_UNSIGNED_LONG,this->MPI_communicator.rootProcess,TAG_NP2_ELECTRO,MPI_COMM_WORLD);
        printf("\t> From MPI process %d, np2 sent to root.\n",this->MPI_communicator.getRank());

        // Send np1_Thermal:
        MPI_Send(&np1_Thermal,3,MPI_UNSIGNED_LONG,this->MPI_communicator.rootProcess,TAG_NP1_THERMAL,MPI_COMM_WORLD);

        // Send np2_thermal:
        MPI_Send(&np2_Thermal,3,MPI_UNSIGNED_LONG,this->MPI_communicator.rootProcess,TAG_NP2_THERMAL,MPI_COMM_WORLD);

        this->mygrid_Electro.id = this->MPI_communicator.getRank();
        this->mygrid_Thermal.id = this->MPI_communicator.getRank();

        for(int k = 0 ; k < 3 ; k ++){
            this->mygrid_Electro.np1[k] = np1_Electro[k];
            this->mygrid_Electro.np2[k] = np2_Electro[k];
            this->mygrid_Thermal.np1[k] = np1_Thermal[k];
            this->mygrid_Thermal.np2[k] = np2_Thermal[k];
        }
    }else{
        fprintf(stderr,"InterfaceToParaviewer::initializeAll::ERROR\n");
        fprintf(stderr,"Should not end up here. Aborting.\n");
        fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
        abort();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // The subgrid has been initialized on the root process.
    // Tell the object mygrid which are the vector and scalar fields:
    this->mygrid_Electro.vectors["ElectricField"] = NULL;
    this->mygrid_Electro.vectors["MagneticField"] = NULL;
    this->mygrid_Thermal.scalars["Temperature"]   = NULL;

    if(this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess){
        printf("InterfaceToParaviewer::initializeAll::OUT\n");
    }
}

/**
 * 
 * @brief Write the desired data to ouput files via vtl.
 * 
 * Detailed explanation goes here.
 */
void InterfaceToParaviewer::convertAndWriteData(unsigned long currentStep,
            std::string type /*"thermal" or "electro", case sensitive*/){

    /* FETCH THE OUPUT FILE NAME */
    map<std::string,std::string> outputFileNames;
    std::string outputName;

    if(this->is_grid_creator_new == false){
        // Is of type GridCreator:
        outputFileNames = this->grid_Creator.input_parser.get_outputNames();
        outputName = outputFileNames["output"];
    }else{
        // Is of type GridCreator_NEW:
        outputFileNames = this->grid_Creator_NEW.input_parser.get_outputNames();
        outputName = outputFileNames["output"];
    }

    cout << "Output VTK files will be named by " + outputName << endl;

    /* DETERMINE WHICH GRID TO SAVE */
    if(strcmp(type.c_str(),"thermal") == 0){
        /* SAVE THERMAL GRID */
        try{
            if(this->is_grid_creator_new == false){
                export_spoints_XML_custom_GridCreator(
                                "THERMAL",
                                outputName,
                                currentStep,
                                this->grid_Thermal, 
                                this->mygrid_Thermal,
                                this->grid_Creator, 
                                vtl::ZIPPED);
            }else{
                export_spoints_XML_custom_GridCreator_NEW(
                                "THERMAL",
                                outputName,
                                currentStep,
                                this->grid_Thermal, 
                                this->mygrid_Thermal,
                                this->grid_Creator_NEW, 
                                vtl::ZIPPED)
            }
        }catch(...){
            fprintf(stderr,"InterfaceToParaviewer::convertAndWriteData::ERROR\n");
            fprintf(stderr,"Error while writing output on process MPI %d. Aborting.\n",
                this->MPI_communicator.getRank());
            fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
            std::abort();
        }

        /* ONLY THE ROOT PROCESS CALLS THE FOLLOWING FUNCTION */
        if (this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess)
        {
            cout << ">>> MPI ROOT WRITING XMLP (THERMAL GRID)" << endl;
            if(this->is_grid_creator_new == false){
                export_spoints_XMLP_custom_GridCreator(
                                    "THERMAL",
                                    outputName, 
                                    currentStep, 
                                    this->grid_Thermal, 
                                    this->mygrid_Thermal, 
                                    this->sgrids_Thermal, 
                                    this->grid_Creator,
                                    vtl::ZIPPED);
            }else{
                export_spoints_XMLP_custom_GridCreator_NEW(
                                    "ELECTRO",
                                    outputName, 
                                    currentStep, 
                                    this->grid_Electro, 
                                    this->mygrid_Electro, 
                                    this->sgrids_Electro, 
                                    this->grid_Creator_NEW,
                                    vtl::ZIPPED);
            }
        }

        /* END OF THERMAL GRID SAVING */
    }else if(strcmp(type.c_str(),"electro") == 0){
        /* SAVE ELECTROMAGNETIC GRID */
        try{
            if(this->is_grid_creator_new == false){
                export_spoints_XML_custom_GridCreator(
                                "ELECTRO",
                                outputName,
                                currentStep,
                                this->grid_Electro, 
                                this->mygrid_Electro,
                                this->grid_Creator, 
                                vtl::ZIPPED);
            }else{
                export_spoints_XML_custom_GridCreator_NEW(
                                "ELECTRO",
                                outputName,
                                currentStep,
                                this->grid_Electro, 
                                this->mygrid_Electro,
                                this->grid_Creator_NEW, 
                                vtl::ZIPPED)
            }
        }catch(...){
            fprintf(stderr,"InterfaceToParaviewer::convertAndWriteData::ERROR\n");
            fprintf(stderr,"Error while writing output on process MPI %d. Aborting.\n",
                this->MPI_communicator.getRank());
            fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
            std::abort();
        }

        /* ONLY THE ROOT PROCESS CALLS THE FOLLOWING FUNCTION */
        if (this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess)
        {
            cout << ">>> MPI ROOT WRITING XMLP (ELECTRO GRID)" << endl;
            if(this->is_grid_creator_new == false){
                export_spoints_XMLP_custom_GridCreator(
                                    "ELECTRO",
                                    outputName, 
                                    currentStep, 
                                    this->grid_Electro, 
                                    this->mygrid_Electro, 
                                    this->sgrids_Electro, 
                                    this->grid_Creator,
                                    vtl::ZIPPED);
            }else{
                export_spoints_XMLP_custom_GridCreator_NEW(
                                    "ELECTRO",
                                    outputName, 
                                    currentStep, 
                                    this->grid_Electro, 
                                    this->mygrid_Electro, 
                                    this->sgrids_Electro, 
                                    this->grid_Creator_NEW,
                                    vtl::ZIPPED);
            }
        }

        /* END OF ELECTROMAGNETIC GRID SAVING */
    }else{
        /* WRONG PARAMETER */
        fprintf(stderr,"InterfaceToParaviewer::convertAndWriteData::ERROR\n");
        fprintf(stderr,"\t>>> Wrong saving type. Either 'thermal' or 'electro' (case insensitive).\n");
        fprintf(stderr,"\t>>> Received type=%s. Aborting.\n",type.c_str());
        fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
        abort();
    }

    cout << "WRITING DONE (MPI " << this->MPI_communicator.getRank() << endl;
}