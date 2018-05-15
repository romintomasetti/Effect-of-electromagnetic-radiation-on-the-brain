#include "InterfaceToParaviewer.h"
#include "vtlSPoints_romin.h"

#include "unistd.h"

#include <string>
#include <string.h>
#include <cstring>
#include <map>

#define TAG_NP1_ELECTRO 20
#define TAG_NP2_ELECTRO 21
#define TAG_NP1_THERMAL 22
#define TAG_NP2_THERMAL 23

// Include some libraries for directory support:
#ifdef _WIN32
#include <direct.h>
// MSDN recommends against using getcwd & chdir names
#define cwd _getcwd
#define cd _chdir
#else
#include "unistd.h"
#define cwd getcwd
#define cd chdir
#endif

#include<sys/stat.h>

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


    dx_Electro = this->grid_Creator_NEW.input_parser.deltaX_Electro;
    dy_Electro = this->grid_Creator_NEW.input_parser.deltaY_Electro;
    dz_Electro = this->grid_Creator_NEW.input_parser.deltaZ_Electro;

    delta_thermal = this->grid_Creator_NEW.input_parser.delta_Thermal;

    /* CHECK DELTAS */
    if(dx_Electro <= 0.0 || dy_Electro <= 0.0 || dz_Electro <= 0.0 || delta_thermal <= 0.0){
        fprintf(stderr,"InterfaceToParaviewer::initializeAll::ERROR :: One of the spatial step for the electromagnetic/thermal mesh is negative or zeros.\n");
        fprintf(stderr,"InterfaceToParaviewer::initializeAll::ABORTING (line %d, file %s).\n",__LINE__,__FILE__);
        abort();
    }

    /* INITIALIZE 'grid_thermal'/'grid_Electro' for each MPI process (even for non-root ones) */
    this->grid_Electro.o = this->grid_Creator_NEW.originOfWholeSimulation_Electro;
    this->grid_Thermal.o = this->grid_Creator_NEW.originOfWholeSimulation_Thermal;

    this->grid_Electro.dx = vtl_romin::Vec3d(dx_Electro,dy_Electro,dz_Electro);
    this->grid_Thermal.dx = vtl_romin::Vec3d(delta_thermal,delta_thermal,delta_thermal);
    
    this->grid_Electro.np1 = vtl_romin::Vec3i(0,0,0); 
    this->grid_Thermal.np1 = vtl_romin::Vec3i(0,0,0);

    double length_X_whole_dom_electro = 0.0;
    double length_Y_whole_dom_electro = 0.0;
    double length_Z_whole_dom_electro = 0.0;
    double length_x_whole_dom_thermal = 0.0;
    double length_y_whole_dom_thermal = 0.0;
    double length_z_whole_dom_thermal = 0.0;

    length_X_whole_dom_electro = this->grid_Creator_NEW.input_parser.lengthX_WholeDomain_Electro;
    length_Y_whole_dom_electro = this->grid_Creator_NEW.input_parser.lengthY_WholeDomain_Electro;
    length_Z_whole_dom_electro = this->grid_Creator_NEW.input_parser.lengthZ_WholeDomain_Electro;
    length_x_whole_dom_thermal = this->grid_Creator_NEW.input_parser.lengthX_WholeDomain_Thermal;
    length_y_whole_dom_thermal = this->grid_Creator_NEW.input_parser.lengthY_WholeDomain_Thermal;
    length_z_whole_dom_thermal = this->grid_Creator_NEW.input_parser.lengthZ_WholeDomain_Thermal;

    size_t nodesWholeDom_X_Electro = length_X_whole_dom_electro / dx_Electro +1;
    size_t nodesWholeDom_Y_Electro = length_Y_whole_dom_electro / dy_Electro +1;
    size_t nodesWholeDom_Z_Electro = length_Z_whole_dom_electro / dz_Electro +1;

    this->grid_Electro.np2 = vtl_romin::Vec3i(nodesWholeDom_X_Electro,nodesWholeDom_Y_Electro,nodesWholeDom_Z_Electro);

    size_t nodesWholeDom_X_Thermal = length_x_whole_dom_thermal / delta_thermal +1;
    size_t nodesWholeDom_Y_Thermal = length_y_whole_dom_thermal / delta_thermal +1;
    size_t nodesWholeDom_Z_Thermal = length_z_whole_dom_thermal / delta_thermal +1;
    
    this->grid_Thermal.np2 = vtl_romin::Vec3i(nodesWholeDom_X_Thermal,nodesWholeDom_Y_Thermal,nodesWholeDom_Z_Thermal);

    // Initialize the subgrid if I am the root process:
    if(this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess){
        #ifndef NDEBUG
            printf("InterfaceToParaviewer::initializeAll\n");
            printf("\t> From MPI %d, I am the root !\n",this->MPI_communicator.getRank());
            printf("\t> Initializing the subgrids (EM and TH) of class-type SPoints with %d processes.\n",
                        nb_MPI);
        #endif

        // Allocating space:
        this->sgrids_Electro.resize(nb_MPI);
        this->sgrids_Thermal.resize(nb_MPI);

        // Loop over each grid and get the starting/ending indices:
        for(unsigned char I = 0 ; I < nb_MPI ; I ++){

            #ifndef NDEBUG
                printf("\t> MPI %d :: Initializing sgrid_electro[%d] and sgrid_thermal[%d]...\n",
                        this->MPI_communicator.getRank(),I,I);
            #endif

            // Giving the MPI ID to the element sgrids[I]:
            this->sgrids_Electro[I].id = I;
            this->sgrids_Thermal[I].id = I;

            // Giving the spacing to the element sgrids[I]:
            this->sgrids_Electro[I].dx = vtl_romin::Vec3d(dx_Electro,dy_Electro,dz_Electro);
            this->sgrids_Thermal[I].dx = vtl_romin::Vec3d(delta_thermal,delta_thermal,delta_thermal);

            // Setting the origin:
            this->sgrids_Electro[I].o = this->grid_Creator_NEW.originOfWholeSimulation_Electro;
            this->sgrids_Thermal[I].o = this->grid_Creator_NEW.originOfWholeSimulation_Thermal;

            // Setting the origin and end indices of each subgrid:
            if(I == this->MPI_communicator.rootProcess){
                this->mygrid_Electro.id = this->MPI_communicator.getRank();
                this->mygrid_Thermal.id = this->MPI_communicator.getRank();
                // Don't communicate:
                for(int k = 0 ; k < 3 ; k ++){
                    
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
                            + this->mygrid_Thermal.np2[k];

                }
            }else{
                unsigned long np1_Electro[3];
                unsigned long np1_Thermal[3];
                unsigned long np2_Electro[3];
                unsigned long np2_Thermal[3];

                // Receive np1_Electro:
                MPI_Recv(&np1_Electro,3,MPI_UNSIGNED_LONG,I,TAG_NP1_ELECTRO,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                // Receive np2_Electro:
                MPI_Recv(&np2_Electro,3,MPI_UNSIGNED_LONG,I,TAG_NP2_ELECTRO,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                // Receive np1_Thermal:
                MPI_Recv(&np1_Thermal,3,MPI_UNSIGNED_LONG,I,TAG_NP1_THERMAL,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                // Receive np2_Thermal:
                MPI_Recv(&np2_Thermal,3,MPI_UNSIGNED_LONG,I,TAG_NP2_THERMAL,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                for(int k = 0 ; k < 3 ; k ++){
                    this->sgrids_Electro[I].np1[k] = np1_Electro[k];
                    this->sgrids_Electro[I].np2[k] = np2_Electro[k];

                    this->sgrids_Thermal[I].np1[k] = np1_Thermal[k];
                    this->sgrids_Thermal[I].np2[k] = np2_Thermal[k];
                }
            }

            //cout << this->sgrids_Electro[I] << endl;
            //cout << this->sgrids_Thermal[I] << endl;
        }

    }else if(this->MPI_communicator.isRootProcess() == INT_MIN){

        #ifndef NDEBUG
            printf("InterfaceToParaviewer::initializeAll\n");
            printf("\t> From MPI process %d, I am not the root.\n",this->MPI_communicator.getRank());
        #endif

        unsigned long np1_Electro[3];
        unsigned long np2_Electro[3];
        unsigned long np1_Thermal[3];
        unsigned long np2_Thermal[3];

        for(int k = 0 ; k  <  3 ; k ++ ){
                
            np1_Electro[k] = this->grid_Creator_NEW.originIndices_Electro[k];
            np2_Electro[k] = np1_Electro[k] + this->grid_Creator_NEW.sizes_EH[k];

            np1_Thermal[k] = this->grid_Creator_NEW.originIndices_Thermal[k];
            np2_Thermal[k] = np1_Thermal[k] + this->grid_Creator_NEW.size_Thermal[k];
        }   

        // Send np1_Electro:
        MPI_Send(&np1_Electro,3,MPI_UNSIGNED_LONG,this->MPI_communicator.rootProcess,TAG_NP1_ELECTRO,MPI_COMM_WORLD);

        // Send np2_Electro:
        MPI_Send(&np2_Electro,3,MPI_UNSIGNED_LONG,this->MPI_communicator.rootProcess,TAG_NP2_ELECTRO,MPI_COMM_WORLD);

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

    #ifndef NDEBUG
        if(this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess){
            printf("InterfaceToParaviewer::initializeAll::OUT\n");
        }
    #endif
}

/**
 * @brief This function returns a folder name contained inside a string and the output file's name.
 */
std::vector<std::string> get_folder_from_name(std::string outputName){
    
    // First element is the folder name, second is the file name.
    std::vector<std::string> returned_folder_name = {string(),outputName};

    if(outputName.find('/') != std::string::npos){
        /* The output name specifies an output folder name */

        // Folder name:
        returned_folder_name[0] = outputName.substr(0,outputName.find('/'));

        // File name:
        returned_folder_name[1] = outputName.substr(outputName.find('/')+1);
    }

    return returned_folder_name;
}

/**
 * @brief This function creates the gien folder and goes into it (changes the current working directory).
 * 
 *  It also returns the parent folder's name.
 */
std::string InterfaceToParaviewer::create_folder_and_go_in(std::string folderName){

    char buf[4096];

    // Get current working directory:
    std::string currentWorkingDir = cwd(buf, sizeof buf);

    #ifndef NDEBUG
        printf("Current working directory is %s.\n",currentWorkingDir.c_str());
    #endif

    // Try to go into 'folderName' directory:
    if (0 == cd(folderName.c_str())) {
        #ifndef NDEBUG
            std::cout << "CWD changed to: " << cwd(buf, sizeof buf) << std::endl;
        #endif
    }else{
        if(this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess){
            #ifndef NDEBUG
                printf("Must create the directory %s...\n",folderName.c_str());
            #endif
            #if defined(_WIN32)
                _mkdir(folderName.c_str());
            #else 
                mkdir(folderName.c_str(), 0700); 
            #endif
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(0 == cd(folderName.c_str())){
            #ifndef NDEBUG
                std::cout << "CWD changed to: " << cwd(buf, sizeof(buf)) << std::endl;
            #endif
        }else{
            DISPLAY_ERROR_ABORT_CLASS(
                "Cannot create/change directory %s !",
                folderName.c_str()
            );
        }
	}

    return currentWorkingDir;

}


/**
 * 
 * @brief Write the desired data to ouput files via vtl.
 * 
 * Detailed explanation goes here.
 */
void InterfaceToParaviewer::convertAndWriteData(unsigned long currentStep,
            std::string type /*"thermal" or "electro", case sensitive*/,
            algoElectroToVtlRomin pmlVersSauvegarde){

    /* FETCH THE OUPUT FILE NAME */
    map<std::string,std::string> outputFileNames;
    std::string outputName;

    outputFileNames = this->grid_Creator_NEW.input_parser.get_outputNames();
    outputName = outputFileNames["output"];

    /// Determine if the output files should be stored in a folder:
    std::vector<std::string> folderAndFileNames = get_folder_from_name(outputName);

    outputName = folderAndFileNames[1];

    // Will contain the current folder:
    std::string parentFolder = string();

    if(folderAndFileNames[0] != std::string()){
        /* A folder name was specified, go into it !*/
        parentFolder = this->create_folder_and_go_in(folderAndFileNames[0]);
    }

    /* DETERMINE WHICH GRID TO SAVE */
    if(strcmp(type.c_str(),"THERMAL") == 0){
        /* SAVE THERMAL GRID */

        export_spoints_XML_custom_GridCreator_NEW(
                        "THERMAL",
                        outputName,
                        currentStep,
                        this->grid_Thermal, 
                        this->mygrid_Thermal,
                        this->grid_Creator_NEW, 
                        vtl_romin::ZIPPED,
                        pmlVersSauvegarde);

        /* ONLY THE ROOT PROCESS CALLS THE FOLLOWING FUNCTION */
        if (this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess)
        {
            
            export_spoints_XMLP_custom_GridCreator_NEW(
                                "THERMAL",
                                outputName, 
                                currentStep, 
                                this->grid_Thermal, 
                                this->mygrid_Thermal, 
                                this->sgrids_Thermal, 
                                //this->grid_Creator_NEW,
                                vtl_romin::ZIPPED);
        }

        /* END OF THERMAL GRID SAVING */
    }else if(strcmp(type.c_str(),"ELECTRO") == 0){
        /* SAVE ELECTROMAGNETIC GRID */
            
        export_spoints_XML_custom_GridCreator_NEW(
                        "ELECTRO",
                        outputName,
                        currentStep,
                        this->grid_Electro, 
                        this->mygrid_Electro,
                        this->grid_Creator_NEW, 
                        vtl_romin::ZIPPED,
                        pmlVersSauvegarde);

        /* ONLY THE ROOT PROCESS CALLS THE FOLLOWING FUNCTION */
        if (this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess)
        {            
            export_spoints_XMLP_custom_GridCreator_NEW(
                                "ELECTRO",
                                outputName, 
                                currentStep, 
                                this->grid_Electro, 
                                this->mygrid_Electro, 
                                this->sgrids_Electro, 
                                //this->grid_Creator_NEW,
                                vtl_romin::ZIPPED);
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

    /// Revert to parent working directory:
    if( 0 != cd(parentFolder.c_str())){
        fprintf(stderr,"In %s :: could not revert to parent directory %s. Aborting.\n",
            __FUNCTION__,parentFolder.c_str());
        fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
        #ifdef MPI_COMM_WORLD
            MPI_Abort(MPI_COMM_WORLD,-1);
        #else
            abort();
        #endif
    }
}