#include "InterfaceToParaviewer.h"
#include "vtl/vtlSPoints.h"

#include "unistd.h"

#include <string>
#include <map>

#define TAG_NP1 20
#define TAG_NP2 21

// Initilize everything:
void InterfaceToParaviewer::initializeAll(void){
    /* MAIN GOAL : create a vector of subgrids */

    // Retrieve the number of MPI processes:
    int nb_MPI = this->MPI_communicator.getNumberOfMPIProcesses();

    double dx = this->grid_Creator.input_parser.deltaX;
    double dy = this->grid_Creator.input_parser.deltaY;
    double dz = this->grid_Creator.input_parser.deltaZ;

    /* INITIALIZE 'grid' for each MPI process (even for non-root ones) */
    this->grid.o = this->grid_Creator.originOfWholeSimulation;
    this->grid.dx = vtl::Vec3d(dx,dy,dz);
    this->grid.np1 = vtl::Vec3i(0,0,0);
    size_t nodesWholeDom_X = (size_t) this->grid_Creator.input_parser.lengthX /
                    this->grid_Creator.input_parser.deltaX +1;
    size_t nodesWholeDom_Y = (size_t) this->grid_Creator.input_parser.lengthY /
                    this->grid_Creator.input_parser.deltaY +1;
    size_t nodesWholeDom_Z = (size_t) this->grid_Creator.input_parser.lengthZ /
                    this->grid_Creator.input_parser.deltaZ +1;
    this->grid.np2 = vtl::Vec3i(nodesWholeDom_X,nodesWholeDom_Y,nodesWholeDom_Z);



    MPI_Barrier(MPI_COMM_WORLD);
    // Initialize the subgrid if I am the root process:
    if(this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess){
        printf("InterfaceToParaviewer::initializeAll\n");
        printf("\t> From MPI process %d, I am the root !\n",this->MPI_communicator.getRank());
        printf("\t> Initializing the subgrid of class-type SPoints with %d processes.\n",
                    nb_MPI);

        // Resize the sgrid field:
        this->sgrids.resize(nb_MPI);
        // Loop over each grid and get the starting/ending indices:
        for(unsigned char I = 0 ; I < nb_MPI ; I ++){
            printf("\t> Initializing sgrid[%d]...\n",I);
            // Giving the MPI ID to the element sgrids[I]:
            this->sgrids[I].id = I;
            // Giving the spacing to the element sgrids[I]:
            this->sgrids[I].dx = vtl::Vec3d(dx,dy,dz);
            // Setting the origin:
            this->sgrids[I].o  = this->grid_Creator.originOfWholeSimulation;
            // Setting the origin and end inidices of each subgrid:
            if(I == this->MPI_communicator.rootProcess){
                // Don't communicate:
                for(int k = 0 ; k < 3 ; k ++){
                    this->sgrids[I].np1[k] = this->grid_Creator.originIndices[k];
                    this->sgrids[I].np2[k] = this->sgrids[I].np1[k] + this->grid_Creator.numberOfNodesInEachDir[k];
                }
            }else{
                unsigned long np1[3], np2[3];
                MPI_Recv(&np1,3,MPI_UNSIGNED_LONG,I,TAG_NP1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                printf("\t> From root, received np1 from %d.\n",I);
                MPI_Recv(&np2,3,MPI_UNSIGNED_LONG,I,TAG_NP2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                printf("\t> From root, received np2 from %d.\n",I);
                for(int k = 0 ; k < 3 ; k ++){
                    this->sgrids[I].np1[k] = np1[k];
                    this->sgrids[I].np2[k] = np2[k];
                }
            }

            cout << this->sgrids[I] << endl;

            printf("\t> MPI %d goes from (%ld,%ld,%ld) to (%ld,%ld,%ld)\n",I,this->sgrids[I].np1[0],
                                                                            this->sgrids[I].np1[1],
                                                                            this->sgrids[I].np1[2],
                                                                            this->sgrids[I].np2[0],
                                                                            this->sgrids[I].np2[1],
                                                                            this->sgrids[I].np2[2]);
        }
    }else if(this->MPI_communicator.isRootProcess() == INT_MIN){
        printf("InterfaceToParaviewer::initializeAll\n");
        printf("\t> From MPI process %d, I am not the root.\n",this->MPI_communicator.getRank());
        unsigned long np1[3],np2[3];
        for(int k = 0 ; k  <  3 ; k ++ ){
            np1[k] = this->grid_Creator.originIndices[k];
            np2[k] = np1[k] + this->grid_Creator.numberOfNodesInEachDir[k];
            cout << "np2[2]=" << np2[2] << endl;
        }
        MPI_Send(&np1,3,MPI_UNSIGNED_LONG,this->MPI_communicator.rootProcess,TAG_NP1,MPI_COMM_WORLD);
        printf("\t> From MPI process %d, np1 sent to root.\n",this->MPI_communicator.getRank());
        MPI_Send(&np2,3,MPI_UNSIGNED_LONG,this->MPI_communicator.rootProcess,TAG_NP2,MPI_COMM_WORLD);
        printf("\t> From MPI process %d, np2 sent to root.\n",this->MPI_communicator.getRank());
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // The subgrid has been initialized on the root process.
    // Tell the object mygrid which are the vector and scalar fields:
    this->mygrid.vectors["ElectricField"] = NULL;
    this->mygrid.vectors["MagneticField"] = NULL;
    this->mygrid.scalars["Temperature"] = NULL;

    if(this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess){
        printf("InterfaceToParaviewer::initializeAll::OUT\n");
    }
}

// Convert and write output:
void InterfaceToParaviewer::convertAndWriteData(unsigned long currentStep){
    // FETCH THE OUTPUT FILE NAME:
    MPI_Barrier(MPI_COMM_WORLD);
    sleep(1);
    printf("HELL I AM PROCESS %d.\n",this->MPI_communicator.getRank());
    map<std::string,std::string> outputFileNames = this->grid_Creator.input_parser.get_outputNames();
    std::string outputName = outputFileNames["output"];
    cout << "Output VTK files will be named by " + outputName << endl;

    /* ALL MPI PROCESSES CALL THE FOLLOWING SAVING FUNCTION */
    try{
        export_spoints_XML(outputName, currentStep, this->grid, this->mygrid,this->grid_Creator, vtl::ZIPPED);
    }catch(...){
        printf("Error while writing output on process MPI %d. Aborting.\n",
                this->MPI_communicator.getRank());
        std::abort();
    }

    /* ONLY THE ROOT PROCESS CALLS THE FOLLOWING FUNCTION */
    if (this->MPI_communicator.isRootProcess() == this->MPI_communicator.rootProcess)
    {
        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>MPI ROOT WRITING XMLP" << endl;
        export_spoints_XMLP(outputName, currentStep, this->grid, this->mygrid, this->sgrids, vtl::ZIPPED);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}