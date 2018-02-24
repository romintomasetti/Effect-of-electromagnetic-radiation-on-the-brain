#include "AlgoElectro.h"

#include "Materials.h"
#include "GridCreator.h"
#include <fstream>
#include <cstring>
#include "ElectromagneticSource.h"
#include "Node3DField.h"

#include <unistd.h>

#include <omp.h>
#include <mpi.h>

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KNRM  "\x1B[0m"
#define KBLU  "\x1B[34m"
#define KYEL  "\x1B[33m"







void AlgoElectro::communicate(GridCreator& mesh, 
                                MPI_Initializer& MPI_communicator){
/*à faire*/

}

/* convention:Table array contient le maillage du MPI+ les voisins*/
/* convention: nbrElts ne contient pas les voisins*/


double AlgoElectro::Compute_dt(GridCreator &mesh){
    double dx=mesh.deltaX;
    double dy=mesh.deltaY;
    double dz=mesh.deltaZ;
    double dt=0.0;
    double tmp=0.0;
    double c=0.0;
    int i=0;                                                                                
    for (i=0;i<mesh.materials.numberOfMaterials;i++){

            //double temperature,unsigned char material, unsigned char property
            /* !!!!!!!!!!!!!!/*� faire avec T initial*/

            // Get material:
            string material = mesh.materials.materialName_FromMaterialID[i];

            double mu_material = mesh.materials.getProperty(
                    mesh.input_parser.GetInitTemp_FromMaterialName[material],
                    i,4);    

            /* !!!!!!!!!!!!!!/*� faire avec T initial*/
            double epsilon_material = mesh.materials.getProperty(
                   mesh.input_parser.GetInitTemp_FromMaterialName[material],
                    i,4);     

            c=1/(sqrt(mu_material*epsilon_material));

            if(i==0){
                dt=1/(c*sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));
            }
            else{
                tmp=1/(c*sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));
                if(tmp<dt){
                    dt=tmp;
                }
            }
    }
    return dt;
}

void AlgoElectro::update(GridCreator &mesh, InterfaceToParaviewer& interfaceForOutput){   
    /// Fetching final simulation time:
    double t_final   = mesh.input_parser.get_stopTime();
    /// Computing the optimal time step:
    double dt        = this->Compute_dt(mesh);


    /// Current simulation time:
    double t_current = 0.0;

    /// Temperature:
    double T = 0.0;
    /// Permeability:
    double mu_material = 0.0;
    /// Permeability:
    double epsilon_material = 0.0;
    /// Electrical conductivity:
    double elec_conduct_mat = 0.0;

    // Local and global numbers of the current node:
    unsigned long local[3];
    unsigned long global[3];


    ////////////////////////////////////////////////////////////
	/// WE NEED AT LEAST 6 OMP THREADS !!!                   ///
    /// In fact, we have 'parallélipipède rectangle' so      ///
    /// we have 6 simultaneous communications (one per face) ///
	////////////////////////////////////////////////////////////
	/*
    #pragma omp master
	{
		if(omp_get_max_threads() < 6){
			omp_set_num_threads(6);
			#if DEBUG > 0
				printf("%sFrom MPI P %d: I changed OMP_NUM_THREADS to %d !%s\n",
					KRED,mesh.MPI_communicator.getRank(),omp_get_max_threads(),KNRM);
			#endif
		}
	}
    */

    /////////////////////////////////////////////////////////////////////
    /// The following lines are usefull for MPI_MULTIPLE_THREAD usage ///
    /////////////////////////////////////////////////////////////////////
    /* REQ is the number of requests to check for. Taht is, 6 sends and 
     * 6 receives so 12 checks. */
    /*
    unsigned int REQ = 12;
    MPI_Request *requests = (MPI_Request*) calloc(REQ,sizeof(MPI_Request));
	int checkRequest[REQ];
	for(int cc = 0 ; cc < REQ ; cc ++ )
		checkRequest[cc] = -1;
	unsigned short int counterReq = 0;
    */
    //////////////////////////////////////
    /// BEGINNING OF THE OPENMP REGION ///
    //////////////////////////////////////
    /*
    #pragma omp parallel default(none)\
        shared(mesh)\
        shared(t_current,dt,t_final)\
        private(local,global)\
        private(T,mu_material,epsilon_material)
    {
        // Create for each OMP thread 4 vectors, to send and receive data:
        double *ElectricField_toRecv = NULL;
        double *ElectricField_toSend = NULL;
        double *MagneticField_toRecv = NULL;
        double *MagneticField_toSend = NULL;

        // Each OMP thread should be aware of the position of the MPI process:
        unsigned short int MPI_curr_X = 0;
        unsigned short int MPI_curr_Y = 0;
        unsigned short int MPI_curr_Z = 0;
        // Each OMP thread should know if there is a neighboor:        
        bool isNeighboor_at_UP    = false;
        bool isNeighboor_at_DOWN  = false;
        bool isNeighboor_at_SOUTH = false;
        bool isNeighboor_at_NORTH = false;
        bool isNeighboor_at_WEST  = false;
        bool isNeighboor_at_EAST  = false;
        // Boolean to know if the thread can send and receive:
        bool canSend = false;
        bool canRecv = false;
        // Boolean to know if the thread first sends or receives:
        bool firstSend = true;
        // Each OMP thread should know the 'coordinates' of the neighboor
        // with who it communicates:
        unsigned short int MPI_comm_X = 0;
        unsigned short int MPI_comm_Y = 0;
        unsigned short int MPI_comm_Z = 0;*/

        /* LINK AND INITIALIZE ALL VARIABLES WRITTEN ABOVE */
       /* if(omp_get_thread_num() == 0){
            if(isNeighboor_at_UP == true){
                // 
            }
        }*/

        /////////////////////////////////////////////////
        /// CONVENTION FOR COMMUNICATION              ///
        /// 1) OMP thread(0) communicates with UP.    ///
        /// 2) OMP thread(1) communicates with DOWN.  ///
        /// 3) OMP thread(2) communicates with SOUTH. ///
        /// 3) OMP thread(3) communicates with NORTH. ///
        /// 4) OMP thread(4) communicates with WEST.  ///
        /// 5) OMP thread(5) communicates with EAST.  ///
        /////////////////////////////////////////////////

       /* printf("Inside algoElectro update :: aborting()\n");
        abort();*/
        int counter = 0;
        while(t_current<t_final){

            /* AFFICHER LA TRANCHE Z=8 */
            int tranche = 8;
            printf("\n\t***************AFFICHAGE DE LA TRANCHE Z = %d***********\n\n",tranche);

            

            for(unsigned long i=1; i <= mesh.numberOfNodesInEachDir[0] ;i++){
                    for(unsigned long j=1; j <= mesh.numberOfNodesInEachDir[1]  ;j++){
                        if(mesh.nodesElec(i,j,tranche).field[2] == 0){
                            printf("%s%f%s,",KGRN,mesh.nodesElec(i,j,tranche).field[2],KNRM);
                        }else if(mesh.nodesElec(i,j,tranche).field[2] < 0){
                            printf("%s%f%s,",KBLU,mesh.nodesElec(i,j,tranche).field[2],KNRM);
                        }else{
                            printf("%s%f%s,",KRED,mesh.nodesElec(i,j,tranche).field[2],KNRM);
                        }
                    }
                    cout << "END" << endl;
            }
            ///////////////////////////
            /* UPDATE MAGNETIC FIELD */
            ///////////////////////////
            for(unsigned long k = 0 ; k <= mesh.numberOfNodesInEachDir[2] ; k++ ){

                for(unsigned long j=0;j <= mesh.numberOfNodesInEachDir[1] ; j++){

                    for(unsigned long i=0;i <= mesh.numberOfNodesInEachDir[0]  ; i++){
                        printf("%s\t>>> NODE(%ld,%ld,%ld) <<<%s\n",KGRN,i,j,k,KNRM);
                        /* Get the global indices */
                        local[0] = i;
                        local[1] = j;
                        local[2] = k;
                        mesh.LocalToGlobal(local,global);

                        

                        T = mesh.nodesMagn(i,j,k).Temperature;
                        printf("> Fetching temperature: %f\n",T);

                        printf("> Fetching mu_material for material %d.\n",mesh.nodesMagn(i,j,k).material);

                        mu_material = mesh.materials.getProperty(T,mesh.nodesMagn(i,j,k).material,4);  

                        printf("> Fetched mu_material for material %d is %f.\n",
                                                mesh.nodesMagn(i,j,k).material,mu_material);


                        elec_conduct_mat = mesh.materials.getProperty(T,mesh.nodesMagn(i,j,k).material,6);
                        
                        /* mettre dans gridcreator "numberofProcess" et " myrank"  */
                        //Transfo = mesh.LocalToGlobal(i,j,k,mesh.MPI_communicator.getRank(),mesh.numberofprocess) ;  
                        printf("> isInsideSource...");
                        mesh.input_parser.source.isInsideSource(global[0],global[1],global[2]);
                        printf("Done.\n");

                        if(mesh.input_parser.source.isInsideSource(global[0],global[1],global[2])){                
                            /* Fonction à faire  dans Electromagnetic Source composantes donne l'info si c'est E ou H et x ou y ou z */
                            printf("> Fetching source value...");
                            mesh.input_parser.source.computeSourceValue(mesh,global[0],global[1],global[2],t_current,'H');
                            printf("Done.\n");
                            /*mesh.nodesMagn(i,j,k).field[1]=mesh.input_parser.source.computeSourceValue(mesh,i,j,k,t_current);
                            mesh.nodesMagn(i,j,k).field[2]=mesh.input_parser.source.computeSourceValue(mesh,i,j,k,t_current);*/
                        }
                        else{

                            /* update magnetic field H_x */

                            //////////////////////////////////////
                            //////////////////////////////////////////////////////////////////////////////////////////
                            ///////////////// achtung attention  attentionze /////
                            ////// ajouter SIGNMA /// (see formulas)


                            mesh.nodesMagn(i,j,k).field[0] = mesh.nodesMagn(i,j,k).field[0]  +

                                            (dt/(mu_material*mesh.deltaZ))*(mesh.nodesElec(i,j,k+1).field[1]
                                            -mesh.nodesElec(i,j,k).field[1]) -

                                            (dt/(mu_material*mesh.deltaY))*(mesh.nodesElec(i,j+1,k).field[2]
                                            -mesh.nodesElec(i,j,k).field[2]);

                            /* update magnetic Field H_y */
                            mesh.nodesMagn(i,j,k).field[1]= mesh.nodesMagn(i,j,k).field[1]  +
                                            (dt/(mu_material*mesh.deltaX))*(mesh.nodesElec(i+1,j,k).field[2]
                                            -mesh.nodesElec(i,j,k).field[2]) -

                                            (dt/(mu_material*mesh.deltaZ))*(mesh.nodesElec(i,j,k+1).field[0]
                                            -mesh.nodesElec(i,j,k).field[0]);
                            
                            /* update magnetic Field H_z */
                            mesh.nodesMagn(i,j,k).field[2]= mesh.nodesMagn(i,j,k).field[2]  +
                                            (dt/(mu_material*mesh.deltaY))*(mesh.nodesElec(i,j+1,k).field[0]-
                                            mesh.nodesElec(i,j,k).field[0]) -

                                            (dt/(mu_material*mesh.deltaX))*(mesh.nodesElec(i+1,j,k).field[1]-
                                            mesh.nodesElec(i,j,k).field[1]);
                        }
                    }
                }
            }

            /* update electric field  */
            
            printf("%sELECTRIC FIELD%s\n",KRED,KNRM);

            for(unsigned long k=1 ; k <= mesh.numberOfNodesInEachDir[0];k++){
                for(unsigned long j=1 ; j <= mesh.numberOfNodesInEachDir[1];j++){
                    for(unsigned long i=1; i <= mesh.numberOfNodesInEachDir[2];i++){
                        printf("%s\t>>> NODE(%ld,%ld,%ld) <<<%s\n",KGRN,i,j,k,KNRM);
                        /* Get the global indices */
                        local[0] = i;
                        local[1] = j;
                        local[2] = k;
                        mesh.LocalToGlobal(local,global);

                        T = mesh.nodesElec(i,j,k).Temperature;

                        epsilon_material = mesh.materials.getProperty(T,mesh.nodesElec(i,j,k).material,5);

                        elec_conduct_mat = mesh.materials.getProperty(T,mesh.nodesElec(i,j,k).material,6);

                        if(mesh.input_parser.source.isInsideSource(global[0],global[1],global[2])){

                            mesh.input_parser.source.computeSourceValue(mesh, global[0],global[1],global[2],t_current,'E');
                            //mesh.nodesMagn(i,j,k).field[1]=mesh.input_parser.source.computeSourceValue(mesh, i,j,k,t_current,composants_5);
                            //mesh.nodesMagn(i,j,k).field[2]=mesh.input_parser.source.computeSourceValue(mesh, i,j,k,t_current,composants_6);
                        }else{
                            /* update magnetic field E_x */
                            mesh.nodesElec(i,j,k).field[0]= mesh.nodesElec(i,j,k).field[0] +
                            
                                                 (dt/(epsilon_material*mesh.deltaY))*(mesh.nodesMagn(i,j,k).field[2]-
                                                 mesh.nodesMagn(i,j-1,k).field[2]) - 

                                                 (dt/(epsilon_material*mesh.deltaZ))*(mesh.nodesMagn(i,j,k).field[1]-
                                                 mesh.nodesMagn(i,j,k-1).field[1]);
                            /* update magnetic Field E_y */
                            mesh.nodesElec(i,j,k).field[1]= mesh.nodesElec(i,j,k).field[1] + 
                                                (dt/(epsilon_material*mesh.deltaZ))*(mesh.nodesMagn(i,j,k).field[0]-
                                                mesh.nodesMagn(i,j,k-1).field[0])  - 

                                                (dt/(epsilon_material*mesh.deltaX))*(mesh.nodesMagn(i,j,k).field[2]-
                                                mesh.nodesMagn(i-1,j,k).field[2]);
                            /* update magnetic Field E_z */
                            mesh.nodesElec(i,j,k).field[2]= mesh.nodesElec(i,j,k).field[2] +
                                                 (dt/(epsilon_material*mesh.deltaX))*(mesh.nodesMagn(i,j,k).field[1]-
                                                 mesh.nodesMagn(i-1,j,k).field[1]) -
                                                 
                                                  (dt/(epsilon_material*mesh.deltaY))*(mesh.nodesMagn(i,j,k).field[0]-
                                                  mesh.nodesMagn(i,j-1,k).field[0]);
                        }
                    }
                }
            }   
            

            /* WRITE RESULTS */
            interfaceForOutput.convertAndWriteData(this->currentStep);

            ////////////////////////////////////////
            /// UPDATING CURRENT SIMULATION TIME ///
            ////////////////////////////////////////
            t_current=t_current+dt;
            this->currentStep ++;
        }

        printf("ALGO::END::DT=%f::FINAL_TIME=%f\n",dt,t_current);
   /* }*/

    //// FREEING MEMORY ////
   /* free(requests);*/
}


AlgoElectro::~AlgoElectro(void){
    cout << "Calling AlgoCreator destructor\n";
}


