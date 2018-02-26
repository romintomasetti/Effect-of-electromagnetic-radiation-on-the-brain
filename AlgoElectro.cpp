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

    double deltaX = mesh.deltaX;
    double deltaY = mesh.deltaY;
    double deltaZ = mesh.deltaZ;

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
    /* REQ is the number of requests to check for. That is, 6 sends and 
     * 6 receives so 12 checks. */
    
    /*
    unsigned int REQ_MPI = 12;
    MPI_Request *requests_MPI = (MPI_Request*) calloc(REQ_MPI,sizeof(MPI_Request));
	int checkRequest[REQ_MPI];
	for(int cc = 0 ; cc < REQ_MPI ; cc ++ )
		checkRequest[cc] = -1;
	unsigned short int counterReq = 0;
    */

    /////////////////////////////////////////////////////
    /// INITIALIZING NEIGHBOORS FOR MPI COMMUNICATION ///
    /////////////////////////////////////////////////////
    /*
    RankNeighbour[0] = SOUTH (along the opposite direction of the x-axis)
	RankNeighbour[1] = NORTH (along the direction of the x-axis)
	RankNeighbour[2] = WEST (along the opposite direction of the y-axis)
	RankNeighbour[3] = EAST (along the direction of the y-axis)
	RankNeighbour[4] = DOWN (along the opposite direction of the z-axis)
	RankNeighbour[5] = UP (along the direction of the z-axis) 
    */
    /*
    bool isNeighboor_at_UP    = false;
    bool isNeighboor_at_DOWN  = false;
    bool isNeighboor_at_SOUTH = false;
    bool isNeighboor_at_NORTH = false;
    bool isNeighboor_at_WEST  = false;
    bool isNeighboor_at_EAST  = false;

    if(mesh.MPI_communicator.RankNeighbour[0] != -1){isNeighboor_at_SOUTH = true;}
    if(mesh.MPI_communicator.RankNeighbour[1] != -1){isNeighboor_at_NORTH = true;}
    if(mesh.MPI_communicator.RankNeighbour[2] != -1){isNeighboor_at_WEST  = true;}
    if(mesh.MPI_communicator.RankNeighbour[3] != -1){isNeighboor_at_EAST  = true;}
    if(mesh.MPI_communicator.RankNeighbour[4] != -1){isNeighboor_at_DOWN  = true;}
    if(mesh.MPI_communicator.RankNeighbour[5] != -1){isNeighboor_at_UP    = true;}
    */
    
    //////////////////////////////////////
    /// BEGINNING OF THE OPENMP REGION ///
    //////////////////////////////////////
    /*
    #pragma omp parallel default(none)\
        shared(ompi_mpi_comm_world)\
        shared(mesh,interfaceForOutput)\
        shared(requests_MPI,checkRequest,REQ_MPI,counterReq)\
        shared(t_current,dt,t_final)\
        shared(isNeighboor_at_DOWN,isNeighboor_at_EAST,isNeighboor_at_NORTH)\
        shared(isNeighboor_at_SOUTH,isNeighboor_at_UP,isNeighboor_at_WEST)\
        private(local,global)\
        private(T,mu_material,epsilon_material)\
        private(elec_conduct_mat)
    {
    */

        // Create for each OMP thread 2 vectors, to send and receive data:
        /*
        Node3DField **ElectricNodes_toSend = NULL;
        size_t size1_send, size2_send;
        Node3DField **ElectricNodes_toRecv = NULL;
        size_t size1_recv, size2_recv;
        */

        // Each thread has its direction:
        /*
        char direction;
        */

        // If there is a neighboor at Down, we must send to and receive from Down.
        // So we must initialize both _toSend and _toRecv Node3DField arrays.
        /////////////////////////////////////////////////
        /// CONVENTION FOR COMMUNICATION              ///
        /// 1) OMP thread(0) communicates with SOUTH. ///
        /// 2) OMP thread(1) communicates with NORTH. ///
        /// 3) OMP thread(2) communicates with WEST.  ///
        /// 3) OMP thread(3) communicates with EAST.  ///
        /// 4) OMP thread(4) communicates with DOWN.  ///
        /// 5) OMP thread(5) communicates with UP.    ///
        /////////////////////////////////////////////////
        /*
        if(omp_get_thread_num() == 0){
            // Communication with SOUTH.
            direction = 'S';
            if(isNeighboor_at_SOUTH == true){
                counterReq += 2;
                checkRequest[0] = 1;
                checkRequest[1] = 1;
                // There is a neighboor at south. Initialize arrays.
                size1_send = mesh.numberOfNodesInEachDir[1];
                size2_send = mesh.numberOfNodesInEachDir[2];
                size1_recv = size1_send;
                size2_recv = size2_send;
                ElectricNodes_toSend = new Node3DField*[size1_send];
                ElectricNodes_toRecv = new Node3DField*[size1_recv];
                for(size_t I = 0 ; I < size1_send ; I ++)
                    ElectricNodes_toSend[I] = new Node3DField[size2_send];
                for(size_t I = 0 ; I < size1_recv ; I ++)
                    ElectricNodes_toRecv[I] = new Node3DField[size2_recv];
            }else{
                counterReq += 2;
            }
        }
        if(omp_get_thread_num() == 1){
            // Communication with NORTH.
            direction = 'N';
            if(isNeighboor_at_NORTH == true){
                counterReq += 2;
                checkRequest[2] = 1;
                checkRequest[3] = 1;
                // There is a neighboor at north. Initialize arrays.
                size1_send = mesh.numberOfNodesInEachDir[1];
                size2_send = mesh.numberOfNodesInEachDir[2];
                size1_recv = size1_send;
                size2_recv = size2_send;
                ElectricNodes_toSend = new Node3DField*[size1_send];
                ElectricNodes_toRecv = new Node3DField*[size1_recv];
                for(size_t I = 0 ; I < size1_send ; I ++)
                    ElectricNodes_toSend[I] = new Node3DField[size2_send];
                for(size_t I = 0 ; I < size1_recv ; I ++)
                    ElectricNodes_toRecv[I] = new Node3DField[size2_recv];
            }else{
                counterReq += 2;
            }
        }
        if(omp_get_thread_num() == 2){
            // Communication with West.
            direction = 'W';
            if(isNeighboor_at_WEST){
                counterReq += 2;
                checkRequest[4] = 1;
                checkRequest[5] = 1;
                // There is a neighboor at West. Initialize arrays.
                size1_send = mesh.numberOfNodesInEachDir[0];
                size2_send = mesh.numberOfNodesInEachDir[2];
                size1_recv = size1_send;
                size2_recv = size2_send;
                ElectricNodes_toSend = new Node3DField*[size1_send];
                ElectricNodes_toRecv = new Node3DField*[size1_recv];
                for(size_t I = 0 ; I < size1_send ; I ++)
                    ElectricNodes_toSend[I] = new Node3DField[size2_send];
                for(size_t I = 0 ; I < size1_recv ; I ++)
                    ElectricNodes_toRecv[I] = new Node3DField[size2_recv];
            }else{
                counterReq += 2;
            }
        }
        if(omp_get_thread_num() == 3){
            // Communication with EAST.
            direction = 'E';
            if(isNeighboor_at_EAST){
                counterReq += 2;
                checkRequest[6] = 1;
                checkRequest[7] = 1;
                // There is a neighboor at EAST. Initialize arrays.
                size1_send = mesh.numberOfNodesInEachDir[0];
                size2_send = mesh.numberOfNodesInEachDir[2];
                size1_recv = size1_send;
                size2_recv = size2_send;
                ElectricNodes_toSend = new Node3DField*[size1_send];
                ElectricNodes_toRecv = new Node3DField*[size1_recv];
                for(size_t I = 0 ; I < size1_send ; I ++)
                    ElectricNodes_toSend[I] = new Node3DField[size2_send];
                for(size_t I = 0 ; I < size1_recv ; I ++)
                    ElectricNodes_toRecv[I] = new Node3DField[size2_recv];
            }else{
                counterReq += 2;
            }
        }
        if(omp_get_thread_num() == 4){
            // Communication with Down.
            direction = 'D';
            if(isNeighboor_at_DOWN){
                counterReq += 2;
                checkRequest[8] = 1;
                checkRequest[9] = 1;
                // There is a neighboor at DOWN. Initialize arrays.
                size1_send = mesh.numberOfNodesInEachDir[0];
                size2_send = mesh.numberOfNodesInEachDir[1];
                size1_recv = size1_send;
                size2_recv = size2_send;
                ElectricNodes_toSend = new Node3DField*[size1_send];
                ElectricNodes_toRecv = new Node3DField*[size1_recv];
                for(size_t I = 0 ; I < size1_send ; I ++)
                    ElectricNodes_toSend[I] = new Node3DField[size2_send];
                for(size_t I = 0 ; I < size1_recv ; I ++)
                    ElectricNodes_toRecv[I] = new Node3DField[size2_recv];
            }else{
                counterReq += 2;
            }
        }
        if(omp_get_thread_num() == 5){
            // Communication with UP.
            direction = 'U';
            if(isNeighboor_at_UP){
                counterReq += 2;
                checkRequest[10] = 1;
                checkRequest[11] = 1;
                // There is a neighboor at UP. Initialize arrays.
                size1_send = mesh.numberOfNodesInEachDir[0];
                size2_send = mesh.numberOfNodesInEachDir[1];
                size1_recv = size1_send;
                size2_recv = size2_send;
                ElectricNodes_toSend = new Node3DField*[size1_send];
                ElectricNodes_toRecv = new Node3DField*[size1_recv];
                for(size_t I = 0 ; I < size1_send ; I ++)
                    ElectricNodes_toSend[I] = new Node3DField[size2_send];
                for(size_t I = 0 ; I < size1_recv ; I ++)
                    ElectricNodes_toRecv[I] = new Node3DField[size2_recv];
            }else{
                counterReq += 2;
            }
        }
        */
        /////////////////////////////////
        /// ALL OMP THREADS WAIT HERE ///
        /////////////////////////////////
        /*
        #pragma omp barrier
        if(REQ_MPI != counterReq){
            printf("AlgoElectro.cpp::ERROR\n");
            printf("\tThe number of MPI requests should be %d but has %d.\n",
                REQ_MPI,counterReq);
            printf("\tAborting.\n");
            std::abort();
        }
        */
       /* printf("Inside algoElectro update :: aborting()\n");
        abort();*/
        int counter = 0;
        while(t_current<t_final){

            /* AFFICHER LA TRANCHE Z=3 */
            int tranche = 3;
            #pragma omp master
            {
                printf("\n\t***************AFFICHAGE DE LA TRANCHE Z = %d***********\n\n",tranche);
                printf("ELECTYRIC FIELD\n");
                for(unsigned long i=1; i <= mesh.numberOfNodesInEachDir[0] ;i++){
                        for(unsigned long j=1; j <= mesh.numberOfNodesInEachDir[1]  ;j++){
                            if(mesh.nodesElec(i,j,tranche).field[2] == 0){
                               // printf("%f,",mesh.nodesElec(i,j,tranche).field[2]);
                               #ifdef _WIN32
                                    printf("%f,",mesh.nodesElec(i,j,tranche).field[2]);
                                #elif __linux__
                                    printf("%s%f%s,",KGRN,mesh.nodesElec(i,j,tranche).field[2],KNRM);
                                #endif
                            }else if(mesh.nodesElec(i,j,tranche).field[2] < 0){
                                #ifdef _WIN32
                                    printf("%f,",mesh.nodesElec(i,j,tranche).field[2]);
                                 #elif __linux__
                                        printf("%s%f%s,",KBLU,mesh.nodesElec(i,j,tranche).field[2],KNRM);
                                 #endif
                            }else{
                                #ifdef _WIN32
                                    printf("%f,",mesh.nodesElec(i,j,tranche).field[2]);
                                #elif __linux__
                                    printf("%s%f%s,",KRED,mesh.nodesElec(i,j,tranche).field[2],KNRM);
                                #endif
                            }
                        }
                        printf("\n");
                }
                printf("Magn field\n");
                for(unsigned long i=0; i <= mesh.numberOfNodesInEachDir[0]+1 ;i++){
                        for(unsigned long j=0; j <= mesh.numberOfNodesInEachDir[1]+1  ;j++){
                            if(mesh.nodesMagn(i,j,tranche).field[2] == 0){
                               #ifdef _WIN32
                                    printf("%f,",mesh.nodesMagn(i,j,tranche).field[2]);
                                #elif __linux__
                                    printf("%s%f%s,",KGRN,mesh.nodesMagn(i,j,tranche).field[2],KNRM);
                                #endif
                            }else if(mesh.nodesMagn(i,j,tranche).field[2] < 0){
                                #ifdef _WIN32
                                    printf("%f,",mesh.nodesMagn(i,j,tranche).field[2]);
                                 #elif __linux__
                                        printf("%s%f%s,",KBLU,mesh.nodesMagn(i,j,tranche).field[2],KNRM);
                                 #endif
                            }else{
                                #ifdef _WIN32
                                    printf("%f,",mesh.nodesMagn(i,j,tranche).field[2]);
                                #elif __linux__
                                    printf("%s%f%s,",KRED,mesh.nodesMagn(i,j,tranche).field[2],KNRM);
                                #endif
                            }
                        }
                        printf("\n");
                }
                //char temp;
                //cin >> temp;  // pour appuyer sur enter à chaque affichage de tableau
            }
            #pragma omp barrier
            

            ///////////////////////////
            /* UPDATE MAGNETIC FIELD */
            ///////////////////////////
            for(unsigned long k = 0 ; k <= mesh.numberOfNodesInEachDir[2] ; k++ ){

                for(unsigned long j=0;j <= mesh.numberOfNodesInEachDir[1] ; j++){

                    for(unsigned long i=0;i <= mesh.numberOfNodesInEachDir[0]  ; i++){
                        printf("%s\t>>> NODE MAGN(%ld,%ld,%ld) <<<%s\n",KGRN,i,j,k,KNRM);
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
                        
                        /* COEFFICIENTS MAGNETIC FIELD */
                        double COEF_H  = elec_conduct_mat*dt/(2*mu_material);

                        double C_hxh   = (1-COEF_H) / (1+COEF_H);
                        double C_hxe_1 = 1 / ( 1 + COEF_H) * dt/(mu_material*deltaZ);
                        double C_hxe_2 = 1 / ( 1 + COEF_H) * dt/(mu_material*deltaY);
                        
                        double C_hyh   = (1-COEF_H) / (1+COEF_H);
                        double C_hye_1 = 1 / ( 1 + COEF_H) * dt/(mu_material*deltaX);
                        double C_hye_2 = 1 / ( 1 + COEF_H) * dt/(mu_material*deltaZ);
                        
                        double C_hzh   = (1-COEF_H) / (1+COEF_H);
                        double C_hze_1 = 1 / ( 1 + COEF_H) * dt/(mu_material*deltaY); 
                        double C_hze_2 = 1 / ( 1 + COEF_H) * dt/(mu_material*deltaX);

                        /* mettre dans gridcreator "numberofProcess" et " myrank"  */
                        //Transfo = mesh.LocalToGlobal(i,j,k,mesh.MPI_communicator.getRank(),mesh.numberofprocess) ;  
                        printf("> isInsideSource...");
                        mesh.input_parser.source.isInsideSource(global[0],global[1],global[2]);
                        printf("Done.\n");

                        if(false && mesh.input_parser.source.isInsideSource(global[0],global[1],global[2])){                
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
                            

                            double temp = C_hxh * mesh.nodesMagn(i,j,k).field[0]
                                                            + C_hxe_1 * (mesh.nodesElec(i,j,k+1).field[1]-
                                                                            mesh.nodesElec(i,j,k).field[1]) 
                                                            - C_hxe_2 * (mesh.nodesElec(i,j+1,k).field[2]-
                                                                            mesh.nodesElec(i,j,k).field[2]);
                                                
                            mesh.nodesMagn(i,j,k).field[0] = (((1-elec_conduct_mat*dt/(2*mu_material)))/((1+elec_conduct_mat*dt/(2*mu_material))))*mesh.nodesMagn(i,j,k).field[0]  +

                                            1/((1+elec_conduct_mat*dt/(2*mu_material)))* ((dt/(mu_material*mesh.deltaZ))*(mesh.nodesElec(i,j,k+1).field[1]
                                            -mesh.nodesElec(i,j,k).field[1]) -

                                            (dt/(mu_material*mesh.deltaY))*(mesh.nodesElec(i,j+1,k).field[2]
                                            -mesh.nodesElec(i,j,k).field[2]));

                            if(abs(mesh.nodesMagn(i,j,k).field[0] - temp) > 1E-15){
                                cout << "THERE ARE DIFFERENT[0]" << temp << "!=" << mesh.nodesMagn(i,j,k).field[0] << endl;
                                abort();
                            }

                            temp = C_hyh * mesh.nodesMagn(i,j,k).field[1]
                                    + C_hye_1 * (mesh.nodesElec(i+1,j,k).field[2]-
                                                 mesh.nodesElec(i,j,k).field[2])
                                    - C_hye_2 * (mesh.nodesElec(i,j,k+1).field[0]-
                                                 mesh.nodesElec(i,j,k).field[0]);

                            /* update magnetic Field H_y */
                            mesh.nodesMagn(i,j,k).field[1]= (((1-elec_conduct_mat*dt/(2*mu_material)))/((1+elec_conduct_mat*dt/(2*mu_material))))*mesh.nodesMagn(i,j,k).field[1]  +
                                            1/((1+elec_conduct_mat*dt/(2*mu_material)))* ((dt/(mu_material*mesh.deltaX))*(mesh.nodesElec(i+1,j,k).field[2]
                                            -mesh.nodesElec(i,j,k).field[2]) -

                                            (dt/(mu_material*mesh.deltaZ))*(mesh.nodesElec(i,j,k+1).field[0]
                                            -mesh.nodesElec(i,j,k).field[0]));
                            
                            if(abs(mesh.nodesMagn(i,j,k).field[1] - temp) > 1E-15){
                                cout << "THERE ARE DIFFERENT[1]" << temp << "!=" << mesh.nodesMagn(i,j,k).field[1] << endl;
                                abort();
                            }
                            
                            temp = C_hzh * mesh.nodesMagn(i,j,k).field[2]
                                    + C_hze_1 * (mesh.nodesElec(i,j+1,k).field[0]-
                                                 mesh.nodesElec(i,j,k).field[0])
                                    - C_hze_2 * (mesh.nodesElec(i+1,j,k).field[1]-
                                                 mesh.nodesElec(i,j,k).field[1]);

                            /* update magnetic Field H_z */
                            mesh.nodesMagn(i,j,k).field[2]= (((1-elec_conduct_mat*dt/(2*mu_material)))/((1+elec_conduct_mat*dt/(2*mu_material))))*mesh.nodesMagn(i,j,k).field[2]  +
                                            1/((1+elec_conduct_mat*dt/(2*mu_material)))* ((dt/(mu_material*mesh.deltaY))*(mesh.nodesElec(i,j+1,k).field[0]-
                                            mesh.nodesElec(i,j,k).field[0]) -

                                            (dt/(mu_material*mesh.deltaX))*(mesh.nodesElec(i+1,j,k).field[1]-
                                            mesh.nodesElec(i,j,k).field[1]));

                            if(abs(mesh.nodesMagn(i,j,k).field[2] - temp) > 1E-15){
                                cout << "THERE ARE DIFFERENT[2]" << temp << "!=" << mesh.nodesMagn(i,j,k).field[2] << endl;
                                abort();
                            }
                        }
                    }
                }
            }
            
            
            ///////////////////////////
            /* UPDATE ELECTRIC FIELD */
            ///////////////////////////
            printf("%sELECTRIC FIELD%s\n",KRED,KNRM);
            for(unsigned long k=1 ; k <= mesh.numberOfNodesInEachDir[0];k++){
                for(unsigned long j=1 ; j <= mesh.numberOfNodesInEachDir[1];j++){
                    for(unsigned long i=1; i <= mesh.numberOfNodesInEachDir[2];i++){
                        printf("%s\t>>> NODE ELEC(%ld,%ld,%ld) <<<%s\n",KGRN,i,j,k,KNRM);
                        /* Get the global indices */
                        local[0] = i;
                        local[1] = j;
                        local[2] = k;
                        mesh.LocalToGlobal(local,global);

                        T = mesh.nodesElec(i,j,k).Temperature;

                        epsilon_material = mesh.materials.getProperty(T,mesh.nodesElec(i,j,k).material,5);

                        elec_conduct_mat = mesh.materials.getProperty(T,mesh.nodesElec(i,j,k).material,6);

                        /* COEFFICIENTS ELECTRIC FIELD */
                        double COEF_E  = elec_conduct_mat*dt/(2*epsilon_material);
                        
                        double C_exe   = (1-COEF_E) / (1+COEF_E);
                        double C_exh_1 = 1 / ( 1 + COEF_E) * dt/(epsilon_material*deltaY);
                        double C_exh_2 = 1 / ( 1 + COEF_E) * dt/(epsilon_material*deltaZ);

                        double C_eye   = (1-COEF_E) / (1+COEF_E);
                        double C_eyh_1 = 1 / ( 1 + COEF_E) * dt/(epsilon_material*deltaZ);
                        double C_eyh_2 = 1 / ( 1 + COEF_E) * dt/(epsilon_material*deltaX);

                        double C_eze   = (1-COEF_E) / (1+COEF_E);
                        double C_ezh_1 = 1 / ( 1 + COEF_E) * dt/(epsilon_material*deltaX);
                        double C_ezh_2 = 1 / ( 1 + COEF_E) * dt/(epsilon_material*deltaY);


                        if(mesh.input_parser.source.isInsideSource(global[0],global[1],global[2])){

                            mesh.input_parser.source.computeSourceValue(mesh, global[0],global[1],global[2],t_current,'E');
                            //mesh.nodesMagn(i,j,k).field[1]=mesh.input_parser.source.computeSourceValue(mesh, i,j,k,t_current,composants_5);
                            //mesh.nodesMagn(i,j,k).field[2]=mesh.input_parser.source.computeSourceValue(mesh, i,j,k,t_current,composants_6);
                        }else{
                            /* update electric field E_x */
                            double temp = C_exe * mesh.nodesElec(i,j,k).field[0]
                                        + C_exh_1 * (mesh.nodesMagn(i,j,k).field[2]-
                                                     mesh.nodesMagn(i,j-1,k).field[2])
                                        - C_exh_2 * (mesh.nodesMagn(i,j,k).field[1]-
                                                     mesh.nodesMagn(i,j,k-1).field[1]);
                            mesh.nodesElec(i,j,k).field[0]= (((1-elec_conduct_mat*dt/(2*epsilon_material)))/((1+elec_conduct_mat*dt/(2*epsilon_material))))*mesh.nodesElec(i,j,k).field[0] +
                            
                                                 1/((1+elec_conduct_mat*dt/(2*epsilon_material)))* ((dt/(epsilon_material*mesh.deltaY))*(mesh.nodesMagn(i,j,k).field[2]-
                                                 mesh.nodesMagn(i,j-1,k).field[2]) - 

                                                 (dt/(epsilon_material*mesh.deltaZ))*(mesh.nodesMagn(i,j,k).field[1]-
                                                 mesh.nodesMagn(i,j,k-1).field[1]));

                            if(abs(temp - mesh.nodesElec(i,j,k).field[0]) > 1E-15){
                                cout << "They are different nodesElec[0]" << endl;
                                abort();
                            }

                            /* update electric Field E_y */
                            temp = C_eye * mesh.nodesElec(i,j,k).field[1]
                                    + C_eyh_1 * (mesh.nodesMagn(i,j,k).field[0]-
                                                 mesh.nodesMagn(i,j,k-1).field[0])
                                    - C_eyh_2 * (mesh.nodesMagn(i,j,k).field[2]-
                                                 mesh.nodesMagn(i-1,j,k).field[2]);

                            double TEMP_1 = (((1-elec_conduct_mat*dt/(2*epsilon_material)))/((1+elec_conduct_mat*dt/(2*epsilon_material))));
                            
                            mesh.nodesElec(i,j,k).field[1]= TEMP_1*mesh.nodesElec(i,j,k).field[1] + 
                                                1/((1+elec_conduct_mat*dt/(2*epsilon_material)))*( ((dt/(epsilon_material*mesh.deltaZ))*(mesh.nodesMagn(i,j,k).field[0]-
                                                mesh.nodesMagn(i,j,k-1).field[0])  - 

                                                (dt/(epsilon_material*mesh.deltaX))*(mesh.nodesMagn(i,j,k).field[2]-
                                                mesh.nodesMagn(i-1,j,k).field[2])));

                            if(abs(temp - mesh.nodesElec(i,j,k).field[1]) > 1E-15){
                                cout << "They are different nodesElec[1]" << endl;
                                cout << temp << "!=" << mesh.nodesElec(i,j,k).field[1] << endl;
                                abort();
                            }
                                
                            /* update electric Field E_z */
                            temp = C_eze * mesh.nodesElec(i,j,k).field[2]
                                    + C_ezh_1 * (mesh.nodesMagn(i,j,k).field[1]-
                                                 mesh.nodesMagn(i-1,j,k).field[1])
                                    - C_ezh_2 * (mesh.nodesMagn(i,j,k).field[0]-
                                                 mesh.nodesMagn(i,j-1,k).field[0]);


                            mesh.nodesElec(i,j,k).field[2]= (((1-elec_conduct_mat*dt/(2*epsilon_material)))/((1+elec_conduct_mat*dt/(2*epsilon_material))))*mesh.nodesElec(i,j,k).field[2] +
                                                 1/((1+elec_conduct_mat*dt/(2*epsilon_material)))* ((dt/(epsilon_material*mesh.deltaX))*(mesh.nodesMagn(i,j,k).field[1]-
                                                 mesh.nodesMagn(i-1,j,k).field[1]) -
                                                 
                                                  (dt/(epsilon_material*mesh.deltaY))*(mesh.nodesMagn(i,j,k).field[0]-
                                                  mesh.nodesMagn(i,j-1,k).field[0]));

                            if(abs(temp - mesh.nodesElec(i,j,k).field[2]) > 1E-15){
                                cout << "They are different nodesElec[2]" << endl;
                                cout << temp << "!=" << mesh.nodesElec(i,j,k).field[2] << endl;
                                abort();
                            }
                        }
                    }
                }
            }   
            
            /////////////////////////////////////////
            /* COMMUNICATION BETWEEN MPI PROCESSES */
            /////////////////////////////////////////
            /*
            if(omp_get_thread_num() >= 0 && omp_get_thread_num() <= 5){
                // Request the electric field to send:
                mesh.GetVecSend(direction,&ElectricNodes_toSend,
                    size1_send,size2_send);
                mesh.communicateWithNeighboor(direction,
                    size1_send,size2_send,
                    size1_recv,size2_recv,
                    &ElectricNodes_toSend,
                    &ElectricNodes_toRecv,
                    &requests_MPI);
            }
            */
            ///////////////////////////////////////////////
            /* CHECK ALL THE MPI REQUESTS ARE FULLFILLED */
            ///////////////////////////////////////////////
            /*
            #pragma omp master
            {
                for(int iii = 0 ; iii < REQ_MPI ; iii++){
                    if(checkRequest[iii] == 1){
                        MPI_Wait(&requests_MPI[iii],MPI_STATUS_IGNORE);
                    }
                }
            }
            #pragma omp barrier
            */

            ////////////////////////////////////////////////
            /* PUT RECEIVED FIELDS INSIDE THE RIGHT ARRAY */
            ////////////////////////////////////////////////
            /*
            if(omp_get_thread_num() >=0 && omp_get_thread_num() <=5){
                mesh.SetVecRecv(direction,&ElectricNodes_toRecv,
                    size1_recv,size2_recv);
            }
            #pragma omp barrier
            */
            /* WRITE RESULTS - ONLY MAIN OMP THREAD */
            ////////////////////////////////////////
            /// UPDATING CURRENT SIMULATION TIME ///
            /// ONLY MAIN OMP THREAD !!!         ///
            ////////////////////////////////////////
            /*
            #pragma omp master
            {
                */
                interfaceForOutput.convertAndWriteData(this->currentStep);
                
                t_current=t_current+dt;
                this->currentStep ++;
                /*
            }
            #pragma omp barrier
            */
            //MPI_Barrier(MPI_COMM_WORLD);
        }

        printf("ALGO::END::DT=%f::FINAL_TIME=%f\n",dt,t_current);

        //////////////////////////////////////
        /// FREE MEMORY OF EACH OMP THREAD ///
        //////////////////////////////////////
        /*
        if(ElectricNodes_toRecv != NULL){
            for(size_t I = 0 ; I < size1_recv ; I ++)
                delete[] ElectricNodes_toRecv[I];
            delete[] ElectricNodes_toRecv;
        }
        if(ElectricNodes_toSend != NULL){
            for(size_t I = 0 ; I < size1_send ; I ++)
                delete[] ElectricNodes_toSend[I];
            delete[] ElectricNodes_toSend;
        }
        */
    /*
    }
    */

    //// FREEING MEMORY ////
    /*free(requests_MPI);*/
}


AlgoElectro::~AlgoElectro(void){
    cout << "Calling AlgoCreator destructor\n";
}


