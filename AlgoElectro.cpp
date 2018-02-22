#include "AlgoElectro.h"

#include "Materials.h"
#include "GridCreator.h"
#include <fstream>
#include <cstring>
#include "ElectromagneticSource.h"
#include "Node3DField.h"







void AlgoElectro::communicate(GridCreator mesh){
/*à faire*/

}

/* convention:Table array contient le maillage du MPI+ les voisins*/
/* convention: nbrElts ne contient pas les voisins*/


double AlgoElectro::Compute_dt(GridCreator mesh){
    double dx=mesh.deltaX;
    double dy=mesh.deltaY;
    double dz=mesh.deltaZ;
    double dt=0.0;
    double tmp=0.0;
    double c=0.0;
    int i=0;                                                                                
    for (i=0;i<mesh.materials.numberOfMaterials;i++){
            double mu_material=mesh.materials.getProperty(Materials_object[i].T_in,i,4);         /* !!!!!!!!!!!!!!/*� faire avec T initial*/     
            double epsilon_material=mesh.materials.getProperty(Materials_object[i].T_in,i,5);    /* !!!!!!!!!!!!!!/*� faire avec T initial*/
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

void AlgoElectro::update(GridCreator mesh,double dt,double t_current){   
    unsigned long i,j,k;

    double T = 0.0;
    double mu_material = 0.0;

    unsigned long local[3];
    unsigned long global[3];

    for(k = 1 ; k < mesh.numberOfNodesInEachDir[0] ; k++ ){

        for(j=1;j < mesh.numberOfNodesInEachDir[1];j++){

            for(i=1;i < mesh.numberOfNodesInEachDir[2];i++){
                
                /* Get the global indices */
                local[0] = i;
                local[1] = j;
                local[2] = k;
                mesh.LocalToGlobal(&local,&global);

                T = mesh.nodesMagn(i,j,k).Temperature;

                mu_material = mesh.materials.getProperty(T,4,mesh.nodesMagn(i,j,k).material);  

                /* mettre dans gridcreator "numberofProcess" et " myrank"  */
                //Transfo = mesh.LocalToGlobal(i,j,k,mesh.MPI_communicator.getRank(),mesh.numberofprocess) ;  

                if(mesh.input_parser.source.isInsideSource(global[0],global[1],global[2])){                
                    /* Fonction à faire  dans Electromagnetic Source composantes donne l'info si c'est E ou H et x ou y ou z */
                    mesh.input_parser.source.computeSourceValue(mesh,i,j,k,t_current,'H');
                    /*mesh.nodesMagn(i,j,k).field[1]=mesh.input_parser.source.computeSourceValue(mesh,i,j,k,t_current);
                    mesh.nodesMagn(i,j,k).field[2]=mesh.input_parser.source.computeSourceValue(mesh,i,j,k,t_current);*/
                }
                else{

                    /* update magnetic field H_x */
                    mesh.nodesMagn(i,j,k).field[0]= mesh.nodesMagn(i,j,k).field[0]  +
                                     (dt/(mu_material*mesh.deltaZ))*(mesh.nodesElec(i,j,k).field[1]
                                     -mesh.nodesElec(i,j,k-1).field[1]) -
                                     (dt/(mu_material*mesh.deltaY))*(mesh.nodesElec(i,j,k).field[2]
                                     -mesh.nodesElec(i,j-1,k).field[2]);
                    /* update magnetic Field H_y */
                    mesh.nodesMagn(i,j,k).field[1]= mesh.nodesMagn(i,j,k).field[1]  +
                                     (dt/(mu_material*mesh.deltaX))*(mesh.nodesElec(i,j,k).field[2]-mesh.nodesElec(i-1,j,k).field[2]) -
                                     (dt/(mu_material*mesh.deltaZ))*(mesh.nodesElec(i,j,k).field[0]-mesh.nodesElec(i,j,k-1).field[0]);
                    /* update magnetic Field H_z */
                    mesh.nodesMagn(i,j,k).field[2]= mesh.nodesMagn(i,j,k).field[2]  +
                                     (dt/(mu_material*mesh.deltaY))*(mesh.nodesElec(i,j,k).field[0]-mesh.nodesElec(i,j-1,k).field[0]) -
                                     (dt/(mu_material*mesh.deltaX))*(mesh.nodesElec(i,j,k).field[1]-mesh.nodesElec(i-1,j,k).field[1]);
                }
            }
        }
    }

    /* update electric field  */
    double epsilon_material = 0.0;
    for(k=1;k<mesh.numberOfNodesInEachDir[0];k++){
        for(j=1;j<mesh.numberOfNodesInEachDir[1];j++){
            for(i=1;i<mesh.numberOfNodesInEachDir[2];i++){
                
                /* Get the global indices */
                local[0] = i;
                local[1] = j;
                local[2] = k;
                mesh.LocalToGlobal(&local,&global);

                T = mesh.nodesElec(i,j,k).Temperature;

                epsilon_material = mesh.materials.getProperty(T,4,mesh.nodesElec(i,j,k).material);

                if(mesh.input_parser.source.isInsideSource(global[0],global[1],global[2])){
                    mesh.input_parser.source.computeSourceValue(mesh, i,j,k,t_current,'E');
                    //mesh.nodesMagn(i,j,k).field[1]=mesh.input_parser.source.computeSourceValue(mesh, i,j,k,t_current,composants_5);
                    //mesh.nodesMagn(i,j,k).field[2]=mesh.input_parser.source.computeSourceValue(mesh, i,j,k,t_current,composants_6);
                }else{
                    /* update magnetic field E_x */
                    mesh.nodesElec(i,j,k).field[0]= mesh.nodesElec(i,j,k).field[0] + (dt/(epsilon_material*mesh.deltaY))*(mesh.nodesMagn(i,j+1,k).field[2]-mesh.nodesMagn(i,j,k).field[2]) - (dt/(epsilon_material*mesh.deltaZ))*(mesh.nodesMagn(i,j,k+1).field[1]-mesh.nodesMagn(i,j,k).field[1]);
                    /* update magnetic Field E_y */
                    mesh.nodesElec(i,j,k).field[1]= mesh.nodesElec(i,j,k).field[1] + (dt/(epsilon_material*mesh.deltaZ))*(mesh.nodesMagn(i,j,k+1).field[0]-mesh.nodesMagn(i,j,k).field[0])  - (dt/(epsilon_material*mesh.deltaX))*(mesh.nodesMagn(i+1,j,k).field[2]-mesh.nodesMagn(i,j,k).field[2]);
                    /* update magnetic Field E_z */
                    mesh.nodesElec(i,j,k).field[2]= mesh.nodesElec(i,j,k).field[2] + (dt/(epsilon_material*mesh.deltaX))*(mesh.nodesMagn(i+1,j,k).field[1]-mesh.nodesMagn(i,j,k).field[1]) - (dt/(epsilon_material*mesh.deltaY))*(mesh.nodesMagn(i,j+1,k).field[0]-mesh.nodesMagn(i,j,k).field[0]);
                }
            }
        }
    }   
}


/*Fonction principale*/
void AlgoElectro::run(GridCreator mesh){
    double t_current=0.0;
    /*rajouter le temps final dans le fichier*/
    double t_final = mesh.input_parser.get_stopTime();                                      
    /*calcul condition CFL*/
    double dt=Compute_dt(mesh);
    /*loop over time*/
    while(t_current<t_final){        
        this->communicate(mesh);                              /* à faire*/
        this->update(mesh,dt,t_current);
        t_current=t_current+dt;
    }
}




