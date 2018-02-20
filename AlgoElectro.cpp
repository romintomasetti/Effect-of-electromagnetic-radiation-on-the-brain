#include "Materials.h"
#include "GridCreator.h"
#include <fstream>
#include <cstring>
#include "ElectromagneticSource.h"
#include "Node3DField.h"


#define dx GridCreator.deltaX
#define dy GridCreator.deltaY
#define dz GridCreator.deltaY
#define frequency ElectromagneticSource.getfrequency()

void communicate(GridCreator){


}

/*Table array contient le maillage du MPI+ les voisins*/
/*nbrElts ne contient pas les voisins*/
/* le maillage correspond à celui du champ electrique*/


double Compute_dt(double dx,double dy, double dz){
    /*probleme comment avoir les différents matériaux */
double dt=0;
double tmp=0;
double c=0;
for (i=0;i<Materials.numberOfMaterials;i++){
        double mu_material=Materials.getProperty(Materials[i].T_in,i,4);              /* !!!!!!!!!!!!!!/*à faire avec T initial*/
        double epsilon_material=Materials.getProperty(Materials[i].T_in,i,5);
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

void update(GridCreator mesh){
    int i,j,k;

    double dt=Compute_dt(dx,dy,dz);

    double t_final=mesh.InputParser.tempsfinal; /*rajouter le temps final*/

    double t_current=0;

    /*loop over time*/
    while(t_current<t_final){
        /*to make*/
        void communicate(GridCreator mesh)
        /* update magnetic Field */
        /* start to i,j,k=2 pas prendre le voisin en compte */
        for (i=1;i<mesh.nbrElts_X;i++){
            for(j=1;j<mesh.nbrElts_Y;j++){
                for(k=1;j<mesh.nbrElts_Z;k++){
                    if(isinsightsource(i,j,k)){
                        mesh.nodeMagn(i,j,k).field[0]=sin(frequency*2*pi*t)
                        mesh.nodeMagn(i,j,k).field[1]=
                        mesh.nodeMagn(i,j,k).field[2]=
                    }
                    else{

                        /* update magnetic field H_x */
                        mesh.nodeMagn(i,j,k).field[0]= mesh.nodeMagn(i,j,k).field[0] + (dt/(mu*dz))*(mesh.nodeElec(i,j,k).field[1]-mesh.nodeElec(i,j,k-1).field[1]) - (dt/(mu*dy))*(mesh.nodeElec(i,j,k).field[2]-mesh.nodeElec(i,j-1,k).field[2]);
                        /* update magnetic Field H_y */
                        mesh.nodeMagn(i,j,k).field[1]= mesh.nodeMagn(i,j,k).field[1]  + (dt/(mu*dx))*(mesh.nodeElec(i,j,k).field[2]-mesh.nodeElec(i-1,j,k).field[2]) -  (dt/(mu*dz))*(mesh.nodeElec(i,j,k).field[0]-mesh.nodeElec(i,j,k-1).field[0]);
                        /* update magnetic Field H_z */
                        mesh.nodeMagn(i,j,k).field[2]= mesh.nodeMagn(i,j,k).field[2]  + (dt/(mu*dy))*(mesh.nodeElec(i,j,k).field[0]-mesh.nodeElec(i,j-1,k).field[0]) -  (dt/(mu*dx))*(mesh.nodeElec(i,j,k).field[1]-mesh.nodeElec(i-1,j,k).field[1]);

                        }
                    }
                }
            }
        /* update electric field  */
            for (i=2;i<mesh.nbrElts_X;i++){
                for(j=2;j<mesh.nbrElts_Y;j++){
                    for(k=2;j<mesh.nbrElts_Z;k++){
                    /* update magnetic field H_x */
                    mesh.nodeElec(i,j,k).field[0]= mesh.nodeElec(i,j,k).field[0] + (dt/(epsilon*dy))*(mesh.nodeMagn(i,j+1,k).field[2]-mesh.nodeMagn(i,j,k).field[2]) - (dt/(epsilon*dz))*(mesh.nodeMagn(i,j,k+1).field[1]-mesh.nodeMagn(i,j,k).field[1]);
                    /* update magnetic Field H_y */
                    mesh.nodeElec(i,j,k).field[1]= mesh.nodeElec(i,j,k).field[1] + (dt/(epsilon*dz))*(mesh.nodeMagn(i,j,k+1).field[0]-mesh.nodeMagn(i,j,k).field[0])  - (dt/(epsilon*dx))*(mesh.nodeMagn(i+1,j,k).field[2]-mesh.nodeMagn(i,j,k).field[2];
                    /* update magnetic Field H_z */
                    mesh.nodeElec(i,j,k).field[2]= mesh.nodeElec(i,j,k).field[2] + (dt/(epsilon*dx))*(mesh.nodeMagn(i+1,j,k).field[1]-mesh.nodeMagn(i,j,k).field[1]) - (dt/(epsilon*dy))*(mesh.nodeMagn(i,j+1,k).field[0]-mesh.nodeMagn(i,j,k).field[0]);
                        }
                    }
                }
            t_current=t_current+dt;
    }
}



void run(GridCreator mesh){
/* rajouter dans classe nodes dx (left and right),dy (left and right),dz(left and right) => 6 parameters to put*/
/* rajouter propriétés des préfacteurs pour calculer E et H */

/* to make  */
void update()
}
//this.udpdate*

