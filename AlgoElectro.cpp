#include "Materials.h"
#include "GridCreator.h"
#include <fstream>
#include <cstring>


#define dt mesh.deltaT
#define dx mesh.deltaX
#define dy mesh.deltaY
#define dz mesh.deltaZ
#define mu mesh.materials.getproperties()   /* !!!!!!!!!à finir*/
#define epsilon mesh.materials.getpr         /* !!!!! à finir */



void communicate(GridCreator){


}

/*Table array contient le maillage du MPI+ les voisins*/
/*nbrElts ne contient pas les voisins*/
/* le maillage correspond à celui du champ electrique*/

void udpate(GridCreator mesh){
    int i,j,k;


    /* to make */
    void communicate(GridCreator mesh)
    /* update magnetic Field */
    /* start to i,j,k=2 pas prendre le voisin en compte */

    /* Verifier les indices */
    for (i=1;i<mesh.nbrElts_X;i++){
        for(j=1;j<mesh.nbrElts_Y;j++){
            for(k=1;j<mesh.nbrElts_Z;k++){
        /* update magnetic field H_x */
        mesh.nodeMagn(i,j,k).field[0]= mesh.nodeMagn(i,j,k).field[0] + (dt/(mu*dz))*(mesh.nodeElec(i,j,k).field[1]-mesh.nodeElec(i,j,k-1).field[1]) - (dt/(mu*dy))*(mesh.nodeElec(i,j,k).field[2]-mesh.nodeElec(i,j-1,k).field[2]);
        /* update magnetic Field H_y */
        mesh.nodeMagn(i,j,k).field[1]= mesh.nodeMagn(i,j,k).field[1]  + (dt/(mu*dx))*(mesh.nodeElec(i,j,k).field[2]-mesh.nodeElec(i-1,j,k).field[2]) -  (dt/(mu*dz))*(mesh.nodeElec(i,j,k).field[0]-mesh.nodeElec(i,j,k-1).field[0]);
        /* update magnetic Field H_z */
        mesh.nodeMagn(i,j,k).field[2]= mesh.nodeMagn(i,j,k).field[2]  + (dt/(mu*dy))*(mesh.nodeElec(i,j,k).field[0]-mesh.nodeElec(i,j-1,k).field[0]) -  (dt/(mu*dx))*(mesh.nodeElec(i,j,k).field[1]-mesh.nodeElec(i-1,j,k).field[1]);
            }
        }
    }
    /* update electric field  */
    for (i=2;i<mesh.nbrElts_X;i++){
        for(j=2;j<mesh.nbrElts_Y;j++){
            for(k=2;j<mesh.nbrElts_Z;k++){
        /* update magnetic field H_x */
        mesh.nodeElec(i,j,k).field[0]= mesh.nodeElec(i,j,k).field[0] + (dt/(epsilon*dy))*(mesh.nodeMagn(i,j,k).field[2]-mesh.nodeMagn(i,j-1,k).field[2]) - (dt/(epsilon*dz))*(mesh.nodeMagn(i,j,k).field[1]-mesh.nodeMagn(i,j,k-1).field[1]);
        /* update magnetic Field H_y */
        mesh.nodeElec(i,j,k).field[1]= mesh.nodeElec(i,j,k).field[1] + (dt/(epsilon*dz))*(mesh.nodeMagn(i,j,k).field[0]-mesh.nodeMagn(i,j,k-1).field[0])  - (dt/(epsilon*dx))*(mesh.nodeMagn(i,j,k).field[2]-mesh.nodeMagn(i-1,j,k).field[2];
        /* update magnetic Field H_z */
        mesh.nodeElec(i,j,k).field[2]= mesh.nodeElec(i,j,k).field[2] + (dt/(epsilon*dx))*(mesh.nodeMagn(i,j,k).field[1]-mesh.nodeMagn(i-1,j,k).field[1]) - (dt/(epsilon*dy))*(mesh.nodeMagn(i,j,k).field[0]-mesh.nodeMagn(i,j-1,k).field[0]);
            }
        }
    }



}



void run(GridCreator mesh){
/* rajouter dans classe nodes dx (left and right),dy (left and right),dz(left and right) => 6 parameters to put*/
/* rajouter propriétés des préfacteurs pour calculer E et H */

/* to make  */
void udpate()
}
//this.udpdate*

