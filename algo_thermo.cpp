#include "algo_thermo.hpp"

using namespace vtl;


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Function get_my_rank
int get_my_rank()
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    return myid;
}

//Function check_MUMPS
void check_MUMPS(DMUMPS_STRUC_C &id)
{
    if (id.infog[0] < 0)
    {
        std::cout << "[" << get_my_rank() << "] MUMPS Error:\n";
        std::cout << "\tINFOG(1)=" << id.infog[0] << '\n';
        std::cout << "\tINFOG(2)=" << id.infog[1] << std::endl;
    }
}



// Function init_MUMPS
void init_MUMPS(DMUMPS_STRUC_C &id)
{
    id.comm_fortran = -987654; //USE_COMM_WORLD;
    id.par = 1;                // 1=host involved in factorization phase
    id.sym = 0;                // 0=unsymmetric
    id.job = -1;
    std::cout << "[" << get_my_rank() << "] Init MUMPS package." << std::endl;
    dmumps_c(&id);
    check_MUMPS(id);
}



//Function end_MUMPS
void end_MUMPS(DMUMPS_STRUC_C &id)
{
    id.job = -2;
    std::cout << "[" << get_my_rank() << "] Terminate MUMPS instance." << std::endl;
    dmumps_c(&id);
    check_MUMPS(id);
}

//Function Wall
void wall_geometry(unsigned int *material_at_nodes,unsigned int N_x,unsigned int N_y, unsigned int N_z){
    for(unsigned int parcours=0 ; parcours < N_x*N_y*N_z ; parcours++){
        unsigned int tmp = floor(parcours/(N_x*N_y));
        if(tmp<(unsigned int) (0.5*(double)N_z)){
            material_at_nodes[parcours]=0;
        }else{
            material_at_nodes[parcours]=1;
        }
    
    }
}
//Function création geometry
void remplissage_material_at_nodes(unsigned int *material_at_nodes,unsigned int N_x,unsigned int N_y,unsigned int N_z){
    for(unsigned int parcours=0; parcours < N_x*N_y*N_z ; parcours++ ){
        material_at_nodes[parcours]=0;
    }
}



//Function convection around brain
void convection_brain(double *A, MUMPS_INT *Indices_line_A, MUMPS_INT *Indices_row_A,unsigned int *counter_nonvalue_A,unsigned int *position_equation,unsigned  int parcours,double Delta,unsigned int N_x,unsigned int N_y,unsigned int N_z,double h,double T_infiny,double *k,double *convection_contribution,unsigned int dir,unsigned int wichmaterial){
    
    //Face x avant ""i==1"
    if (dir==0){
        A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+2;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours;     
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =h;
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;     
        (*counter_nonvalue_A)++;
        convection_contribution[*position_equation]=h*T_infiny;
        (*position_equation)++;
    }

    //Face x après
    if(dir==1){
        A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+2;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours;     
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =h;
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;     
        (*counter_nonvalue_A)++;
        convection_contribution[*position_equation]=h*T_infiny;
        (*position_equation)++;

    }

    //Face y gauche
    if(dir==2){
        A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =h;
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;     
        (*counter_nonvalue_A)++;
        convection_contribution[*position_equation]=h*T_infiny;
        (*position_equation)++;
    }

    //Face y droite
    if(dir==3){
        A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =h;
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;     
        (*counter_nonvalue_A)++;
        convection_contribution[*position_equation]=h*T_infiny;
        (*position_equation)++;
    }
    //Face z  bas 
    if(dir==4){
       A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x*N_y;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x*N_y;     
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =h;
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;     
        (*counter_nonvalue_A)++;
        convection_contribution[*position_equation]=h*T_infiny;
        (*position_equation)++; 
    }
    //Face z haut 
    if(dir==5){
        A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x*N_y;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x*N_y;     
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =h;
        Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;     
        (*counter_nonvalue_A)++;
        convection_contribution[*position_equation]=h*T_infiny;
        (*position_equation)++;
    }
}


//conditions limits fonction
void condition_limit_equation(double *B, MUMPS_INT *Indices_line_B, MUMPS_INT *Indices_row_B,double *A, MUMPS_INT *Indices_line_A, MUMPS_INT *Indices_row_A,unsigned int *counter_nonvalue_B,unsigned int *counter_nonvalue_A,unsigned int *position_equation,unsigned int parcours,unsigned int *Stateofeachface,double Delta,unsigned int N_x,unsigned int N_y,unsigned int N_z,unsigned int *neighbooroutside,double h,double T_infiny,double *k,double *convection_contribution,unsigned int wichmaterial)
{


    //recoit le tableau k => savoir dans quelle matérieau on est 

    // Face i==1
    if(neighbooroutside[0]==1){

        if(Stateofeachface[0]==0){

            B[*counter_nonvalue_B] = 1;
            Indices_line_B[*counter_nonvalue_B] = (*position_equation)+1;
            Indices_row_B[*counter_nonvalue_B] = parcours;
            A[*counter_nonvalue_A] =1;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours;     
            (*counter_nonvalue_A)++;
            (*counter_nonvalue_B)++;
            (*position_equation)++;

        }else if(Stateofeachface[0]==1){
            
            A[*counter_nonvalue_A] =-1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+2;     
            (*counter_nonvalue_A)++;
            (*position_equation)++;
        

        }else if(Stateofeachface[0]==2){
            A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+2;
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =h;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1;     
            (*counter_nonvalue_A)++;
            convection_contribution[*position_equation]=h*T_infiny;
            (*position_equation)++;
        
        }else{
            printf("Problem for Face ???, Line %d\n",__LINE__);
            abort();
        }
    }

    // Face i==N_x-2
    if(neighbooroutside[1]==1){

        if(Stateofeachface[1]==0){

            B[*counter_nonvalue_B] = 1;
            Indices_line_B[*counter_nonvalue_B] = (*position_equation)+1;
            Indices_row_B[*counter_nonvalue_B] = parcours+2;
            A[*counter_nonvalue_A] =1;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+2;     
            (*counter_nonvalue_A)++;
            (*counter_nonvalue_B)++;
            (*position_equation)++;

        }else if(Stateofeachface[1]==1){

            A[*counter_nonvalue_A] =1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+2;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =-1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours;     
            (*counter_nonvalue_A)++;
            (*position_equation)++;
        
        }else if(Stateofeachface[1]==2){
            A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+2;
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =h;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1;     
            (*counter_nonvalue_A)++;
            convection_contribution[*position_equation]=h*T_infiny;
            (*position_equation)++;
        
        }else{
            printf("Problem for Face ???, Line %d\n",__LINE__);
            abort();
        }
    }
    
    //Face j==1
    if(neighbooroutside[2]==1){

        if(Stateofeachface[2]==0){

            B[*counter_nonvalue_B] = 1;
            Indices_line_B[*counter_nonvalue_B] = (*position_equation)+1;
            Indices_row_B[*counter_nonvalue_B] = parcours+1-N_x;
            A[*counter_nonvalue_A] =1;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
            (*counter_nonvalue_A)++;
            (*counter_nonvalue_B)++;
            (*position_equation)++;

        }else if(Stateofeachface[2]==1){

            A[*counter_nonvalue_A] =1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =-1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;     
            (*counter_nonvalue_A)++;
            (*position_equation)++;
        

        }else if(Stateofeachface[2]==2){
            A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =h;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1;     
            (*counter_nonvalue_A)++;
            convection_contribution[*position_equation]=h*T_infiny;
            (*position_equation)++;
            
        }else{
            printf("Problem for Face ????, Line %d\n",__LINE__);
            abort();
        }
    }

    //Face j==N_y-2
    if(neighbooroutside[3]==1){

        if(Stateofeachface[3]==0){

            B[*counter_nonvalue_B] = 1;
            Indices_line_B[*counter_nonvalue_B] = (*position_equation)+1;
            Indices_row_B[*counter_nonvalue_B] = parcours+1+N_x;
            A[*counter_nonvalue_A] =1;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;     
            (*counter_nonvalue_A)++;
            (*counter_nonvalue_B)++; 
            (*position_equation)++;

        }else if(Stateofeachface[3]==1){
            
            A[*counter_nonvalue_A] =1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =-1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;     
            (*counter_nonvalue_A)++;
            (*position_equation)++;

        }else if(Stateofeachface[3]==2){
            A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =h;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1;     
            (*counter_nonvalue_A)++;
            convection_contribution[*position_equation]=h*T_infiny;
            (*position_equation)++;
          
        }else{
            printf("Problem for Face ???, Line %d\n",__LINE__);
            abort();
        }
    }
    //Face k==1
    if(neighbooroutside[4]==1){

        if(Stateofeachface[4]==0){
            B[*counter_nonvalue_B] = 1;
            Indices_line_B[*counter_nonvalue_B] =(*position_equation)+1;
            Indices_row_B[*counter_nonvalue_B] = parcours+1-(N_x*N_y);
            A[*counter_nonvalue_A] =1;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-(N_x*N_y);     
            (*counter_nonvalue_A)++;
            (*counter_nonvalue_B)++; 
            (*position_equation)++;

        }else if(Stateofeachface[4]==1){

            A[*counter_nonvalue_A] =1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x*N_y;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =-1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x*N_y;     
            (*counter_nonvalue_A)++;
            (*position_equation)++;

        }else if(Stateofeachface[4]==2){
            A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x*N_y;
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x*N_y;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =h;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1;     
            (*counter_nonvalue_A)++;
            convection_contribution[*position_equation]=h*T_infiny;
            (*position_equation)++;
            
        }else{
            printf("Problem for Face ????, Line %d\n",__LINE__);
            abort();
        }
    }
    //Face k==N_z-2
    if(neighbooroutside[5]==1){
        if(Stateofeachface[5]==0){

            B[*counter_nonvalue_B] = 1;
            Indices_line_B[*counter_nonvalue_B] = (*position_equation)+1;
            Indices_row_B[*counter_nonvalue_B] = parcours+1+(N_x*N_y);
            A[*counter_nonvalue_A] =1;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+(N_x*N_y);     
            (*counter_nonvalue_A)++;
            (*counter_nonvalue_B)++; 
            (*position_equation)++;

        }else if(Stateofeachface[5]==1){

            A[*counter_nonvalue_A] =1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x*N_y;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =-1/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x*N_y;     
            (*counter_nonvalue_A)++;
            (*position_equation)++;

        }else if(Stateofeachface[5]==2){
            A[*counter_nonvalue_A] =k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x*N_y;
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =-k[wichmaterial]/(2*Delta);
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x*N_y;     
            (*counter_nonvalue_A)++;
            A[*counter_nonvalue_A] =h;
            Indices_line_A[*counter_nonvalue_A]=(*position_equation)+1;
            Indices_row_A[*counter_nonvalue_A]=parcours+1;     
            (*counter_nonvalue_A)++;
            convection_contribution[*position_equation]=h*T_infiny;
            (*position_equation)++;

        }else{
            printf("Problem for Face ???, Line %d\n",__LINE__);
            abort();
        }
    }

}















//Function to give the number of non zero in matrice A and B
 void Numberofnon_zero_function(MKL_INT *numberofnon_nullvalue_A,MKL_INT *numberofnon_nullvalue_B, unsigned int *Stateofeachface,unsigned int N_x,unsigned int N_y,unsigned int N_z){
   
    // corners
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x*2+(N_y-2)*2)*2+(N_z-2)*4;
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x*2+(N_y-2)*2)*2+(N_z-2)*4;
    
    // Face 1 i==0
    if(Stateofeachface[0]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_y-2)*(N_z-2);
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_y-2)*(N_z-2);
    }else if(Stateofeachface[0]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_y-2)*(N_z-2)*2;
    }else if(Stateofeachface[0]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_y-2)*(N_z-2)*3;       
    }
    // Face 2 i==N_x-1
    if(Stateofeachface[1]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_y-2)*(N_z-2);
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_y-2)*(N_z-2);
    }else if(Stateofeachface[1]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_y-2)*(N_z-2)*2;
    }else if(Stateofeachface[1]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_y-2)*(N_z-2)*3;
    }

    // Face 3 j==0
    if(Stateofeachface[2]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_z-2);
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x-2)*(N_z-2);
    }else if(Stateofeachface[2]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_z-2)*2;
    }else if(Stateofeachface[2]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_z-2)*3;
    }

    // Face 4 j==N_y-1
    if(Stateofeachface[3]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_z-2);
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x-2)*(N_z-2);
    }else if(Stateofeachface[3]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_z-2)*2;
    }else if(Stateofeachface[3]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_z-2)*3;
    }

    // Face 5 k==0
    if(Stateofeachface[4]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2);
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x-2)*(N_y-2);
    }else if(Stateofeachface[4]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2)*2;
    }else if(Stateofeachface[4]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2)*3;
    }

    // Face 6 k==N_z-1
    if(Stateofeachface[5]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2);
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x-2)*(N_y-2);
    }else if(Stateofeachface[5]==1){        
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2)*2;
     }else if(Stateofeachface[5]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2)*3;
    }
    
    
    // Middle Contribution
    unsigned int contribution_middle=7*(N_x-2)*(N_y-2)*(N_z-2);
    
    *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+contribution_middle;
    *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+contribution_middle;
 }








//Function to give the Temperature initial
void InitializeTemperature(double *InitialTemperature,unsigned int Number_total,unsigned int N_x,unsigned int N_y,unsigned int N_z, double T_infiny, double Delta,double *temperature_initial,unsigned int *material_at_nodes, double *temperature_initial_face,unsigned int *Stateofeachface, unsigned int type_simulation_value, unsigned int thermal_distribution,double  amplitude_thermal_distribution){

        unsigned int count_y=0;
        
        for(unsigned int parcours=0;parcours<Number_total;parcours++){
            //Thermal distibution
            if(thermal_distribution==1){
                InitialTemperature[parcours]=amplitude_thermal_distribution*sin(M_PI*(double)count_y/(double)(N_x-1));
            }else{
                //Analytic case
                if(type_simulation_value==0){
                    unsigned int tmp = floor(parcours/(N_x*N_y));
                    unsigned int value_testx= parcours-tmp*N_x*N_y;
                    if(count_y==0){
                        if(Stateofeachface[0]==0){
                            InitialTemperature[parcours]=temperature_initial_face[0];
                        }
                        else {
                            InitialTemperature[parcours]=temperature_initial[material_at_nodes[parcours]];
                        }
                    }else if(count_y==N_x-1){
                        if(Stateofeachface[1]==0){
                            InitialTemperature[parcours]=temperature_initial_face[1];
                        }else{
                            InitialTemperature[parcours]=temperature_initial[material_at_nodes[parcours]];
                        }
                    }else if(value_testx < N_x){
                        if(Stateofeachface[2]==0){
                            InitialTemperature[parcours]=temperature_initial_face[2];
                        }else{
                            InitialTemperature[parcours]=temperature_initial[material_at_nodes[parcours]];
                        }
                    }else if(value_testx >= (N_x*N_y)-N_x){
                        if(Stateofeachface[3]==0){
                            InitialTemperature[parcours]=temperature_initial_face[3];
                        }else{
                            InitialTemperature[parcours]=temperature_initial[material_at_nodes[parcours]];
                        }
                    }else if(parcours<N_x*N_y){
                        if(Stateofeachface[4]==0){
                            InitialTemperature[parcours]=temperature_initial_face[4];
                        }else{
                            InitialTemperature[parcours]=temperature_initial[material_at_nodes[parcours]];
                        }
                    }else if(parcours >= Number_total -N_x*N_y){
                        if(Stateofeachface[5]==0){
                            InitialTemperature[parcours]=temperature_initial_face[5];
                            printf("InitialTemperature[%d]=%lf",parcours,InitialTemperature[parcours]);
                            printf("hello\n");
                            //abort();
                        }else{
                            InitialTemperature[parcours]=temperature_initial[material_at_nodes[parcours]];
                        }
                    }else{
                        InitialTemperature[parcours]=temperature_initial[material_at_nodes[parcours]];
                    }

                //Brain case
                }else if(type_simulation_value==1){
                    InitialTemperature[parcours]=temperature_initial[material_at_nodes[parcours]];
                //Problemd
                }else{
                    printf("Type_simulation_value is not 0 or 1.This error comes from Line %d \n",__LINE__);
                    abort();
                }
            }
        count_y=count_y+1;
        if(count_y==N_x){
            count_y = 0;
        }
	    }                 
}


//Function to read the file containing the heat source  in each node
// Rmodifier après !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void ReadFile(int Number_total, double *Q,double dt){
	
	int i;
	FILE *fichier = NULL;
    fichier = fopen("Value_Q_thermique.txt","r");

	if(fichier != NULL){
        double tempp;

        for(i=0; i<Number_total; i++){
            fscanf(fichier,"%lf",&tempp);
            Q[i] = tempp*dt;
        }

        fclose(fichier);

	}else{
        printf("FDP ton file n'est pas fopen. Cette erreur vient de la ligne %d \n",__LINE__);
		abort();
    }
}


//Function to compute B*Temperature_before+Q
void daxpy_call(int Number_total,double *Q,double *BY_results){
    double coeff_saxpy =  1;

    cblas_daxpy(Number_total,
	            coeff_saxpy, 
	            Q,
	            1,
	            BY_results,
	            1);
}



//Function to compute B*Temperature_before
void mkl_call(int Number_total,double *B, int *Indices_line_B, int *Colonn_line_B, int numberofnon_nullvalue_B, double *Temperature_before,double *BY_results){
   char DoNotTranspose = 'N';
    mkl_dcoogemv(&DoNotTranspose, // Parce-que la matrice B n'est pas transposée
		        &Number_total, // La taille de B
		        B, // La matrice B
		        Indices_line_B, // L'index des lignes qui contiennent des valeurs
		        Colonn_line_B, // L'index des colonnes qui contiennent des valeurs
		        &numberofnon_nullvalue_B, // Le nombre d'élément(s) non nul(s)
		        Temperature_before, // Le vecteur x
		        BY_results); // Le product B*Y
}



//Function resolve
void resolve(DMUMPS_STRUC_C &id, vtl::SPoints &grid,unsigned int Number_total, unsigned int *Stateofeachface, unsigned int N_x, unsigned int N_y, unsigned int N_z, double Delta, double dt,double theta,double h,double T_infiny, unsigned int type_simulation_value, double *temperature_initial,double t_final_thermal, double * temperature_initial_face, unsigned int rate_save_thermo,char *geometry_material, unsigned int wall_thermo, unsigned int thermal_distribution,double  amplitude_thermal_distribution,unsigned int heat_distribution, double amplitude_heat_distribution,GridCreator_NEW & gridElectro)
{
    // Read the file and put it inside a vector:
    unsigned int* material_at_nodes = NULL;
    //Brain
    if(type_simulation_value==1){
        size_t size_mat = 0;
        
        if( NULL != (material_at_nodes = read_input_geometry_file(geometry_material,&size_mat))){
            printf("%s\n",geometry_material);
            if(size_mat != N_x * N_y * N_z){
                printf("size_mat %zu\n",size_mat);
                printf("Number total =%u\n",N_x*N_y*N_z);
                printf("Size error!! %d",__LINE__);
                abort();
            }
            printf("size_mat =%zu\n",size_mat);
            printf("Number total =%u\n",N_x*N_y*N_z);
        }else{        
            printf("There is an abort() at line %d. Sorry for that.\n",__LINE__);
            abort();
        }
        // Case analytic
    }else{
        material_at_nodes = (unsigned int *) calloc(Number_total,sizeof(unsigned int));
        if(wall_thermo==1){
            wall_geometry(material_at_nodes,N_x,N_y,N_z);
        }else{            
            remplissage_material_at_nodes(material_at_nodes,N_x,N_y,N_z);
        }

    }

    
   
   


    
    //counter y
    unsigned int count_y = 0;
    unsigned int *nodes_nearbrain_prim =(unsigned int *) calloc(Number_total,sizeof(unsigned int));
    if(nodes_nearbrain_prim == NULL){	 
        printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
        abort();
    }

    //Number of non zero values in matrix A
  	MKL_INT  numberofnon_nullvalue_A=0;
    MKL_INT  numberofnon_nullvalue_B=0;
    unsigned int nodes_surface=0;

    printf("calcul du nombre de zero pour A et B \n");

    if(type_simulation_value==0){
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analytic case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Numberofnon_zero_function(&numberofnon_nullvalue_A,&numberofnon_nullvalue_B, Stateofeachface, N_x, N_y, N_z);
    }else{
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case of the brain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for (unsigned int parcours=0; parcours<Number_total; parcours++){
            //variable to known the position
            unsigned int tmp = floor(parcours/(N_x*N_y));
            unsigned int value_testx = parcours-tmp*N_x*N_y;           
            unsigned int b=0;  

            // In the air
            if(material_at_nodes[parcours]==0){
                //Face i==0
                if(count_y==0 && nodes_nearbrain_prim[parcours]==0){
                    nodes_nearbrain_prim[parcours]=0;
                    b=1;
                }
                //Face i==N_x-1
                if(count_y==N_x-1 && nodes_nearbrain_prim[parcours]==0){
                    nodes_nearbrain_prim[parcours]=0;
                    b=1;
                }
                //Face j==0
                if(value_testx < N_x && nodes_nearbrain_prim[parcours]==0){
                    nodes_nearbrain_prim[parcours]=0;
                    b=1;
                }
                //Face j==N_y-1
                if(value_testx >= (N_x*N_y)-N_x && nodes_nearbrain_prim[parcours]==0){
                    nodes_nearbrain_prim[parcours]=0;
                    b=1;
                }
                //Face z==0
                if(parcours <= N_x*N_y && nodes_nearbrain_prim[parcours]==0){
                    nodes_nearbrain_prim[parcours]=0;
                    b=1;  
                }
                //Face z==N_z-1
                if(parcours >= Number_total -N_x*N_y && nodes_nearbrain_prim[parcours]==0){
                   nodes_nearbrain_prim[parcours]=0;
                    b=1; 
                }
                if(b==0){
                    // neighbour at x+1 if no air ? 
                    if(material_at_nodes[parcours+1]!=0 && nodes_nearbrain_prim[parcours]==0 && b==0){
                        nodes_nearbrain_prim[parcours]=1;
                        nodes_surface++;
                    }
                    // neighbour at x-1 if no air ? 
                    
                    if(material_at_nodes[parcours-1]!=0 && nodes_nearbrain_prim[parcours]==0 && b==0){
                        nodes_nearbrain_prim[parcours]=1;
                        nodes_surface++;
                    }
                    // neighbour at y+1 if no air ? 
                    if(material_at_nodes[parcours+N_x]!=0 && nodes_nearbrain_prim[parcours]==0 && b==0){
                        nodes_nearbrain_prim[parcours]=1;
                        nodes_surface++;
                    }
                    // neighbour at y-1 if no air ? 
                    if(material_at_nodes[parcours-N_x]!=0 && nodes_nearbrain_prim[parcours]==0 && b==0){
                        nodes_nearbrain_prim[parcours]=1;
                        nodes_surface++;
                    }
                    // neighbour at z+1 if no air ? 
                    if(material_at_nodes[parcours+N_x*N_y]!=0 && nodes_nearbrain_prim[parcours]==0 && b==0){
                        nodes_nearbrain_prim[parcours]=1;
                        nodes_surface++;
                    }
                    // neighbour at z-1 if no air ? 
                    if(material_at_nodes[parcours-N_x*N_y]!=0 && nodes_nearbrain_prim[parcours]==0 && b==0){
                        nodes_nearbrain_prim[parcours]=1;
                        nodes_surface++;
                    }
                }
            } 
            // Position along the x-axis
            count_y=count_y+1;
            if(count_y==N_x){
                    count_y = 0;
            }
        }
        //Reset
        nodes_surface=0;

        printf("1er step for matrix A and B est fais\n");
        //abort();

        unsigned int tmp_numb=0;
        unsigned int tmp_numb2=0;
        unsigned int count=0;
        // Determin the number of no-zero values for matrix A and B 
        for(unsigned int parcours=0; parcours<Number_total; parcours++){
            if(nodes_nearbrain_prim[parcours]==1){
                tmp_numb2++;
            }
            count++;


            if(nodes_nearbrain_prim[parcours]==1){
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_A++;
                tmp_numb++;
            }
            if(nodes_nearbrain_prim[parcours]==0  && material_at_nodes[parcours]==0){
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_B++;
                tmp_numb++;
            }
            if(material_at_nodes[parcours]!=0){
                 //update
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_B++;
                        
                //selon x

                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_B++;
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_B++;

                    //selony
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_B++;
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_B++;

                //selonz
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_B++;
                numberofnon_nullvalue_A++;
                numberofnon_nullvalue_B++;
                tmp_numb++;
            }
        }
        
    printf("tmp_numb=%u\n",tmp_numb);
    printf("count=%u\n",count);
    printf("nodes_surface:%u\n",nodes_surface);
    printf("tmp_numb2=%u\n",tmp_numb2);
    }
    printf("Fin du calcul du nombre de zero pour A et B \n");

    //%%%%%%%%%%%%%%%%%%%%%%%% For the 2 cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // Variable for the matrix A , B and convection
    // Format MKL ici pour indices!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double *B = (double *) calloc(numberofnon_nullvalue_B,sizeof(double));    
    MUMPS_INT *Indices_line_B = (MUMPS_INT *) calloc(numberofnon_nullvalue_B,sizeof(MUMPS_INT));
    MUMPS_INT*Indices_row_B = (MUMPS_INT *) calloc(numberofnon_nullvalue_B,sizeof(MUMPS_INT));
    double *A = (double *) calloc(numberofnon_nullvalue_A,sizeof(double));    
    MUMPS_INT *Indices_line_A = (MUMPS_INT *) calloc(numberofnon_nullvalue_A,sizeof(MUMPS_INT));
    MUMPS_INT *Indices_row_A = (MUMPS_INT *) calloc(numberofnon_nullvalue_A,sizeof(MUMPS_INT));
    double *convection_contribution = (double *) calloc(Number_total,sizeof(double));
    
    int step=0;


    // Verification of the calloc
     if(B == NULL || Indices_line_B == NULL || Indices_row_B == NULL || A == NULL || Indices_line_A == NULL || Indices_row_A == NULL || convection_contribution == NULL)
    {
        printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
        abort();
    }
    

    // count the number of non-value
    unsigned int counter_nonvalue_A=0;
    unsigned int counter_nonvalue_B=0;




    
    // Vector of the heat source
    double *Q = (double *) calloc(Number_total,sizeof(double));
    if(Q == NULL){
        printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
        abort();
    }


    double Cst_A;
    double Cst_B;
    unsigned int position_equation=0;

    //nodes at the extern
    unsigned int *nodes_extern =(unsigned int *) calloc(Number_total,sizeof(unsigned int));
    // nodes around the brain (nodes in air)
    unsigned int *nodes_nearbrain = (unsigned int*) calloc(Number_total,sizeof(unsigned int));
    unsigned int *nodes_normal_1= (unsigned int*) calloc(Number_total,sizeof(unsigned int));
    unsigned int *nodes_normal_2= (unsigned int*) calloc(Number_total,sizeof(unsigned int));
    unsigned int *nodes_normal_3= (unsigned int*) calloc(Number_total,sizeof(unsigned int));
    unsigned int *nodes_normal_4= (unsigned int*) calloc(Number_total,sizeof(unsigned int));
    unsigned int *nodes_normal_5= (unsigned int*) calloc(Number_total,sizeof(unsigned int));
    unsigned int *nodes_normal_6= (unsigned int*) calloc(Number_total,sizeof(unsigned int));

    //Verification calloc

    if(nodes_extern == NULL || nodes_nearbrain == NULL || nodes_normal_1 == NULL || nodes_normal_2 == NULL || nodes_normal_3 == NULL || nodes_normal_4 == NULL || nodes_normal_5 == NULL || nodes_normal_6 == NULL ){
        printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
        abort();
    }



    unsigned int tmp_nodes_normales=0;
    unsigned int tmp_nodes_normales_exter=0;
    unsigned int countxxx=0;

    //Variable for Paraview
    std::vector<double> test(Number_total);
    if(get_my_rank() == 0)
    {
        //%%%%%%%%%%%%%%%%%%%%  Nodes externes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        count_y=0;
        for (unsigned int parcours=0; parcours<Number_total; parcours++){
            
            //Variables use for located
            unsigned int tmp = floor(parcours/(N_x*N_y));
            unsigned int value_testx = parcours-tmp*N_x*N_y;

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case analytique %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            if(type_simulation_value==0){
                if( parcours<N_x*N_y && count_y==0){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //x=N_x
                else if( parcours<N_x*N_y && count_y==N_x-1){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //y=0
                else if(parcours<N_x*N_y && value_testx < N_x){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //y=N_y
                else if(parcours<N_x*N_y && value_testx >= (N_x*N_y)-N_x){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //Plan z=N_z
                //x=0
                else if(parcours >= Number_total -N_x*N_y && count_y==0){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //x=N_x
                else if(parcours >= Number_total -N_x*N_y && count_y==N_x-1){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //y=0
                else if(parcours >= Number_total -N_x*N_y && value_testx < N_x){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //y=N_y
                else if(parcours >= Number_total -N_x*N_y && value_testx >= (N_x*N_y)-N_x){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //corner sur verticalement
                //x==0 et y==0
                else if(count_y==0 && value_testx < N_x){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }

                //x==N_x et y==0
                else if(count_y==N_x-1 && value_testx < N_x){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }

                //x==0 et y==N_y
                else if(count_y==0 && value_testx >= (N_x*N_y)-N_x){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }
                //x==N_x et y==N_y
                else if(count_y==N_x-1 && value_testx >= (N_x*N_y)-N_x){
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    nodes_extern[parcours]=1;
                    position_equation++;
                }  

                //Nodes extern
                //Face 1  i==0
                else if(count_y == 0){
                    nodes_extern[parcours]=1;
                    
                }

                //Nodes extern
                //Face 2   i==N_x-1
                else if (count_y == N_x-1){
                    nodes_extern[parcours]=1;
                    
                }
                
                //nodes extern
                //Face 3 j==0
                else if (value_testx < N_x){
                    nodes_extern[parcours]=1;
                    
                }

                //nodes extern
                //Face 4  j==N_y-1
                else if (value_testx >= (N_x*N_y)-N_x){
                    nodes_extern[parcours]=1;
                    
                }


                //nodes extern
                //Face 5  k==0
                else if(parcours  < N_x*N_y){
                    nodes_extern[parcours]=1;
                    
                }

                //nodes extern
                //Face 6 k==N_z-1
                else if(parcours >= Number_total -N_x*N_y){
                    nodes_extern[parcours]=1;
                    
                }
            }

            //%%%%%%%%%%%%%%%%%%%%%%%%% Case brain boundary of the brain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(type_simulation_value==1){
                unsigned int a=0;
                if(material_at_nodes[parcours]==0){

                    if(nodes_nearbrain_prim[parcours]==1){
                        countxxx++;
                        if(material_at_nodes[parcours+1]>0){
                            tmp_nodes_normales_exter++;
                            nodes_nearbrain[parcours]=1;
                            if(material_at_nodes[parcours+2]>0){
                                nodes_normal_2[parcours]=1;
                                tmp_nodes_normales++;
                                a=1;
                            }
                        }
                        if(material_at_nodes[parcours-1]>0 && a==0){
                            nodes_nearbrain[parcours]=1;
                            tmp_nodes_normales_exter++;
                            if(material_at_nodes[parcours-2]>0){
                                nodes_normal_1[parcours]=1;
                                tmp_nodes_normales++;
                                a=1;   
                            }
                        }
                        if(material_at_nodes[parcours+N_x]>0 && a==0){
                            nodes_nearbrain[parcours]=1;
                            tmp_nodes_normales_exter++;
                            if(material_at_nodes[parcours+2*N_x]>0){
                                nodes_normal_4[parcours]=1;
                                tmp_nodes_normales++; 
                                a=1;
                            }
                        }
                        if(material_at_nodes[parcours-N_x]>0 && a==0){
                            nodes_nearbrain[parcours]=1;
                            tmp_nodes_normales_exter++;
                            if(material_at_nodes[parcours-2*N_x]>0){
                                nodes_normal_3[parcours]=1;
                                tmp_nodes_normales++;
                                a=1;
                            }
                        }
                        if(material_at_nodes[parcours+N_x*N_y]>0 && a==0){
                            nodes_nearbrain[parcours]=1;
                            tmp_nodes_normales_exter++;
                            if(material_at_nodes[parcours+2*N_x*N_y]>0){
                                nodes_normal_6[parcours]=1;
                                tmp_nodes_normales++;
                                a=1;
                            }
                        }
                        if(material_at_nodes[parcours-N_x*N_y]>0 && a==0){
                            nodes_nearbrain[parcours]=1;
                            tmp_nodes_normales_exter++;
                            if(material_at_nodes[parcours-2*N_x*N_y]>0){
                                nodes_normal_5[parcours]=1;
                                tmp_nodes_normales++;
                                a=1;
                            }
                        }
                        if(a==0){
                            printf("Un noeud trouve pas pour faire la convection!!!!!\n");
                            printf("numero noeud %d\n",parcours);
                            abort();
                        }
                    }

                }

            }
            count_y=count_y+1;
            if(count_y==N_x){
                    count_y = 0;
            }
        }

        //Verification
        unsigned int counter_nodes=0;
        for(unsigned int parcours=0; parcours < Number_total ; parcours++){
            if(nodes_nearbrain[parcours]==1){
                counter_nodes++;
            }
        }
        printf("counter_nodes=%u\n",counter_nodes);
        printf("tmp_nodes_normales=%u\n",tmp_nodes_normales);
        printf("tmp_nodes_normales_exter=%u\n",tmp_nodes_normales_exter);
        printf("countxxx=%u",countxxx);        
        //abort();

        //%%%%%%%%%%%%%%%%%%%%  Nodes insides   %%%%%%%%%%%%%%%%%
        unsigned int c=0;
        count_y=0;
        unsigned int nodes_inter_count=0;
        for (unsigned int parcours=0; parcours < Number_total; parcours++){
            unsigned int dir;
            unsigned int wichmaterial;
            unsigned int tmp = floor(parcours/(N_x*N_y));
            unsigned int value_testx = parcours-tmp*N_x*N_y;
            unsigned int neighbooroutside[6]={0,0,0,0,0,0};


            //%%%%%%%%%%%%%%%%%%%%% Case Brain %%%%%%%%%%%%%%%%%%%%%%%%ù
            if(type_simulation_value==1){

                //Value paramater
                double *k=(double *) calloc(gridElectro.materials.unified_material_list.size(),sizeof(double));
                double *c_p=(double *)calloc(gridElectro.materials.unified_material_list.size(),sizeof(double));
                double *rho=(double *)calloc(gridElectro.materials.unified_material_list.size(),sizeof(double));
                for(unsigned int i=0; i<gridElectro.materials.unified_material_list.size();i++){
                    k[i]=gridElectro.materials.unified_material_list[i].properties["THERMALCONDUCTIVITY(W/M/°C)"];
                    c_p[i]=gridElectro.materials.unified_material_list[i].properties["HEATCAPACITY(J/KG/°C)"];
                    rho[i]=gridElectro.materials.unified_material_list[i].properties["DENSITY(KG/M³)"];
                    /*printf("> Materials %s : rho= %lf, c_p = %lf , k=%lf\n"
                            ,gridElectro.materials.unified_material_list[i].name.c_str(),gridElectro.materials.unified_material_list[i].properties["DENSITY(KG/M³)"]
                            ,gridElectro.materials.unified_material_list[i].properties["HEATCAPACITY(J/KG/°C)"],
                            gridElectro.materials.unified_material_list[i].properties["THERMALCONDUCTIVITY(W/M/°C)"]);*/
                }

                //neighbour air near the brain 
                if(nodes_nearbrain_prim[parcours]==1){
                    //convection en avant 
                    if(nodes_normal_2[parcours]==1){                        
                        dir=0;
                        wichmaterial=material_at_nodes[parcours+1];
                        convection_brain( A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, &position_equation, parcours+1, Delta, N_x, N_y, N_z, h, T_infiny, k,convection_contribution,dir,wichmaterial);
                        c++;
                    //convection en arrière
                    }   
                    if(nodes_normal_1[parcours]==1){                        
                        dir=1;
                        wichmaterial=material_at_nodes[parcours-1];
                        convection_brain( A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, &position_equation, parcours-1, Delta, N_x, N_y, N_z, h, T_infiny, k,convection_contribution,dir,wichmaterial);
                        c++;
                    //convection à droite
                    }
                    if(nodes_normal_3[parcours]==1){                        
                        dir=3;
                        wichmaterial=material_at_nodes[parcours-N_x];
                        convection_brain( A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, &position_equation, parcours-N_x, Delta, N_x, N_y, N_z, h, T_infiny, k,convection_contribution,dir,wichmaterial);
                        c++;
                    //convection à gauche
                    }
                    if(nodes_normal_4[parcours]==1){                        
                        dir=2;
                        wichmaterial=material_at_nodes[parcours+N_x];
                        convection_brain( A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, &position_equation, parcours+N_x, Delta, N_x, N_y, N_z, h, T_infiny, k,convection_contribution,dir,wichmaterial);
                        c++;
                    }
                    //convection en haut 
                    if(nodes_normal_5[parcours]==1){
                        dir=5;
                        wichmaterial=material_at_nodes[parcours-(N_x*N_y)];
                        convection_brain( A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, &position_equation, parcours-(N_x*N_y), Delta, N_x, N_y, N_z, h, T_infiny, k,convection_contribution,dir,wichmaterial);
                        c++;

                    }
                    //convection en bas
                    if(nodes_normal_6[parcours]==1){
                        dir=4;
                        wichmaterial=material_at_nodes[parcours+(N_x*N_y)];
                        convection_brain( A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, &position_equation, parcours+(N_x*N_y), Delta, N_x, N_y, N_z, h, T_infiny, k,convection_contribution,dir,wichmaterial);
                        c++;
                    }
                    
                }else if(material_at_nodes[parcours]==0){
                    //In air not near neighbour brain 
                    A[counter_nonvalue_A]=1;
                    B[counter_nonvalue_B]=1;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;                    
                    position_equation++;

                }else if(material_at_nodes[parcours]!=0){
                    // In brain 


                    nodes_inter_count++;
                    //update
                    B[counter_nonvalue_B] = (rho[material_at_nodes[parcours]]*c_p[material_at_nodes[parcours]])-(theta*dt/(2*Delta*Delta))*(6*k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]]
                                    +k[material_at_nodes[parcours+1]]+k[material_at_nodes[parcours+N_x]]+k[material_at_nodes[parcours-N_x]]+k[material_at_nodes[parcours+N_x*N_y]]+k[material_at_nodes[parcours-N_x*N_y]]);
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    A[counter_nonvalue_A]=(rho[material_at_nodes[parcours]]*c_p[material_at_nodes[parcours]])+((1-theta)*dt/(2*Delta*Delta))*(6*k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]]
                                    +k[material_at_nodes[parcours+1]]+k[material_at_nodes[parcours+N_x]]+k[material_at_nodes[parcours-N_x]]+k[material_at_nodes[parcours+N_x*N_y]]+k[material_at_nodes[parcours-N_x*N_y]]);
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                            
                            
                    //selon x
                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+1]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+1]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+2;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+2;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                        //selony
                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B]=Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+N_x+1;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+N_x+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours-N_x+1;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours-N_x+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                    //selonz
                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x*N_y]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x*N_y]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+N_x*N_y+1;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+N_x*N_y+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x*N_y]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x*N_y]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] =position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours-N_x*N_y+1;
                    A[counter_nonvalue_A] = Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours-N_x*N_y+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                    
//partie chaleur à remodifier !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Q[position_equation]=1000*sin(M_PI*count_y/(N_x-1))*dt;
                    
                    position_equation++;
                }else{
                    printf("Erreur Remplisage matrice\n");
                    abort();
                }               
            }else{
            //%%%%%%%%%%%%%%%%%%%% Case analytic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                //Parameter
                double rho[2]={0,0};
                double k[2]={0,0};
                double c_p[2]={0,0};
                if(wall_thermo==1){
                    rho[0]=2300;
                    rho[1]=640;
                    k[0]=1.4;
                    k[1]=0.15;
                    c_p[0]=880;
                    c_p[1]=2805;
                }else{
                    rho[0]=3970;
                    k[0]=46;
                    c_p[0]=765;
                }

                if(nodes_extern[parcours]==1){
                        //nothing to do 
                }else{
                    
                    //update
                    B[counter_nonvalue_B] = (rho[material_at_nodes[parcours]]*c_p[material_at_nodes[parcours]])-(theta*dt/(2*Delta*Delta))*(6*k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]]
                                    +k[material_at_nodes[parcours+1]]+k[material_at_nodes[parcours+N_x]]+k[material_at_nodes[parcours-N_x]]+k[material_at_nodes[parcours+N_x*N_y]]+k[material_at_nodes[parcours-N_x*N_y]]);
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+1;
                    A[counter_nonvalue_A]=(rho[material_at_nodes[parcours]]*c_p[material_at_nodes[parcours]])+((1-theta)*dt/(2*Delta*Delta))*(6*k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]]
                                    +k[material_at_nodes[parcours+1]]+k[material_at_nodes[parcours+N_x]]+k[material_at_nodes[parcours-N_x]]+k[material_at_nodes[parcours+N_x*N_y]]+k[material_at_nodes[parcours-N_x*N_y]]);
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
                            
                            
                    //selon x
                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+1]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+1]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+2;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+2;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                        //selony
                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B]=Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+N_x+1;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+N_x+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours-N_x+1;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours-N_x+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                    //selonz
                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x*N_y]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x*N_y]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] = position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours+N_x*N_y+1;
                    A[counter_nonvalue_A]=Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours+N_x*N_y+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;

                    Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x*N_y]])*dt/(2*Delta*Delta);
                    Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x*N_y]])*dt/(2*Delta*Delta);
                    B[counter_nonvalue_B] = Cst_B;
                    Indices_line_B[counter_nonvalue_B] =position_equation+1;
                    Indices_row_B[counter_nonvalue_B] = parcours-N_x*N_y+1;
                    A[counter_nonvalue_A] = Cst_A;
                    Indices_line_A[counter_nonvalue_A]=position_equation+1;
                    Indices_row_A[counter_nonvalue_A]=parcours-N_x*N_y+1;
                    counter_nonvalue_A++;
                    counter_nonvalue_B++;
            
                    if(heat_distribution==1){                        
                        Q[position_equation]=amplitude_heat_distribution*sin(M_PI*(double)count_y/(double)(N_x-1))*dt;
                    }else{
                        Q[position_equation]=0;
                    }
                    
                        
                    position_equation++;
                    //Conditions frontières domaine

                    wichmaterial=material_at_nodes[parcours];
                    //Face i==0
                    if(count_y == 1){
                        //remplissage du tableau voisin externes
                        if(nodes_extern[parcours-1]==1)
                        neighbooroutside[0]=1;
                        if(nodes_extern[parcours+1]==1)
                        neighbooroutside[1]=1;
                        if(nodes_extern[parcours-N_x]==1)
                        neighbooroutside[2]=1;
                        if(nodes_extern[parcours+N_x]==1)
                        neighbooroutside[3]=1;
                        if(nodes_extern[parcours-N_x*N_y]==1)
                        neighbooroutside[4]=1;
                        if(nodes_extern[parcours+N_x*N_y]==1)
                        neighbooroutside[5]=1;
                        
                        //condition 
                        condition_limit_equation(B,Indices_line_B,Indices_row_B,A,Indices_line_A,Indices_row_A,&counter_nonvalue_B,&counter_nonvalue_A,&position_equation,parcours,Stateofeachface,Delta,N_x,N_y,N_z,neighbooroutside,h,T_infiny,k,convection_contribution, wichmaterial);
                        

                    }
                    
                    
                    //Face 2   i==N_x-1
                    else if (count_y == N_x-2){
                        //remplissage du tableau voisin externes
                        if(nodes_extern[parcours-1]==1)
                        neighbooroutside[0]=1;
                        if(nodes_extern[parcours+1]==1)
                        neighbooroutside[1]=1;
                        if(nodes_extern[parcours-N_x]==1)
                        neighbooroutside[2]=1;
                        if(nodes_extern[parcours+N_x]==1)
                        neighbooroutside[3]=1;
                        if(nodes_extern[parcours-N_x*N_y]==1)
                        neighbooroutside[4]=1;
                        if(nodes_extern[parcours+N_x*N_y]==1)
                        neighbooroutside[5]=1;

                        //conditions
                        condition_limit_equation(B,Indices_line_B,Indices_row_B,A,Indices_line_A,Indices_row_A,&counter_nonvalue_B,&counter_nonvalue_A,&position_equation,parcours,Stateofeachface,Delta,N_x,N_y,N_z,neighbooroutside,h,T_infiny,k,convection_contribution, wichmaterial);
                        

                    }

                    
                    //Face 3 j==0
                    else if (value_testx >= N_x && value_testx < 2*N_x ){

                        //remplissage du tableau voisin externes
                        if(nodes_extern[parcours-1]==1)
                        neighbooroutside[0]=1;
                        if(nodes_extern[parcours+1]==1)
                        neighbooroutside[1]=1;
                        if(nodes_extern[parcours-N_x]==1)
                        neighbooroutside[2]=1;
                        if(nodes_extern[parcours+N_x]==1)
                        neighbooroutside[3]=1;
                        if(nodes_extern[parcours-N_x*N_y]==1)
                        neighbooroutside[4]=1;
                        if(nodes_extern[parcours+N_x*N_y]==1)
                        neighbooroutside[5]=1;

                        //conditions                    
                        condition_limit_equation(B,Indices_line_B,Indices_row_B,A,Indices_line_A,Indices_row_A,&counter_nonvalue_B,&counter_nonvalue_A,&position_equation,parcours,Stateofeachface,Delta,N_x,N_y,N_z,neighbooroutside,h,T_infiny,k,convection_contribution,wichmaterial);
                        
                    }
                    
                    
                    //Face 4  j==N_y-1
                    else if(value_testx >= (N_x*N_y)-(2*N_x)  && value_testx < (N_x*N_y)-N_x){
                        //remplissage du tableau voisin externes
                        if(nodes_extern[parcours-1]==1)
                        neighbooroutside[0]=1;
                        if(nodes_extern[parcours+1]==1)
                        neighbooroutside[1]=1;
                        if(nodes_extern[parcours-N_x]==1)
                        neighbooroutside[2]=1;
                        if(nodes_extern[parcours+N_x]==1)
                        neighbooroutside[3]=1;
                        if(nodes_extern[parcours-N_x*N_y]==1)
                        neighbooroutside[4]=1;
                        if(nodes_extern[parcours+N_x*N_y]==1)
                        neighbooroutside[5]=1;

                        //conditions
                        condition_limit_equation(B,Indices_line_B,Indices_row_B,A,Indices_line_A,Indices_row_A,&counter_nonvalue_B,&counter_nonvalue_A,&position_equation,parcours,Stateofeachface,Delta,N_x,N_y,N_z,neighbooroutside,h,T_infiny,k,convection_contribution,wichmaterial);
                    }

                    //Face 5  k==0
                    else if(parcours  >= N_x*N_y && parcours < 2*N_x*N_y){
                        //remplissage du tableau voisin externes
                        if(nodes_extern[parcours-1]==1)
                        neighbooroutside[0]=1;
                        if(nodes_extern[parcours+1]==1)
                        neighbooroutside[1]=1;
                        if(nodes_extern[parcours-N_x]==1)
                        neighbooroutside[2]=1;
                        if(nodes_extern[parcours+N_x]==1)
                        neighbooroutside[3]=1;
                        if(nodes_extern[parcours-N_x*N_y]==1)
                        neighbooroutside[4]=1;
                        if(nodes_extern[parcours+N_x*N_y]==1)
                        neighbooroutside[5]=1;

                        //conditions                    
                        condition_limit_equation(B,Indices_line_B,Indices_row_B,A,Indices_line_A,Indices_row_A,&counter_nonvalue_B,&counter_nonvalue_A,&position_equation,parcours,Stateofeachface,Delta,N_x,N_y,N_z,neighbooroutside,h,T_infiny,k,convection_contribution,wichmaterial);
                        
                    }
                    
                    
                    //Face 6 k==N_z-1
                    else if(parcours >= Number_total -2*N_x*N_y     && parcours < Number_total -N_x*N_y){
                        //remplissage du tableau voisin externes
                        if(nodes_extern[parcours-1]==1)
                        neighbooroutside[0]=1;
                        if(nodes_extern[parcours+1]==1)
                        neighbooroutside[1]=1;
                        if(nodes_extern[parcours-N_x]==1)
                        neighbooroutside[2]=1;
                        if(nodes_extern[parcours+N_x]==1)
                        neighbooroutside[3]=1;
                        if(nodes_extern[parcours-N_x*N_y]==1)
                        neighbooroutside[4]=1;
                        if(nodes_extern[parcours+N_x*N_y]==1)
                        neighbooroutside[5]=1;

                        //conditions
                    
                        condition_limit_equation(B,Indices_line_B,Indices_row_B,A,Indices_line_A,Indices_row_A,&counter_nonvalue_B,&counter_nonvalue_A,&position_equation,parcours,Stateofeachface,Delta,N_x,N_y,N_z,neighbooroutside,h,T_infiny,k,convection_contribution,wichmaterial);

                    }
                    
                    //NO face, we are in the mesh
                    else{
                        //nothing to do with conditions
                    }
                        
                    
                }
            }    
            count_y=count_y+1;
            if(count_y==N_x){
                 count_y = 0;
            }
        }
        printf("c=%u\n",c);
        printf("counter_nonvalue_B=%u\n",counter_nonvalue_B);
        printf("counter_nonvalue_A=%u\n",counter_nonvalue_A);
        printf("numberofnon_nullvalue_A=%u\n",numberofnon_nullvalue_A);
        printf("numberofnon_nullvalue_B=%u\n",numberofnon_nullvalue_B);
        printf("position_equation=%u\n",position_equation);
        printf("Number total=%u\n",Number_total);
        printf("nodes_inter_count=%u\n",nodes_inter_count);
        //abort();  
        
        

        

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of the filling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù

        // Temperature Initial
        double *InitialTemperature = (double *) calloc(Number_total,sizeof(double));
        if(InitialTemperature == NULL){
            printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
            abort();
        }
        InitializeTemperature(InitialTemperature,Number_total,N_x,N_y,N_z,T_infiny,Delta,temperature_initial,material_at_nodes,temperature_initial_face,Stateofeachface,type_simulation_value,thermal_distribution,amplitude_thermal_distribution);
        

        // Lecture of the file outside
        //ReadFile(Number_total,Q,dt); à remettre!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        // Variable 
        double *BY_results = (double *) calloc(Number_total,sizeof(double));

        if(BY_results == NULL)
        {
            printf("Le vecteur TESTVECTOR n'est pas alloué. Aborting...");
            exit(EXIT_FAILURE);
        }
        // contribution de convection mettre dans Q
        daxpy_call(Number_total,convection_contribution,Q);

        // compute B*T_0
        
        mkl_call(Number_total,B,Indices_line_B,Indices_row_B,numberofnon_nullvalue_B,InitialTemperature,BY_results);
       /* for(i=0;i<Number_eqtosolve;i++){
            printf("%lf\n",BY_results[i]);
        }*/
        //Compute B*T_0+Q
        daxpy_call(Number_total,Q,BY_results);
        /*for(i=0;i<Number_eqtosolve;i++){
            printf("%lf\n",BY_results[i]);
        }*/
        

        id.rhs=BY_results;
        
        id.n=Number_total;

        id.nnz=numberofnon_nullvalue_A;
        
        id.irn=Indices_line_A;
            
        id.jcn=Indices_row_A;    
        
        id.a=A;


        

        for(unsigned int i=0;i<Number_total;i++){
                    test[i]=InitialTemperature[i];
            }
        /* for(i=0;i<Number_total;i++){
                printf("%lf\n",id.rhs[i]);
            }*/
    

        grid.scalars["Temp"] = &test;        
        export_spoints_XML("cas3", step, grid, grid, vtl::Zip::ZIPPED);
    }
    //%%%%%%%%%%%%%% end if get_my_rank() == 0 %%%%%%%%%%%%%%%%%%%%%
    MPI_Barrier(MPI_COMM_WORLD);   	

    
    //Sans-Gridcreator
    double t = 0.0;
     

    //Loop over time
    
    step++;
    step =0;
    
    double *Temperature_temp = (double *) calloc(Number_total,sizeof(double));
    if(Temperature_temp == NULL){
            printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
            abort();
        }


    // Pour vérifier l'augmentation de T°
    double *Temperature_verif = (double *) calloc(Number_total,sizeof(double));
    if(Temperature_verif == NULL){
            printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
            abort();
        }

    InitializeTemperature(Temperature_verif,Number_total,N_x,N_y,N_z,T_infiny,Delta,temperature_initial,material_at_nodes,temperature_initial_face,Stateofeachface,type_simulation_value,thermal_distribution,amplitude_thermal_distribution);
        

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // !!!!!!!!!!!!!!!!! count hard codé !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    unsigned int count =0;
    unsigned int number_nodes_inside=0;
    unsigned int number_increasesof1degree=0;
    while(t < t_final_thermal){       
        if(step==0){
        
            id.job=6;
         
            step++;
        }else{
            id.job=3;
        }
        
        dmumps_c(&id);
        check_MUMPS(id);
        
        
        //save here
        if(get_my_rank() ==0){
            if(count==rate_save_thermo){
                for(unsigned int i=0;i<Number_total;i++){
                        test[i]=id.rhs[i];
                }

                grid.scalars["Temp"] = &test;
                
                export_spoints_XML("cas3", step, grid, grid, vtl::Zip::ZIPPED);
            }            


            for(unsigned int i=0;i<Number_total;i++){
                Temperature_temp[i]=id.rhs[i];
            }

            mkl_call(Number_total,B,Indices_line_B,Indices_row_B,numberofnon_nullvalue_B,id.rhs,Temperature_temp);
            
            daxpy_call(Number_total,Q ,Temperature_temp);

            for(unsigned int i=0;i<Number_total;i++){
                    id.rhs[i]=Temperature_temp[i];
            }
            
        }// end if get_my_rank() == 0 
        // Verify  the material increases of 1 degree
        number_nodes_inside=0;
        number_increasesof1degree=0;
        if(type_simulation_value==1){
            if(count==rate_save_thermo){
                    for(unsigned int parcours=0 ; parcours < Number_total ; parcours++){
                        if(material_at_nodes[parcours]>0){
                            number_nodes_inside++;
                            // It heats
                            if(test[parcours]>Temperature_verif[parcours]+4){
                                number_increasesof1degree++;
                            }
                            //It cools
                            if(test[parcours]<Temperature_verif[parcours]-4){
                                number_increasesof1degree++;
                            }
                        }
                    }
                    double verification_temp = (double) (number_increasesof1degree/number_nodes_inside);
                    // Seuil à changer 
                    if(verification_temp>0.95){
                    break;
                }
            }
        }
        count++;
        if(count==rate_save_thermo+1){
            count=1;   
        }     
        t=t+dt;
        step++;
    }
}





















int algo_thermo(int argc, char **argv, GridCreator_NEW & gridElectro)
{


    //Size of the domain, in nodes
    unsigned int N_x=(gridElectro.input_parser.lengthX_WholeDomain_Thermal/gridElectro.input_parser.delta_Thermal)+1;
    unsigned int N_y=(gridElectro.input_parser.lengthY_WholeDomain_Thermal/gridElectro.input_parser.delta_Thermal)+1;
    unsigned int N_z=(gridElectro.input_parser.lengthZ_WholeDomain_Thermal/gridElectro.input_parser.delta_Thermal)+1;

    //Spatial step [m]
    double Delta = gridElectro.input_parser.delta_Thermal;

    //Time step [s]
    double dt = gridElectro.input_parser.thermal_algo_time_step;

    //Parameter theta [-]
    double theta =gridElectro.input_parser.theta_parameter; 

    // Temperature air [K]
    double T_infiny=gridElectro.input_parser.temperature_convection;   

    // Convection Parameter
    double h=gridElectro.input_parser.convection_parameter;        

    // Total number of nodes    
    unsigned int Number_total = N_x*N_y*N_z;

    //  temps de fin
    double  t_final_thermal = gridElectro.input_parser.t_final_thermal;

    //Type of simulation
    char *type_simulation = gridElectro.input_parser.type_simulation_thermal;

    // Find gemoetry
    char *geometry_material = gridElectro.input_parser.geometry_material_thermo;

    // Case Wall
    unsigned int wall_thermo=gridElectro.input_parser.wall_thermo;

    // Thermal Distribution
    unsigned int thermal_distribution=gridElectro.input_parser.thermal_distribution;

    // Amplitude Thermal Distribution
    double  amplitude_thermal_distribution=gridElectro.input_parser.amplitude_thermal_distribution;

    // Heat Distribution
    unsigned int heat_distribution=gridElectro.input_parser.heat_distribution;

    // Amplitude Heat Distribution
    double amplitude_heat_distribution=gridElectro.input_parser.amplitude_heat_distribution;

    // Boundary  
    //  choose between "Neumann"  homogeneous Neumann condition, "Dirichlet" Dirichlet condition or "Convection" convection conditions
    
    char Dirichlet[]="Dirichlet";
    char Neumann[]="Neumann";
    char Convection[]="Convection";
    char cerveau[]="cerveau";
    char Face1[20];
    char Face2[20];
    char Face3[20];
    char Face4[20];
    char Face5[20];
    char Face6[20];

    // type_simulation=1 (brain) and type_simulation=0 (analytic tests)
    unsigned int type_simulation_value=0;


    //%%%%%%%%%%%%% Temperature initiale %%%%%%%%%%%%%%ù
    // Number of different materials
    unsigned int numberofdiffents_material = gridElectro.materials.unified_material_list.size();

    //Initilization of temperature
    double *temperature_initial= (double*) calloc(numberofdiffents_material,sizeof(double));


    if(temperature_initial == NULL){	 
        printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
        abort();
    }

    //Recuperate the temperature
    for(unsigned int i=0 ;i<numberofdiffents_material;i++){
        temperature_initial[i]=gridElectro.materials.unified_material_list[i].initial_temperature;
    }

    double *temperature_initial_face =(double*) calloc(6,sizeof(double));
    if(temperature_initial_face == NULL){	 
        printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
        abort();
    }

    for(unsigned int i=0; i<6;i++){
        temperature_initial_face[i]=gridElectro.input_parser.THERMAL_FACE_BC_VALUE[i];
    }


    //%%%%%%%%%%%%%%% Case of the brain %%%%%%%%%%%%%%%%%%%%%%%%%
    if(strcmp(type_simulation,cerveau)== 0){  
        // Case of the brain
        type_simulation_value=1;
        // By default

        // Face 1 i==0
        strcpy(Face1,Dirichlet);

        // Face 2  i==N_x-1
        strcpy(Face2,Dirichlet);

        // Face 3  j==0
        strcpy(Face3,Dirichlet);

        // Face 4  j==N_y-1
        strcpy(Face4,Dirichlet);

        // Face 5  k==0
        strcpy(Face5,Dirichlet);

        //Face 6   k==N_z-1
        strcpy(Face6,Dirichlet);

    // Analytic examples
    }else{

        type_simulation_value=0;

        // Choose by the user

        // Face 1 i==0
        const char *Face1_tmp=gridElectro.input_parser.THERMAL_FACE_BC_TYPE[0].c_str();
        strcpy(Face1,Face1_tmp);

        // Face 2  i==N_x-1
        const char *Face2_tmp=gridElectro.input_parser.THERMAL_FACE_BC_TYPE[1].c_str();
        strcpy(Face2,Face2_tmp);

        // Face 3  j==0
        const char *Face3_tmp=gridElectro.input_parser.THERMAL_FACE_BC_TYPE[2].c_str();
        strcpy(Face3,Face3_tmp);

        // Face 4  j==N_y-1
        const char *Face4_tmp=gridElectro.input_parser.THERMAL_FACE_BC_TYPE[3].c_str();
        strcpy(Face4,Face4_tmp);

        // Face 5  k==0
        const char *Face5_tmp=gridElectro.input_parser.THERMAL_FACE_BC_TYPE[4].c_str();
        strcpy(Face5,Face5_tmp);

        //Face 6   k==N_z-1
        const char *Face6_tmp=gridElectro.input_parser.THERMAL_FACE_BC_TYPE[5].c_str();
        strcpy(Face6,Face6_tmp);
    }

    // the value = 0 for Dirichlet, the value = 1 for Neumann and the value=2 for both conditions
    unsigned int *Stateofeachface = (unsigned int *) calloc(6,sizeof(unsigned int));
    if(Stateofeachface == NULL){	 
        printf("The table is not calloc().This error comes from Line %d \n",__LINE__);
        abort();
    }

    // Face1
    if(strcmp(Face1, Dirichlet)== 0){
        Stateofeachface[0]=0;
    }else if(strcmp(Face1, Neumann)== 0){
        Stateofeachface[0]=1;
    }else if(strcmp(Face1, Convection)== 0){
        Stateofeachface[0]=2;
    }else{
        printf("The face 1, the boundary condition are not correct Problem in Line %d",__LINE__);
    }

    // Face2
    if(strcmp(Face2, Dirichlet)== 0){
        Stateofeachface[1]=0;
    }else if(strcmp(Face2, Neumann)== 0){
        Stateofeachface[1]=1;
    }else if(strcmp(Face2, Convection)== 0){
        Stateofeachface[1]=2;
    }else{
        printf("The face 2, the boundary condition are not correct Problem in Line %d",__LINE__);
    }

    // Face 3 
    if(strcmp(Face3, Dirichlet)== 0){
        Stateofeachface[2]=0;
    }else if(strcmp(Face3,Neumann)== 0){
        Stateofeachface[2]=1;
    }else if(strcmp(Face3, Convection)== 0){
        Stateofeachface[2]=2;
    }else{
        printf("The face 3, the boundary condition are not correct Problem in Line %d",__LINE__);
    }

    // Face4
    if(strcmp(Face4, Dirichlet)== 0){
        Stateofeachface[3]=0;
    }else if(strcmp(Face4, Neumann)== 0){
        Stateofeachface[3]=1;
    }else if(strcmp(Face4, Convection)== 0){
        Stateofeachface[3]=2;
    }else{
        printf("The face 4, the boundary condition are not correct Problem in Line %d",__LINE__);
    }

    // Face5
    if(strcmp(Face5, Dirichlet)== 0){
        Stateofeachface[4]=0;
    }else if(strcmp(Face5, Neumann)== 0){
        Stateofeachface[4]=1;
    }else if(strcmp(Face5, Convection)== 0){
        Stateofeachface[4]=2;
    }else{
        printf("The face 5, the boundary condition are not correct Problem in Line %d",__LINE__);
    }

    // Face6
    if(strcmp(Face6, Dirichlet)== 0){
        Stateofeachface[5]=0;
    }else if(strcmp(Face6, Neumann)== 0){
        Stateofeachface[5]=1;
    }else if(strcmp(Face6, Convection)== 0){
        Stateofeachface[5]=2;
    }else{
        printf("The face 6, the boundary condition are not correct Problem in Line %d",__LINE__);
    }


    // for Paraview
    SPoints grid;

    // setup grid
    grid.o = Vec3d(0.0, 0.0, 0.0);     // origin    
    grid.np1 = Vec3i(0, 0, 0);    // first index
    grid.np2 = Vec3i(N_x-1, N_y-1, N_z-1); // last index
    grid.dx = Vec3d(Delta, Delta, Delta); // compute spacing

    //Sampling rate frequency
    unsigned int rate_save_thermo=gridElectro.input_parser.SAMPLING_FREQ_THERMAL;


    // Is the main structure
    DMUMPS_STRUC_C id;


    // Initialization of MUMPS
    init_MUMPS(id);

    // Creation of the system and Resolution

    printf("\n===============================================\n==== CALL TO RESOLVE ====\n===============================================\n");

    resolve(id,grid,Number_total,Stateofeachface, N_x, N_y, N_z, Delta ,dt,theta,h, T_infiny,type_simulation_value,temperature_initial,t_final_thermal,temperature_initial_face,rate_save_thermo,geometry_material,wall_thermo,thermal_distribution,amplitude_thermal_distribution,heat_distribution,amplitude_heat_distribution,
                gridElectro);


    // Terminate the process
    end_MUMPS(id);


    MPI_Finalize();

    return 0;
}
