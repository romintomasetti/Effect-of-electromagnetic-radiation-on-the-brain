#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <mkl.h>
#include <stdio.h>
#include "vtl.h"
#include "vtlSPoints.h"

#include "../CREATE_GEOMETRY/readInputGeometryFile.h"

#include <omp.h>
#include "mpi.h"
#include "dmumps_c.h"
using namespace vtl;


#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif

//Function get_my_rank
int get_my_rank()
{
    int myid;
    int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
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


//Function to Put the Dirichlet conditions
void Dirichlet_Boundary(double *B, int *Indices_line_B, int *Indices_row_B, double *A, int *Indices_line_A, int *Indices_row_A,
                 int *counter_nonvalue_A,int *counter_nonvalue_B, int parcours){
    B[*counter_nonvalue_B] = 1;
    Indices_line_B[*counter_nonvalue_B] = parcours+1;
    Indices_row_B[*counter_nonvalue_B] = parcours+1;
    A[*counter_nonvalue_A] =1;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1;     
    (*counter_nonvalue_A)++;
    (*counter_nonvalue_B)++;
}

//Function to Put the Convection conditions
void Convection_Boundary(double *A,int *Indices_line_A,int * Indices_row_A,int *counter_nonvalue_A,
                                int parcours,int whichface,double Delta,int N_x, int N_y,int N_z, double h, double T_infiny, double k, double *convection_contribution){
    //!!!!!!!!!!!!!!! regarder le nombre de valeur non-nulle
    //Face 1 i==0
    if(whichface==0){
        A[*counter_nonvalue_A] =-k/Delta+h;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k/Delta;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+2;     
        (*counter_nonvalue_A)++;
        convection_contribution[parcours]=h*T_infiny;
    }

    // Face 2 i==N_x-1
    if(whichface==1){
        A[*counter_nonvalue_A] =-k/Delta+h;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k/Delta;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours;     
        (*counter_nonvalue_A)++;
        convection_contribution[parcours]=h*T_infiny;
    }
    // Face 3 j==0
    if(whichface==2){
        A[*counter_nonvalue_A] =-k/Delta+h;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k/Delta;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;     
        (*counter_nonvalue_A)++;
        convection_contribution[parcours]=h*T_infiny;
    }

    // Face 4 j==N_y-1
     if(whichface==3){
        A[*counter_nonvalue_A] =-k/Delta+h;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k/Delta;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
        (*counter_nonvalue_A)++;
        convection_contribution[parcours]=h*T_infiny;
    }

    // Face 5 k==0
    if(whichface==4){
        A[*counter_nonvalue_A] =-k/Delta+h;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k/Delta;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x*N_y;     
        (*counter_nonvalue_A)++;
        convection_contribution[parcours]=h*T_infiny;
    }
    // Face 6 k==N_z-1
    if(whichface==5){
        A[*counter_nonvalue_A] =-k/Delta+h;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1;
        (*counter_nonvalue_A)++;
        A[*counter_nonvalue_A] =k/Delta;
        Indices_line_A[*counter_nonvalue_A]=parcours+1;
        Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x*N_y;     
        (*counter_nonvalue_A)++;
        convection_contribution[parcours]=h*T_infiny;
    }
}

//Function to Put the Neumann conditions
void Neumann_Boundary(double *A, int *Indices_line_A, int *Indices_row_A, int *counter_nonvalue_A, int parcours,int whichface,double Delta,int N_x,int N_y,int N_z){
    //Face 1 i==0
    if(whichface==0){
    A[*counter_nonvalue_A] =1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1;     
    (*counter_nonvalue_A)++;
    A[*counter_nonvalue_A] =-1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+2;     
    (*counter_nonvalue_A)++;
    }

    //Face 2  i==N_x-1
    if(whichface==1){
    A[*counter_nonvalue_A] =1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1;     
    (*counter_nonvalue_A)++;
    A[*counter_nonvalue_A] =-1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours;     
    (*counter_nonvalue_A)++;
    }

    //Face 3  j==0
    if(whichface==2){
    A[*counter_nonvalue_A] =1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1;     
    (*counter_nonvalue_A)++;
    A[*counter_nonvalue_A] =-1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x;     
    (*counter_nonvalue_A)++;
    }

    //Face 4  j==N_y-1
    if(whichface==3){
    A[*counter_nonvalue_A] =1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1;     
    (*counter_nonvalue_A)++;
    A[*counter_nonvalue_A] =-1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x;     
    (*counter_nonvalue_A)++;
    }

    //Face 5  k==0
    if(whichface==4){
    A[*counter_nonvalue_A] =1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1;     
    (*counter_nonvalue_A)++;
    A[*counter_nonvalue_A] =-1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1+N_x*N_y;     
    (*counter_nonvalue_A)++;
    }

    //Face 6  k==N_z-1
    if(whichface==5){
    A[*counter_nonvalue_A] =1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1;     
    (*counter_nonvalue_A)++;
    A[*counter_nonvalue_A] =-1/Delta;
    Indices_line_A[*counter_nonvalue_A]=parcours+1;
    Indices_row_A[*counter_nonvalue_A]=parcours+1-N_x*N_y;     
    (*counter_nonvalue_A)++;
    }
}

//Function to give the number of non zero in matrice A and B
 void Numberofnon_zero_function(MKL_INT *numberofnon_nullvalue_A,MKL_INT *numberofnon_nullvalue_B, double *Stateofeachface, int N_x,int N_y,int N_z){
   
    

   // Face 1 i==0
    if(Stateofeachface[0]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+N_y*N_z;
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+N_y*N_z;
    }else if(Stateofeachface[0]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+N_y*N_z*2;
     }else if(Stateofeachface[0]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+N_y*N_z*2;
       
    }
    // Face 2 i==N_x-1
    if(Stateofeachface[1]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+N_y*N_z;
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+N_y*N_z;
    }else if(Stateofeachface[1]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+N_y*N_z*2;
     }else if(Stateofeachface[1]==2){
         *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+N_y*N_z*2;
    }

    // Face 3 j==0
    if(Stateofeachface[2]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*N_z;
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x-2)*N_z;
    }else if(Stateofeachface[2]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*N_z*2;
     }else if(Stateofeachface[2]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*N_z*2;
    }

    // Face 4 j==N_y-1
    if(Stateofeachface[3]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*N_z;
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x-2)*N_z;
    }else if(Stateofeachface[3]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*N_z*2;
     }else if(Stateofeachface[3]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*N_z*2;
    }

    // Face 5 k==0
    if(Stateofeachface[4]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2);
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x-2)*(N_y-2);
    }else if(Stateofeachface[4]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2)*2;
     }else if(Stateofeachface[4]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2)*2;
    }

    // Face 6 k==N_z-1
    if(Stateofeachface[5]==0){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2);
        *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+(N_x-2)*(N_y-2);
    }else if(Stateofeachface[5]==1){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2)*2;
     }else if(Stateofeachface[5]==2){
        *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+(N_x-2)*(N_y-2)*2;
    }
    
    
    // Middle Contribution
    int contribution_middle=7*(N_x-2)*(N_y-2)*(N_z-2);
    
    *numberofnon_nullvalue_A=*numberofnon_nullvalue_A+contribution_middle;
    *numberofnon_nullvalue_B=*numberofnon_nullvalue_B+contribution_middle;
 }








//Function to give the Temperature initial
void InitializeTemperature(double *InitialTemperature, int Number_eqtosolve, int N_x, int N_y, int N_z, double T_infiny, double Delta){

    
    // A DEFINIR


        //Boundary conditions
        //Face 1 i==0
        int Number_total=N_x*N_y*N_z;
        /*for(int parcours=0;parcours<Number_total;parcours++){
            int position=parcours;
            int tmp= position/(N_x*N_y);
            position=position-tmp*(N_x*N_y);
            int tmp_2=position/(N_x);
            position=position-tmp_2*N_x;
            InitialTemperature[parcours]=sin(M_PI*(position*Delta)/((N_x-1)*Delta));
        }*/
        for(int parcours=0;parcours<Number_total;parcours++){
            int tmp = floor(parcours/(N_x*N_y));
            int value_testx= parcours-tmp*N_x*N_y;
             if (value_testx < N_x){
                InitialTemperature[parcours]=313.15;
             }else{
		        InitialTemperature[parcours] = 298.15;
             }
	    }
}
//Function to read the file containing the heat source  in each node
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
void mkl_call(int Number_total,double *B, int *Indices_line_B, int *Colonn_line_B,int numberofnon_nullvalue_B, double *Temperature_before,double *BY_results){
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
void resolve(DMUMPS_STRUC_C &id, vtl::SPoints &grid, int Number_eqtosolve, double *Stateofeachface, int N_x, int N_y, int N_z, double Delta, double dt,double theta,double h,double T_infiny)
{
 
    
    
    // setup grid
    
    //std::vector<MUMPS_INT> irn;
    //std::vector<MUMPS_INT> jcn;
    //std::vector<double> A;
    

    // Sans GridCreator	
    int i;
	
	int Number_total = N_x * N_y * N_z;

	
    // Avec GridCreator
	// Tant qu'on a pas GridCreator, mettre ce bloc en commentaire (et l'argument de la fonction aussi)
    // int N_x=mesh.Numberofnodes_x;
	// int N_y=mesh.Numberofnodes_y;
	// int N_z=mesh.NUmberofnodes_z;
	// unsigned int parcours=0;
	// double dt=mesh.deltat_Thermal;
	// double Delta=mesh.delta_Thermal;
	// double* alpha_temp=mesh.thermal_diffusivity;
	// int Number_total = N_x * N_y * N_z;

    

    //Number of non zero values in matrix A
  	MKL_INT  numberofnon_nullvalue_A=0;
    MKL_INT  numberofnon_nullvalue_B=0;

    Numberofnon_zero_function(&numberofnon_nullvalue_A,&numberofnon_nullvalue_B, Stateofeachface, N_x, N_y, N_z);
    
    // Variable for the matrix A and B
    

    //counter y
    double count_y = 0;

    // Decleration for Matrix B
    double *B = (double *) calloc(numberofnon_nullvalue_B,sizeof(double));    
    int *Indices_line_B = (int *) calloc(numberofnon_nullvalue_B,sizeof(int));
    int *Indices_row_B = (int *) calloc(numberofnon_nullvalue_B,sizeof(int));
    double *A = (double *) calloc(numberofnon_nullvalue_A,sizeof(double));    
    int *Indices_line_A = (int *) calloc(numberofnon_nullvalue_A,sizeof(int));
    int *Indices_row_A = (int *) calloc(numberofnon_nullvalue_A,sizeof(int));
    int whichface=40;
    double *convection_contribution = (double *) calloc(Number_eqtosolve,sizeof(double));
    int step=0;
    /// Read the file and put it inside a vector:
    size_t size_mat = 0;
	unsigned int* material_at_nodes = read_input_geometry_file("test.geometry",&size_mat);
    if(size_mat != N_x * N_y * N_z){
        printf("%d",size_mat);
        printf("taille erreur!! %d",__LINE__);
         abort();
    }
       



    // Verification of the calloc
     if(B == NULL || Indices_line_B == NULL || Indices_row_B == NULL || A == NULL || Indices_line_A == NULL || Indices_row_A == NULL || convection_contribution == NULL)
    {
        printf("FDP au moins un tableau n'est pas alloué. Aborting...");
        exit(EXIT_FAILURE);
    }

    // count the number of non-value
    int counter_nonvalue_A=0;
    int counter_nonvalue_B=0;

    //variable to run through all nodes
    int parcours;
    double k_conv=46;
    // Test wall 2 : material 1) concrete   ; 2) wood cross grain yellow pine  3) glass fiber,coated;duct liner
    double rho[3]={2300,640,32};
    double k[3]={1.4,0.15,0.038};
    double c_p[3]={880,2805,835};
    // Material aluminium
    /*double k[2]{46,46};
    double rho[2]={3970,3970};
    double c_p[2]={765,765};*/

    /*double k[2]={1,1};
    
    double c_p[2]={1,1};
    double rho[2]={1,1};*/
   // double alpha[3]={20E-6,1.203E-6,0.58E-6};
    // Filling matrix A and B
    std::vector<double> test(Number_total);
    // Vector of the heat source
    double *Q = (double *) calloc(Number_eqtosolve,sizeof(double));
    if(Q == NULL){
        printf("FDP ton tableau n'est pas calloc(). Cette erreur vient de la ligne %d \n",__LINE__);
        abort();
    }
    double Cst_A;
    double Cst_B;
    if(get_my_rank() == 0)
    {
        for (parcours=0; parcours<Number_eqtosolve; parcours++){
            
            //Variables use for Face 1 and Face 2
            int tmp = floor(parcours/(N_x*N_y));
            int value_testx = parcours-tmp*N_x*N_y;

            //Boundary conditions
            
            //Face 1  i==0
            if(count_y == 0){
                if(Stateofeachface[0]==0){
                    Dirichlet_Boundary(B, Indices_line_B, Indices_row_B,
                                        A, Indices_line_A, Indices_row_A, &counter_nonvalue_A,&counter_nonvalue_B, parcours);

                }else if(Stateofeachface[0]==1){
                    whichface=0;
                    Neumann_Boundary(A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z);
                }else if(Stateofeachface[0]==2){
                    whichface=0;
                    Convection_Boundary(A, Indices_line_A, Indices_row_A,&counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z,h,T_infiny,k_conv,convection_contribution);
                }else{
                    printf("Problem for Face 3, Line %d\n",__LINE__);
                    abort();
                }
            }
            
            //Face 2   i==N_x-1
            else if (count_y == N_x-1){
                if(Stateofeachface[1]==0){
                    Dirichlet_Boundary(B, Indices_line_B, Indices_row_B,
                                        A, Indices_line_A, Indices_row_A, &counter_nonvalue_A,&counter_nonvalue_B, parcours);

                }else if(Stateofeachface[1]==1){
                    whichface=1;
                    Neumann_Boundary(A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z);
                }else if(Stateofeachface[1]==2){
                    whichface=1;
                    Convection_Boundary( A, Indices_line_A, Indices_row_A,&counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z,h,T_infiny,k_conv,convection_contribution);
                }else{
                    printf("Problem for Face 4, Line %d\n",__LINE__);
                    abort();
                }
            }
            
            //Face 3 j==0
            else if (value_testx < N_x){
                //Dirichlet Condition to make
                if(Stateofeachface[2]==0){
                    Dirichlet_Boundary(B, Indices_line_B, Indices_row_B,
                                        A, Indices_line_A, Indices_row_A, &counter_nonvalue_A,&counter_nonvalue_B, parcours);

                }else if(Stateofeachface[2]==1){
                    whichface=2;
                    Neumann_Boundary(A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z);

                }else if(Stateofeachface[2]==2){
                    whichface=2;
                    Convection_Boundary( A, Indices_line_A, Indices_row_A,&counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z,h,T_infiny,k_conv,convection_contribution);
                }else{
                    printf("Problem for Face 1, Line %d\n",__LINE__);
                    abort();
                }
            }
                
            //Face 4  j==N_y-1
            else if(value_testx >= (N_x*N_y)-N_x){
                if(Stateofeachface[3]==0){
                    Dirichlet_Boundary(B, Indices_line_B, Indices_row_B,
                                        A, Indices_line_A, Indices_row_A, &counter_nonvalue_A,&counter_nonvalue_B, parcours);

                }else if(Stateofeachface[3]==1){
                    whichface=3;
                    Neumann_Boundary(A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z);
                }else if(Stateofeachface[3]==2){
                    whichface=3;
                    Convection_Boundary( A, Indices_line_A, Indices_row_A,&counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z,h,T_infiny,k_conv,convection_contribution);
                }else{
                    printf("Problem for Face 2, Line %d\n",__LINE__);
                    abort();
                }
            }

            //Face 5  k==0
            else if(parcours <= N_x*N_y){
                if(Stateofeachface[4]==0){
                    Dirichlet_Boundary(B, Indices_line_B, Indices_row_B,
                                        A, Indices_line_A, Indices_row_A, &counter_nonvalue_A,&counter_nonvalue_B, parcours);

                }else if(Stateofeachface[4]==1){
                    whichface=4;
                    Neumann_Boundary(A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z);
                }else if(Stateofeachface[4]==2){
                    whichface=4;
                    Convection_Boundary( A, Indices_line_A, Indices_row_A,&counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z,h,T_infiny,k_conv,convection_contribution);
                }else{
                    printf("Problem for Face 5, Line %d\n",__LINE__);
                    abort();
                }
            }

            //Face 6 k==N_z-1
            else if(parcours >= Number_total -N_x*N_y){
                if(Stateofeachface[5]==0){
                    Dirichlet_Boundary(B, Indices_line_B, Indices_row_B,
                                        A, Indices_line_A, Indices_row_A, &counter_nonvalue_A,&counter_nonvalue_B, parcours);

                }else if(Stateofeachface[5]==1){
                    whichface=5;
                    Neumann_Boundary(A, Indices_line_A, Indices_row_A, &counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z);
                }else if(Stateofeachface[5]==2){
                    whichface=5;
                    Convection_Boundary( A, Indices_line_A, Indices_row_A,&counter_nonvalue_A, parcours,whichface,Delta,N_x,N_y,N_z,h,T_infiny,k_conv,convection_contribution);
                }else{
                    printf("Problem for Face 6, Line %d\n",__LINE__);
                    abort();
                }
            }
            
            //NO face, we are in the mesh
            else{
                B[counter_nonvalue_B] = (rho[material_at_nodes[parcours]]*c_p[material_at_nodes[parcours]])-(theta*dt/(2*Delta*Delta))*(6*k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]]
                            +k[material_at_nodes[parcours+1]]+k[material_at_nodes[parcours+N_x]]+k[material_at_nodes[parcours-N_x]]+k[material_at_nodes[parcours+N_x*N_y]]+k[material_at_nodes[parcours-N_x*N_y]]);
                Indices_line_B[counter_nonvalue_B] = parcours+1;
                Indices_row_B[counter_nonvalue_B] = parcours+1;
                A[counter_nonvalue_A]=(rho[material_at_nodes[parcours]]*c_p[material_at_nodes[parcours]])+((1-theta)*dt/(2*Delta*Delta))*(6*k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]]
                            +k[material_at_nodes[parcours+1]]+k[material_at_nodes[parcours+N_x]]+k[material_at_nodes[parcours-N_x]]+k[material_at_nodes[parcours+N_x*N_y]]+k[material_at_nodes[parcours-N_x*N_y]]);
                Indices_line_A[counter_nonvalue_A]=parcours+1;
                Indices_row_A[counter_nonvalue_A]=parcours+1;
                counter_nonvalue_A++;
                counter_nonvalue_B++;
                    
                    
                //selon x
                Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+1]])*dt/(2*Delta*Delta);
                Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+1]])*dt/(2*Delta*Delta);
                B[counter_nonvalue_B] = Cst_B;
                Indices_line_B[counter_nonvalue_B] = parcours+1;
                Indices_row_B[counter_nonvalue_B] = parcours+2;
                A[counter_nonvalue_A]=Cst_A;
                Indices_line_A[counter_nonvalue_A]=parcours+1;
                Indices_row_A[counter_nonvalue_A]=parcours+2;
                counter_nonvalue_A++;
                counter_nonvalue_B++;

                Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]])*dt/(2*Delta*Delta);
                Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-1]])*dt/(2*Delta*Delta);
                B[counter_nonvalue_B] = Cst_B;
                Indices_line_B[counter_nonvalue_B] = parcours+1;
                Indices_row_B[counter_nonvalue_B] = parcours;
                A[counter_nonvalue_A]=Cst_A;
                Indices_line_A[counter_nonvalue_A]=parcours+1;
                Indices_row_A[counter_nonvalue_A]=parcours;
                counter_nonvalue_A++;
                counter_nonvalue_B++;

                //selony
                Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x]])*dt/(2*Delta*Delta);
                Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x]])*dt/(2*Delta*Delta);
                B[counter_nonvalue_B]=Cst_B;
                Indices_line_B[counter_nonvalue_B] = parcours+1;
                Indices_row_B[counter_nonvalue_B] = parcours+N_x+1;
                A[counter_nonvalue_A]=Cst_A;
                Indices_line_A[counter_nonvalue_A]=parcours+1;
                Indices_row_A[counter_nonvalue_A]=parcours+N_x+1;
                counter_nonvalue_A++;
                counter_nonvalue_B++;

                Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x]])*dt/(2*Delta*Delta);
                Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x]])*dt/(2*Delta*Delta);
                B[counter_nonvalue_B] = Cst_B;
                Indices_line_B[counter_nonvalue_B] = parcours+1;
                Indices_row_B[counter_nonvalue_B] = parcours-N_x+1;
                A[counter_nonvalue_A]=Cst_A;
                Indices_line_A[counter_nonvalue_A]=parcours+1;
                Indices_row_A[counter_nonvalue_A]=parcours-N_x+1;
                counter_nonvalue_A++;
                counter_nonvalue_B++;

                //selonz
                Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x*N_y]])*dt/(2*Delta*Delta);
                Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours+N_x*N_y]])*dt/(2*Delta*Delta);
                B[counter_nonvalue_B] = Cst_B;
                Indices_line_B[counter_nonvalue_B] = parcours+1;
                Indices_row_B[counter_nonvalue_B] = parcours+N_x*N_y+1;
                A[counter_nonvalue_A]=Cst_A;
                Indices_line_A[counter_nonvalue_A]=parcours+1;
                Indices_row_A[counter_nonvalue_A]=parcours+N_x*N_y+1;
                counter_nonvalue_A++;
                counter_nonvalue_B++;

                Cst_A = -theta*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x*N_y]])*dt/(2*Delta*Delta);
                Cst_B=(1-theta)*(k[material_at_nodes[parcours]]+k[material_at_nodes[parcours-N_x*N_y]])*dt/(2*Delta*Delta);
                B[counter_nonvalue_B] = Cst_B;
                Indices_line_B[counter_nonvalue_B] = parcours+1;
                Indices_row_B[counter_nonvalue_B] = parcours-N_x*N_y+1;
                A[counter_nonvalue_A]=Cst_A;
                Indices_line_A[counter_nonvalue_A]=parcours+1;
                Indices_row_A[counter_nonvalue_A]=parcours-N_x*N_y+1;
                counter_nonvalue_A++;
                counter_nonvalue_B++;
            }
                
            count_y=count_y+1;
            if(count_y==N_x){
                    count_y = 0;
            }
            
            
        }
    
        

        // End of the filling

        // Temperature Initial
        double *InitialTemperature = (double *) calloc(Number_eqtosolve,sizeof(double));
        InitializeTemperature(InitialTemperature,Number_eqtosolve,N_x,N_y,N_z,T_infiny,Delta);
            
        

        // Lecture of the file outside
        ReadFile(Number_total,Q,dt);


        // Variable 
        double *BY_results = (double *) calloc(Number_eqtosolve,sizeof(double));

        if(BY_results == NULL)
        {
            printf("Le vecteur TESTVECTOR n'est pas alloué. Aborting...");
            exit(EXIT_FAILURE);
        }
        // contribution de convection mettre dans Q
        daxpy_call(Number_eqtosolve,convection_contribution,Q);


        // compute B*T_0
        
        mkl_call(Number_eqtosolve,B,Indices_line_B,Indices_row_B,numberofnon_nullvalue_B,InitialTemperature,BY_results);
        for(i=0;i<Number_eqtosolve;i++){
            printf("%lf\n",BY_results[i]);
        }
        //Compute B*T_0+Q
        daxpy_call(Number_eqtosolve,Q,BY_results);
        for(i=0;i<Number_eqtosolve;i++){
            printf("%lf\n",BY_results[i]);
        }
        id.rhs=BY_results;
        
        id.n=Number_eqtosolve;

        id.nnz=numberofnon_nullvalue_A;
        
        id.irn=Indices_line_A;
            
        id.jcn=Indices_row_A;    
        
        id.a=A;

        for(i=0;i<Number_total;i++){
                    test[i]=InitialTemperature[i];
            }
        /* for(i=0;i<Number_total;i++){
                printf("%lf\n",id.rhs[i]);
            }*/
    

        grid.scalars["Temp"] = &test;        
        export_spoints_XML("sphere", step, grid, grid, vtl::Zip::ZIPPED);
        }// end if get_my_rank() == 0 
        MPI_Barrier(MPI_COMM_WORLD);  
  
    
    
    
    
    
   	

    
    //Sans-Gridcreator
    double t_final = 15000.0;
    double t = 0.0;
     
    //With Grid
    /*double t_final=mesh.temps_Final;
    double dt=mesh.deltat_Thermique;
    double t=0.0*/

    //step
  



    
    

    
    //Loop over time
    
    step++;
    step =0;
    
    double *Temperature_temp = (double *) calloc(Number_eqtosolve,sizeof(double));
   
    // LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int count =0;
    while(t < t_final){       
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
            if(count==10){
                for(i=0;i<Number_total;i++){
                        test[i]=id.rhs[i];
                }

                grid.scalars["Temp"] = &test;
                
                export_spoints_XML("sphere", step, grid, grid, vtl::Zip::ZIPPED);
            }
            for(i=0;i<Number_eqtosolve;i++){
                Temperature_temp[i]=id.rhs[i];
            }
            mkl_call(Number_eqtosolve,B,Indices_line_B,Indices_row_B,numberofnon_nullvalue_B,id.rhs,Temperature_temp);
            
            daxpy_call(Number_eqtosolve,Q,Temperature_temp);

            for(i=0;i<Number_eqtosolve;i++){
                    id.rhs[i]=Temperature_temp[i];
            }
            count++;
            if(count==11)
            count=1;
        }// end if get_my_rank() == 0         
        t=t+dt;
        step++;
    }
}





















int main(int argc, char **argv)
{

  MPI_Init(&argc, &argv);
  
  int N_x = 51;
  int N_y = 51;
  int N_z = 51;
  double Delta = 0.01;
  double dt = 0.5;
  double theta =0.5001;
  double T_infiny=298.15;
  
  //double k=0.1;
  double h=10.0;
  int i;
  int Number_total = N_x*N_y*N_z;

  // Frontiere   
  //  choose between "Neumann"  Neumann condition, "Dirichlet" Dirichlet condition or "Neu_diri" Neumann and Dirichlet conditions

  char Dirichlet[]="Dirichlet";
  char Neumann[]="Neumann";
  char convection[]="Convection";
  // Face 1 i==0
  
  char Face1[]="Neumann";

  // Face 2  i==N_x-1
  char Face2[]="Neumann";

  // Face 3  j==0
  char Face3[]="Dirichlet";

  // Face 4  j==N_y-1
  char Face4[]="Dirichlet";

  // Face 5  k==0

  char Face5[]="Neumann";

  //Face 6   k==N_z-1

  char Face6[]="Neumann";

  int Number_eqtosolve=Number_total;

  // the value = 0 for Dirichlet, the value = 1 for Neumann and the value=2 for both conditions
  double *Stateofeachface = (double *) calloc(6,sizeof(double));


  // Compute the number total of equation

  // Face1
  if(strcmp(Face1, Dirichlet)== 0){
      Stateofeachface[0]=0;
  }else if(strcmp(Face1, Neumann)== 0){
      Stateofeachface[0]=1;
  }else if(strcmp(Face1, convection)== 0){
      Stateofeachface[0]=2;
  }else{
      printf("The face 1, the boundary condition are not correct Problem in Line %d",__LINE__);
  }
  
  // Face2
   if(strcmp(Face2, Dirichlet)== 0){
      Stateofeachface[1]=0;
  }else if(strcmp(Face2, Neumann)== 0){
      Stateofeachface[1]=1;
  }else if(strcmp(Face2, convection)== 0){
        Stateofeachface[1]=2;
  }else{
      printf("The face 2, the boundary condition are not correct Problem in Line %d",__LINE__);
  }

  // Face3
   if(strcmp(Face3, Dirichlet)== 0){
      Stateofeachface[2]=0;
  }else if(strcmp(Face3,Neumann)== 0){
      Stateofeachface[2]=1;
  }else if(strcmp(Face3, convection)== 0){
      Stateofeachface[2]=2;
  }else{
      printf("The face 3, the boundary condition are not correct Problem in Line %d",__LINE__);
  }

  // Face4
   if(strcmp(Face4, Dirichlet)== 0){
      Stateofeachface[3]=0;
  }else if(strcmp(Face4, Neumann)== 0){
      Stateofeachface[3]=1;
  }else if(strcmp(Face4, convection)== 0){
      Stateofeachface[3]=2;
  }else{
      printf("The face 4, the boundary condition are not correct Problem in Line %d",__LINE__);
  }

  // Face5
  if(strcmp(Face5, Dirichlet)== 0){
      Stateofeachface[4]=0;
  }else if(strcmp(Face5, Neumann)== 0){
      Stateofeachface[4]=1;
  }else if(strcmp(Face5, convection)== 0){
      Stateofeachface[4]=2;
  }else{
      printf("The face 5, the boundary condition are not correct Problem in Line %d",__LINE__);
  }

  // Face6
  if(strcmp(Face6, Dirichlet)== 0){
      Stateofeachface[5]=0;
  }else if(strcmp(Face6, Neumann)== 0){
      Stateofeachface[5]=1;
  }else if(strcmp(Face6, convection)== 0){
      Stateofeachface[5]=2;
  }else{
      printf("The face 6, the boundary condition are not correct Problem in Line %d",__LINE__);
  }
  //int *ValuesTemperature = (int *) calloc(Number_total,sizeof(int));
  


  //MKL_INT  numberofnon_nullvalue = (N_x-2)*(N_y-2)*(N_z-2)*7 + N_x*N_y*N_z - (N_x-2)*(N_y-2)*(N_z-2);
  

  SPoints grid;

  // setup grid
  grid.o = Vec3d(0.0, 0.0, 0.0);     // origin    
  grid.np1 = Vec3i(0, 0, 0);    // first index
  grid.np2 = Vec3i(N_x-1, N_y-1, N_z-1); // last index
  grid.dx = Vec3d(Delta, Delta, Delta); // compute spacing

  /*if(ValuesTemperature == NULL){	  printf("FDP ton tableau n'est pas calloc(). Cette erreur vient de la ligne %d \n",__LINE__);
		exit(EXIT_FAILURE);
  }*/


  // Is the main structure
  DMUMPS_STRUC_C id;


  // Initialization of MUMPS
  init_MUMPS(id);

  // Creation of the system and Resolution
  
  resolve(id,grid,Number_eqtosolve,Stateofeachface, N_x, N_y, N_z, Delta ,dt,theta,h, T_infiny);
  

  // Terminate the process
  end_MUMPS(id);


  MPI_Finalize();

  return 0;
}
