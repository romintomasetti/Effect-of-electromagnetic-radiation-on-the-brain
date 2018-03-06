#include "AlgoElectro_NEW.hpp"

#include <ctime>

double AlgoElectro_NEW::Compute_dt(GridCreator_NEW &mesh){
    // Retrieve the spatial step in each direction:
    double dx = mesh.delta_Electromagn[0];
    double dy = mesh.delta_Electromagn[1];
    double dz = mesh.delta_Electromagn[2];
    // Time step:
    double dt = 0.0;
    // Temporary variable:
    double tmp= 0.0;
    // Iterator:
    unsigned char i=0;                                      

    // Iterate over the number of materials. For each material, compute the required time step.
    // At the end, the smallest time step is chosen.                                          
    for (i = 0 ; i < mesh.materials.numberOfMaterials ; i++ ){

            // Get material:
            string material = mesh.materials.materialName_FromMaterialID[i];
            // Get permeability:
            double mu_material = mesh.materials.getProperty(
                    mesh.input_parser.GetInitTemp_FromMaterialName[material],
                    i,4);    

            // Get permittivity:
            double epsilon_material = mesh.materials.getProperty(
                   mesh.input_parser.GetInitTemp_FromMaterialName[material],
                    i,5);     
            // Compute speed of light:
            double c = 1/(sqrt(mu_material*epsilon_material));
            // Take the smallest time step:
            if( i == 0 ){
                dt = 1 / ( c * sqrt( 1 / ( dx * dx ) + 1 / ( dy * dy ) + 1 / (dz *dz) ) );
            }
            else{
                tmp = 1/(c*sqrt(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)));
                if( tmp < dt ){
                    dt = tmp;
                }
            }
    }
    return dt;
}

void AlgoElectro_NEW::update(GridCreator_NEW &grid){

    // Start clock for monitoring CPU time:
    std::clock_t start_algo_update_CPU_TIME;
    std::clock_t end___algo_update_CPU_TIME;
    grid.profiler.addTimingInputToDictionnary("AlgoElectro_NEW_UPDATE");
    start_algo_update_CPU_TIME = std::clock();

    // Retrieve the time step:
    double dt = this->Compute_dt(grid);
    std::cout << "AlgoElectro_NEW :: dt is " << dt << std::endl;

    /*
     *  TEMPERATURE WILL NEVER CHANGE IN THIS ALGORITHM.
     *  INITIALIZE ALL THE COEFFICIENTS NEEDED FOR THE UPDATE EQUATIONS.
     */
    // In the object grid, set the properties mu, eps, magnetic cond. and electric cond. :
    grid.Initialize_Electromagnetic_Properties("AIR_AT_INIT_TEMP");



    // End clock for monitoring CPU time
    end___algo_update_CPU_TIME = std::clock();
    double elapsedTimeSec = (end___algo_update_CPU_TIME - start_algo_update_CPU_TIME)
                                 / (double)(CLOCKS_PER_SEC);
    grid.profiler.incrementTimingInput("AlgoElectro_NEW_UPDATE",elapsedTimeSec);
    std::cout << "AlgoElectro_NEW_UPDATE => Time: " << elapsedTimeSec << " s" << std::endl;
}