/* This class implements some parsing functions for input files */
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>

#include "SetOnceVariable_Template.h"
#include "ElectromagneticSource.h"

#include "header_with_all_defines.hpp"

#include <boost/filesystem.hpp>

#include <UTILS/directory_searching.hpp>

using namespace std;

typedef struct probed_point{
	std::string type_field;
	std::vector<double> coordinates;
	std::string at_which_time;
	std::string filename;
}probed_point;

typedef struct probed_line{
	std::string              type_field;
	std::vector<double>      coords;
    std::vector<std::string> ALL;
	std::string              filename;
    
    void print(void){
        printf("\t> Field:    %s.\n",type_field.c_str());
        printf("\t> Coords:   [%lf,%lf,%lf].\n",coords[0],coords[1],coords[2]);
        printf("\t> ALL:      {%s,%s,%s}.\n",ALL[0].c_str(),ALL[1].c_str(),ALL[2].c_str());
        printf("\t> Filename: %s.\n",filename.c_str());
    }
    
}probed_line;

enum stringDollar_Header1{
    INFOS,
	MESH,
	RUN_INFOS,
	POST_PROCESSING
};
enum stringDollar_Header2{
	NAME,
	REMOVE_EXISTING_FILES,
	DELTAS,
	DOMAIN_SIZE,
	SOURCE,
	TIME_STEP,
	OUTPUT_SAVING,
	STOP_SIMUL_AFTER,
	TEMP_INIT,
	BOUNDARY_CONDITIONS_THERMO,
	BOUNDARY_CONDITIONS_ELECTRO,
	MATERIALS,
	ORIGINS,
	PROBING_POINTS,
	ELECTRO_STEADY_STATE,
    ALGORITHM_TO_APPLY
};

class InputParser{
	private:
		// File name of the input file. Should be a .input file.
		string filename;
		// Check a file exists:
		bool is_file_exist(const string &filename);
		// Parsing function:
		void basicParsing(const string &filename);
		// Check that the line is not a comment:
		bool checkLineISNotComment(ifstream &, string &);
		// Read header 1:
		void readHeader(ifstream &,std::string &);

		/**
		 * @brief Determine a vector of double from a string.
		 * 
		 * The size of the final vector should be provided.
		 * 
		 */
		std::vector<double> determineVectorFromStr(
			std::string,
			size_t size_to_verify_for = 0 );

		void readHeader_INFOS          (ifstream &file);
		void readHeader_MESH           (ifstream &file);
		void readHeader_RUN_INFOS      (ifstream &file);
		void readHeader_POST_PROCESSING(ifstream &file);

		void RemoveAnyBlankSpaceInStr(std::string &);


		// Time (in sec) after which the simulation must be stopped):
		double stopTime = -1.0;

		
		SetOnceVariable_Template<string> simulationType;
		

		/* 
		 * All the properties read in the input file:
		 */
		// Contains error, output and profiling files:
		map<std::string,std::string> outputNames;
	public:
		/// Directory containing all material data files:
		boost::filesystem::path material_data_directory;
		/// File containing geometry parameters:
		std::string file_containing_geometry = std::string();

		/// Probed points:
		std::vector<probed_point> points_to_be_probed;
		
		/// Probed lines:
		std::vector<probed_line> lines_to_be_probed;

		/// Linked to the source behaviour:
		std::vector<std::string> source_time;

		/// Name of the file containing the materials' data:
		std::vector<std::string> material_data_files;
		
		/**
		 * Either 'dipole' or 'simple' or 'FACE_EX'
		 */		
		std::vector<std::string> conditionsInsideSources;

		// Thermal algorithm time step:
		double thermal_algo_time_step = -1;

		// Sampling frequency for the electromagnetic algorithm:
		size_t SAMPLING_FREQ_ELECTRO = 0;
		// Sampling frequency for the thermal algorithm:
		size_t SAMPLING_FREQ_THERMAL = 0;

		// Dictionary for delete operations before computing anything:
		map<std::string,bool> removeWhat_dico;

		// Origin of the grids:
		std::vector<double> origin_Electro_grid = {0,0,0};
		std::vector<double> origin_Thermal_grid = {0,0,0};

		// Map that contains, for each material, the initial temperature:
		map<std::string,double> GetInitTemp_FromMaterialName;

		// Set simulationType
		void set_SimulationType(const string str){
			this->simulationType = str;
		}

		string get_SimulationType(void){
			return this->simulationType.get();
		}

		// Source:
		ElectromagneticSource source;

		// Default constructor:
		int MPI_rank = 0;
		InputParser(int MPI_rank){this->MPI_rank = MPI_rank;};
		// Constructor:
		InputParser(string file_name);
		// Destructor:
		~InputParser(void){
			printf("[MPI %d] - InputParser destructor.\n",this->MPI_rank);
		};
		
		// Default parser, using the field 'filename' of the class:
		void defaultParsingFromFile(int MPI_RANK = 0);
		// Parser:
		void defaultParsingFromFile(std::string &filename,int MPI_RANK = 0);
		// Get lengths
		double get_length_WholeDomain(
			unsigned int /*DIRECTION 0, 1 or 2*/,
			std::string /*type: "Electro" or "Thermal"*/);

		stringDollar_Header1 hashit_Header1 (std::string const& inString);
		stringDollar_Header2 hashit_Header2 (std::string const& inString);

		map<std::string,std::string> get_outputNames(void){
			return this->outputNames;
		}

		double get_stopTime(void){return this->stopTime;}

		// Spatial step for the electromagnetic grid:
		double deltaX_Electro = 0.0;
		double deltaY_Electro = 0.0;
		double deltaZ_Electro = 0.0;
		/* 
		 * Spatial step for the thermal grid, considered as homogeneous,
		 * i.e. deltaX_therm = deltaY_therm = deltaZ_therm
		 */
		double delta_Thermal = 0.0;

		// Ratio between thermal and electromagnetic spatial steps:
		double ratio_EM_TH_delta = 0;
		
		/*
		 * Length of the domain in each direction.
		 * The length is given in meters.
		 */
		double lengthX_WholeDomain_Electro = 0.0;
		double lengthY_WholeDomain_Electro = 0.0;
		double lengthZ_WholeDomain_Electro = 0.0;
		double lengthX_WholeDomain_Thermal = 0.0;
		double lengthY_WholeDomain_Thermal = 0.0;
		double lengthZ_WholeDomain_Thermal = 0.0;

		size_t maxStepsForOneCycleOfElectro = 0;
		
		
		// Final time for thermo algo:
		double t_final_thermal = 0.0;

		// Theta parameter thermal:
		double theta_parameter = 0.0;

		// Type of simulation thermal:
		char *type_simulation_thermal = NULL;
        
        // Geometry file for thermal algorithm:
        char *geometry_material_thermo = NULL;

		// Convection parameter of air:
		double convection_parameter = 0.0;

		// Temperature infiny:
		double temperature_convection=0.0;

		// Case of the wall (analytic):
		unsigned int wall_thermo = 0;

		// Thermal distribution:
		unsigned int thermal_distribution=0;

		std::map<std::string,std::string> TEST_PARAVIEW_MPI_ARGS;

		void deleteFiles(int MPI_RANK = 0);

		/// Boundary conditions for the thermal algorithm.
		/// Example: access BC type of face 0 by THERMAL_FACE_BC_TYPE[0].
		map<size_t,std::string> THERMAL_FACE_BC_TYPE;
		map<size_t,double>      THERMAL_FACE_BC_VALUE;
        
        /// For the steady-state detection:
        size_t SteadyState_CheckEveryPoint = 0;
        bool   check_steady_state          = false;

		/////////////////////////////////////////////////
		/// PARAMETER FOR THE ABC BOUNDARY CONDITIONS ///
		/////////////////////////////////////////////////

		bool   apply_ABC_BCs                    = false;

		/////////////////////////////////////////////////
		/// PARAMETERS FOR THE PML BOUNDARY CONDITONS ///
		/////////////////////////////////////////////////

		bool   apply_PML_BCs                    = false;
		size_t thickness_PML_in_number_of_nodes = std::numeric_limits<std::size_t>::max();
		double PML_order                        = std::numeric_limits<double>::max();
		double PML_sigma_M                      = std::numeric_limits<double>::max();

		/////////////////////////////////////////////////
		/// END OF PARAMETERS FOR THE PML BOUNDARY CONDITONS ///
		/////////////////////////////////////////////////
        
        ////// Which algorithm is applied //////
        bool apply_thermo_algo  = false;
        bool apply_electro_algo = false;
};

#endif
