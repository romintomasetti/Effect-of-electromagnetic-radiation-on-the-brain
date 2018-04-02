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

#include "SetOnceVariable_Template.h"
#include "ElectromagneticSource.h"

#include "header_with_all_defines.hpp"

using namespace std;

typedef struct probed_point{
	std::string type_field;
	std::vector<double> coordinates;
	std::string at_which_time;
	std::string filename;
}probed_point;

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
	BOUNDARY_CONDITIONS,
	MATERIALS,
	ORIGINS,
	PROBING_POINTS
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
		/// Probed points:
		std::vector<probed_point> points_to_be_probed;

		/// Linked to the source behaviour:
		std::string source_time = string();

		/// Name of the file containing the materials' data:
		std::string material_data_file = string();
		
		/**
		 * Either 'dipole' or 'simple'
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
		InputParser(){};
		// Constructor:
		InputParser(string file_name);
		// Destructor:
		~InputParser(void){};
		
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
		size_t maxStepsForOneCycleOfThermal = 0;

		std::map<std::string,std::string> TEST_PARAVIEW_MPI_ARGS;

		void deleteFiles(int MPI_RANK = 0);

		/// Boundary conditions for the thermal algorithm.
		/// Example: access BC type of face 0 by THERMAL_FACE_BC_TYPE[0].
		map<size_t,std::string> THERMAL_FACE_BC_TYPE;
		map<size_t,double>      THERMAL_FACE_BC_VALUE;
};

#endif
