/* This class implements some parsing functions for input files */
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>

#include "SetOnceVariable_Template.h"
#include "ElectromagneticSource.h"

using namespace std;

enum stringDollar_Header1{
    INFOS,
	MESH,
	RUN_INFOS
};
enum stringDollar_Header2{
	NAME,
	DELTAS,
	DOMAIN_SIZE,
	SOURCE,
	STOP_SIMUL_AFTER,
	TEMP_INIT,
	MATERIALS
};

class InputParser{
	private:
		// File name of the input file. Should be a .input file.
		string filename;
		// Check a file exists:
		bool is_file_exist(const string filename);
		// Parsing function:
		void basicParsing(const string filename);
		// Check that the line is not a comment:
		bool checkLineISNotComment(ifstream &, string &);
		// Read header 1:
		void readHeader(ifstream &,std::string &);

		std::vector<double> determineVectorFromStr(std::string);

		void readHeader_INFOS(ifstream &file);
		void readHeader_MESH (ifstream &file);
		void readHeader_RUN_INFOS(ifstream &file);

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
		void defaultParsingFromFile(void);
		// Parser:
		void defaultParsingFromFile(string filename);
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
		unsigned int ratio_EM_TH_delta = 0;
		
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
};

#endif