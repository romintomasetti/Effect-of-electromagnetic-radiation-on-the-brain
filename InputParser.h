/* This class implements some parsing functions for input files */
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>

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
	STOP_SIMUL_AFTER
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

		void readHeader_INFOS(ifstream &file);
		void readHeader_MESH (ifstream &file);
		void readHeader_RUN_INFOS(ifstream &file);

		void RemoveAnyBlankSpaceInStr(std::string &);

		// Source:
		ElectromagneticSource source;

		/* 
		 * All the properties read in the input file:
		 */
		// Contains error, output and profiling files:
		map<std::string,std::string> outputNames;
	public:
		// Default constructor:
		InputParser(){};
		// Constructor:
		InputParser(string file_name);
		// Destructor:
		~InputParser(void){};
		// Deltas:
		double deltaX = 0.0, deltaY = 0.0, deltaZ = 0.0;
		// Domain size:
		double lengthX = 0.0, lengthY = 0.0, lengthZ = 0.0;
		// Default parser, using the field 'filename' of the class:
		void defaultParsingFromFile(void);
		// Parser:
		void defaultParsingFromFile(string filename);
		// Get lengths
		double get_length(unsigned int);

		stringDollar_Header1 hashit_Header1 (std::string const& inString);
		stringDollar_Header2 hashit_Header2 (std::string const& inString);

		map<std::string,std::string> get_outputNames(void){
			return this->outputNames;
		}
};

#endif