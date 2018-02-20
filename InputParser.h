/* This class implements some parsing functions for input files */
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include <iostream>
#include <fstream>

#include "SetOnceVariable_Template.h"
#include "ElectromagneticSource.h"

using namespace std;

class InputParser{
	private:
		// File name of the input file. Should be a .input file.
		string filename;
		// Check a file exists:
		bool is_file_exist(const string filename);
		// Parsing function:
		void basicParsing(const string filename);
		// Check that the line is not a comment:
		bool checkLineISNotComment(ifstream &file, string currentLine);
		// Read header 1:
		void readHeader(ifstream &file);
		// INFOS - NAME:
		SetOnceVariable_Template<string> nameOfSimulation;
		SetOnceVariable_Template<string> nameOfErrorLogFile;
		SetOnceVariable_Template<string> nameOfProfileFile;
		// Deltas:
		double deltaX = 0.0, deltaY = 0.0, deltaZ = 0.0;
		// Domain size:
		double lengthX = 0.0, lengthY = 0.0, lengthZ = 0.0;
		// Source:
		ElectromagneticSource source;
	public:
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
};

#endif