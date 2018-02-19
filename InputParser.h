/* This class implements some parsing functions for input files */
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class InputParser{
	private:
		// File name of the input file. Should be a .input file.
		string filename = string();
		// Check a file exists:
		bool is_file_exist(const string filename);
		// Parsing function:
		void basicParsing(void);
	public:
		// Default constructor:
		InputParser(void){};
		// Constructor:
		InputParser(string filename) : filename(filename) {};
		// Destructor:
		~InputParser(void){};
		// Default parser, using the field 'filename' of the class:
		void defaultParsingFromFile(void);
		// Parser:
		void defaultParsingFromFile(string filename);
};

#endif