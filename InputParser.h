/* This class implements some parsing functions for input files */
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>

class InputParser{
	private:
		// File name of the input file. Should be a .input file.
		string filename = string();
	public:
		
		InputParser(void){};
		InputParser(string filename) : filename(filename) {};
		
};

#endif