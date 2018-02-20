#include "InputParser.h"

InputParser::InputParser(string file_name){
		#if DEBUG > 2
		cout << "InputParser::constructor::IN\n";
		#endif
		this->filename = file_name;
		#if DEBUG > 2
		cout << "InputParser::constructor::OUt\n";
		#endif
}

void InputParser::defaultParsingFromFile(void){
	if(this->filename == string()){
		cout << "\nInputParser::defaultParsingFromFile\n\tNo file provided. ABORTING.\n\n";
		abort();
	}else{
		// Check that the file exist:
		if(this->is_file_exist(this->filename)){
			#if DEBUG > 1
			cout << "Input file exist !" << endl;
			#endif
		}else{
			cout << "Input file doesn't exist ! Aborting.\n";
			abort();
		}
		// Go to the parsing function:
		this->basicParsing(this->filename);
	}
}

void InputParser::defaultParsingFromFile(string filename){
	// Check that the file exists:
	if(this->is_file_exist(filename)){
		#if DEBUG > 1
		cout << "Input file exist !" << endl;
		#endif
	}else{
		cout << "Input file doesn't exist ! Aborting.\n";
		abort();
	}
	// Go to the parsing function:
	this->basicParsing(filename);
}

bool InputParser::is_file_exist(const string fileName){
	cout << "InputParser::is_file_exist::IN\n";
    ifstream infile(fileName);
	cout << "InputParser::is_file_exist::OUT\n";
    return infile.good();
}

void InputParser::basicParsing(const string filename){
	/*
	 * Here is what an input file looks like:
	 *	$INFOS <- this is of kind header 1
	 *		$NAME <- this is of kind header 2
	 */
	// Check the extension of the file:
	if(filename.substr(filename.find_last_of(".")+1) == "input"){
		// The extension is correct, proceed.
		// Open the file for reading:
		ifstream inputFile;
		inputFile.open(filename);
		if(inputFile.fail()){
			// Opening failed, aborting.
			cout << "InputParser::basicParsing::Failed to open " + filename + ". Aborting.\n";
			inputFile.clear();
			abort();
		}else if(inputFile.is_open()){
			// Contains the current read line of the input file:
			string currentLine;
			// Tell if we are in header H1
			bool isIn_H1 = false;
			
			while(!inputFile.eof()){
				getline(inputFile,currentLine);
				// Check that the line is not a comment:
				if(this->checkLineISNotComment(inputFile,currentLine)){
					// The line was a comment, ignoring it.
				}
				cout << "Current line is " + currentLine << endl;
				if(currentLine.find("$") != std::string::npos){
					cout << "There is a dollor in ::" + currentLine + "::\n";
					if(!isIn_H1){
						isIn_H1 = true;
						this->readHeader(inputFile);
					}else{
						cout << "Should not end up here !\n";
						printf("Aborting (File %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
			}
		}else{
			cout << "InputParser::basicParsing::Should not end up here ! Complain to the developer.\n";
			cout << "Aborting.\n";
			abort();
		}
		inputFile.close();
	}else{
		cout << "The input file is not under .input format. Please check your input file";
		cout << " " + filename << endl;
	}
}

void InputParser::readHeader(ifstream &file){
	cout << "Not yet Implemented. Exit. Complain to Romin.\n";
	printf("Aborting (File %s at %d)\n",__FILE__,__LINE__);
	abort();
}

// Check that the line is not a comment:
bool InputParser::checkLineISNotComment(ifstream &file, string currentLine){
	string comment1_beg = "/*";
	string comment1_end = "*/";
	string comment2 = "//";
	cout << "check received " + currentLine << endl;
	if(currentLine.find(comment1_beg) != std::string::npos){
		// We have multiple lines comment ! Read all the comment before exiting.
		cout << "Multiple line comment:\n";
		while(!file.eof()){
			if(currentLine.find(comment1_end) != std::string::npos){
				cout << "\t|" + currentLine << "|\n" << endl;
				break;
			}
			cout << "\t|" + currentLine << "|\n" << endl;
			getline(file,currentLine);
		}
		return true;
	}else if(currentLine.find(comment2) != std::string::npos){
		// We have one line comment !
		cout << "This is a one line comment : " + currentLine << endl;
	}else{
		return false;
	}
	return false;	
}