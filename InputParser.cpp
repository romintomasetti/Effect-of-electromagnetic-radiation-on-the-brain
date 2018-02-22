#include "InputParser.h"

#include <algorithm>
#include <cctype>


// Small enums for the dollar strings (see InputParser::readHeader)
stringDollar_Header1 InputParser::hashit_Header1 (std::string const& inString) {
    if (inString == "INFOS") return INFOS;
    if (inString == "MESH") return MESH;
    if (inString == "RUN_INFOS") return RUN_INFOS;
	else {
		printf("In file %s at %d. Complain to Romin. Abort().\n",__FILE__,__LINE__);
		abort();
	}
}
stringDollar_Header2 InputParser::hashit_Header2 (std::string const& inString) {
    if (inString == "NAME") return NAME;
    if (inString == "DELTAS") return DELTAS;
    if (inString == "DOMAIN_SIZE") return DOMAIN_SIZE;
	if (inString == "SOURCE") return SOURCE;
	if (inString == "STOP_SIMUL_AFTER") return STOP_SIMUL_AFTER;
	if (inString == "TEMP_INIT") return TEMP_INIT;
	else {
		printf("In file %s at %d. Complain to Romin. Abort().\n",__FILE__,__LINE__);
		cout << "Faulty string is ::" + inString + "::" << endl;
		abort();
	}
}

// Get lengths
double InputParser::get_length(unsigned int direction){
	// Direction = 0 gives along X, 1 along Y, 2 along Z:
	if(direction == 0){
		return this->lengthX;
	}else if(direction == 1){
		return this->lengthY;
	}else if(direction == 2){
		return this->lengthZ;
	}else{
		printf("InputParser::get_length::ERROR:\n\tDirection should be between 0 and 2.");
		printf("Aborting (file %s at %d).\n",__FILE__,__LINE__);
		abort();
	}
}

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
			
			while(!inputFile.eof()){
				getline(inputFile,currentLine);
				// Check that the line is not a comment:
				if(this->checkLineISNotComment(inputFile,currentLine)){
					// The line was a comment, ignoring it.
				}
				cout << "Current line is " + currentLine << endl;
				if(currentLine.find("$") != std::string::npos){
					cout << "There is a dollor in ::" + currentLine + "::\n";
					this->readHeader(inputFile,currentLine);
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

void InputParser::readHeader(ifstream &file,std::string &currentLine){
	// We are in a Dollar zone.
	// First, detect which dollar zone it is.
	// Get the string after the dollar:
	std::string strHeader1 = currentLine.substr(currentLine.find("$")+1);
	// Remove any space before the switch:
	strHeader1.erase(std::remove_if(strHeader1.begin(),
			 strHeader1.end(), [](unsigned char x){return std::isspace(x);}),
			 strHeader1.end());
	cout << "String analyzed is " + strHeader1 << endl;
	// Go with the switch:
	switch(hashit_Header1(strHeader1)){
		case INFOS    : 
			this->readHeader_INFOS(file);
			cout << "EXITING SECTION INFOS WITH currentLine=" + currentLine << endl;
			break;
		case MESH     : 
			this->readHeader_MESH (file);
			break;
		case RUN_INFOS: 
			this->readHeader_RUN_INFOS(file);
			break;
		default:
			printf("Should not end up here. Complain to Romin. Abort.");
			printf("(in file %s at %d)\n",__FILE__,__LINE__);
			abort();
	}
}

// Check that the line is not a comment:
bool InputParser::checkLineISNotComment(ifstream &file, string &currentLine){
	/* Check for line(s) being comments or blank */
	std::string str = currentLine;
	this->RemoveAnyBlankSpaceInStr(str);
	if(str == string()){
		cout << "It is a blank line !" << endl;
		getline(file,currentLine);
	}
	string comment1_beg = "/*";
	string comment1_end = "*/";
	string comment2 = "//";
	#if DEBUG > 3
	cout << "check received " + currentLine << endl;
	#endif
	if(currentLine.find(comment1_beg) != std::string::npos){
		// We have multiple lines comment ! Read all the comment before exiting.
		#if DEBUG > 3
		cout << "Multiple line comment:\n";
		#endif
		while(!file.eof()){
			if(currentLine.find(comment1_end) != std::string::npos){
				#if DEBUG > 2
				cout << "\t|" + currentLine << "|\n" << endl;
				#endif
				break;
			}
			#if DEBUG > 2
			cout << "\t|" + currentLine << "|\n" << endl;
			#endif
			getline(file,currentLine);
		}
		return true;
	}else if(currentLine.find(comment2) != std::string::npos){
		// We have one line comment !
		#if DEBUG > 3
		cout << "This is a one line comment : " + currentLine << endl;
		#endif
		while(!file.eof()){
			#if DEBUG > 3
			cout << "Comment 2 entering while\n" << endl;
			#endif
			if(currentLine.find(comment2) == std::string::npos){
				#if DEBUG > 2
				cout << "\t|" + currentLine << "|\n" << endl;
				#endif
				break;
			}
			#if DEBUG > 2
			cout << "\t|" + currentLine << "|\n" << endl;
			#endif
			getline(file,currentLine);
		}
	}else{
		return false;
	}
	return false;	
}

void InputParser::RemoveAnyBlankSpaceInStr(std::string &str){
	str.erase(std::remove_if(str.begin(),
			 str.end(), [](unsigned char x){return std::isspace(x);}),
			 str.end());
}

void InputParser::readHeader_INFOS(ifstream &file){
	//Inside the header(1) called 'INFOS', we have fields:
	//		1) NAME - Contains fields:
	//				a) output : name of the output files
	//				b) error  : name of the error file
	//				c) profile: name of the profiling file (cpu time, etc)
	std::string currentLine = string();
	
	while(currentLine != "INFOS"){
		// Read line:
		getline(file,currentLine);
		// Get rid of comments:
		this->checkLineISNotComment(file,currentLine);
		// Remove any blank space:
		this->RemoveAnyBlankSpaceInStr(currentLine);
		// Remove Dollar sign:
		currentLine = currentLine.substr(currentLine.find("$")+1);
		cout << "BEFORE THE SWITCH : " + currentLine << endl;
		if(currentLine == "INFOS"){break;}
		switch(this->hashit_Header2(currentLine)){
			case NAME:
				cout << "Entering case name\n";
				while(!file.eof()){
					cout << "Entering while\n";
					// Note: sections are ended by $the-section-name.
					getline(file,currentLine);
					this->checkLineISNotComment(file,currentLine);
					this->RemoveAnyBlankSpaceInStr(currentLine);
					cout << currentLine << endl;
					if(currentLine == "$NAME"){
						cout << "EXITING NAME\n";
						break;
					}
					if(currentLine == string()){continue;}
					std::size_t posEqual  = currentLine.find("=");
					std::string propName  = currentLine.substr(0,posEqual); 
					std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());
					cout << propName + "=" + propGiven << endl;
					cout << "To compare with " + currentLine << endl;
					if(propName != "output" && propName != "error" && propName != "profile"){
						printf("InputParser::readHeader_INFOS:: You didn't provide a ");
						printf("good member for $INFOS$NAME.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
					this->outputNames[propName] = propGiven;
				}
				break;
			default:
				printf("Should not end up here. Complain to Romin. Abort.");
				printf("(in file %s at %d)\n",__FILE__,__LINE__);
				cout << "Faulty line is " + currentLine << endl << endl;
				abort();
		}
	}
}

void InputParser::readHeader_MESH (ifstream &file){
	/* Inside the section MESH, thee are fields:
	 *		1) DELTAS containing:
	 * 				a) deltaX
	 * 				b) deltaY
	 *				c) deltaZ
	 */
	std::string currentLine = string();
	while(currentLine != "MESH"){
		// Read line:
		getline(file,currentLine);
		// Get rid of comments:
		this->checkLineISNotComment(file,currentLine);
		// Remove any blank space:
		this->RemoveAnyBlankSpaceInStr(currentLine);
		// Remove Dollar sign:
		currentLine = currentLine.substr(currentLine.find("$")+1);
		cout << "BEFORE THE SWITCH : " + currentLine << endl;
		if(currentLine == "MESH"){break;}
		switch(this->hashit_Header2(currentLine)){
			case DELTAS:
				cout << "Entering case DELTAS\n";
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					cout << currentLine << endl;
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$DELTAS"){
						cout << "EXITING DELTAS\n";
						break;
					}
					// If the string is empty, it was just a white space. Continue.
					if(currentLine == string()){continue;}
					// Find the position of the equal sign:
					std::size_t posEqual  = currentLine.find("=");
					// The property we want to set:
					std::string propName  = currentLine.substr(0,posEqual);
					// The property name the user gave:
					std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());

					cout << propName + "=" + propGiven << endl;
					cout << "To compare with " + currentLine << endl;

					if(propName == "deltaX"){
						this->deltaX = std::stod(propGiven);
						cout << "ADDED DELTAX is " << this->deltaX << endl;
					}else if(propName == "deltaY"){
						this->deltaY = std::stod(propGiven);
						cout << "ADDED DELTAY is " << this->deltaX << endl;
					}else if(propName == "deltaZ"){
						this->deltaZ = std::stod(propGiven);
						cout << "ADDED DELTAZ is " << this->deltaX << endl;
					}else if(propName != "deltaX" 
							&& propName != "deltaY" 
							&& propName != "deltaZ"){
						printf("InputParser::readHeader_MESH:: You didn't provide a ");
						printf("good member for $MESH$DELTAS.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;
			
			case DOMAIN_SIZE:
				cout << "Entering case DOMAIN_SIZE\n";
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					cout << currentLine << endl;
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$DOMAIN_SIZE"){
						cout << "EXITING DOMAIN_SIZE\n";
						break;
					}
					// If the string is empty, it was just a white space. Continue.
					if(currentLine == string()){continue;}
					// Find the position of the equal sign:
					std::size_t posEqual  = currentLine.find("=");
					// The property we want to set:
					std::string propName  = currentLine.substr(0,posEqual);
					// The property name the user gave:
					std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());

					cout << propName + "=" + propGiven << endl;
					cout << "To compare with " + currentLine << endl;

					if(propName == "L_X"){
						this->lengthX = std::stod(propGiven);
						cout << "ADDED lengthX is " << this->lengthX << endl;
					}else if(propName == "L_Y"){
						this->lengthY = std::stod(propGiven);
						cout << "ADDED lengthY is " << this->lengthX << endl;
					}else if(propName == "L_Z"){
						this->lengthZ = std::stod(propGiven);
						cout << "ADDED lengthZ is " << this->lengthX << endl;
					}else if(propName != "L_X" 
							&& propName != "L_Y" 
							&& propName != "L_Z"){
						printf("InputParser::readHeader_MESH:: You didn't provide a ");
						printf("good member for $MESH$DOMAIN_SIZE.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;

			case SOURCE:
				{
					bool nbr_Sources_Defined = false;
					std::cout << "Entering case SOURCE\n";
					while(!file.eof()){
						// Note: sections are ended by $the-section-name.
						// Read line:
						getline(file,currentLine);
						// Get rid of comments:
						this->checkLineISNotComment(file,currentLine);
						// Remove any blank in the string:
						this->RemoveAnyBlankSpaceInStr(currentLine);
						// If the string is "$DELTAS" it means the section ends.
						if(currentLine == "$SOURCE"){
							std::cout << "EXITING SOURCE\n";
							break;
						}
						// If the string is empty, it was just a white space. Continue.
						if(currentLine == string()){continue;}
						// Find the position of the equal sign:
						std::size_t posEqual  = currentLine.find("=");
						// The property we want to set:
						std::string propName  = currentLine.substr(0,posEqual);
						// The property name the user gave:
						std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());

						std::cout << propName + "=" + propGiven << endl;
						std::cout << "To compare with " + currentLine << endl;

						if(propName != "NBR_SOURCES" && nbr_Sources_Defined == false){
							printf("InputParser::readHeader_MESH::CASE SOURCE\n");
							printf("You must first set the number of sources. Aborting.\n");
							std::abort();
						}

						if(propName == "NBR_SOURCES"){
							this->source.set_number_of_sources(std::stod(propGiven));
							std::cout << "NBR_SOURCES set to ";
							std::cout << this->source.get_number_of_sources() << endl;
							nbr_Sources_Defined = true;

						}else if(propName == "L_X"){
							std::vector<double> temp = this->determineVectorFromStr(propGiven);
							cout << "Values are : " << temp[0] << temp[1] << endl;
							this->source.setLengthAlongOneDir(0,temp);

						}else if(propName == "L_Y"){
							std::vector<double> temp = this->determineVectorFromStr(propGiven);
							cout << temp[0] << temp[1] << endl;
							this->source.setLengthAlongOneDir(1,temp);

						}else if(propName == "L_Z"){
							std::vector<double> temp = this->determineVectorFromStr(propGiven);
							cout << temp[0] << temp[1] << endl;
							this->source.setLengthAlongOneDir(2,temp);

						}else if(propName == "C_X"){
							std::vector<double> temp = this->determineVectorFromStr(propGiven);
							cout << temp[0] << temp[1] << endl;

						}else if(propName == "C_Y"){
							std::vector<double> temp = this->determineVectorFromStr(propGiven);
							cout << temp[0] << temp[1] << endl;

						}else if(propName == "C_Z"){
							std::vector<double> temp = this->determineVectorFromStr(propGiven);
							cout << temp[0] << temp[1] << endl;

						}else if(propName == "FRQCY"){
							std::vector<double> temp = this->determineVectorFromStr(propGiven);
							cout << temp[0] << temp[1] << endl;
						
						}else if(propName != "NBR_SOURCES" 
								&& propName != "L_X" 
								&& propName != "L_Y"
								&& propName != "L_Z"){
							printf("InputParser::readHeader_MESH:: You didn't provide a ");
							printf("good member for $MESH$SOURCE.\nAborting.\n");
							cout << propName << endl;
							printf("(in file %s at %d)\n",__FILE__,__LINE__);
							abort();
						}
					}
				}
				break;

			default:
				printf("Should not end up here. Complain to Romin. Abort.");
				printf("(in file %s at %d)\n",__FILE__,__LINE__);
				std::cout << "Faulty line is " + currentLine << endl << endl;
				abort();
		}
	}
}

void InputParser::readHeader_RUN_INFOS(ifstream &file){
	std::string currentLine = string();
	while(currentLine != "RUN_INFOS"){
		// Read line:
		getline(file,currentLine);
		// Get rid of comments:
		this->checkLineISNotComment(file,currentLine);
		// Remove any blank space:
		this->RemoveAnyBlankSpaceInStr(currentLine);
		// Remove Dollar sign:
		currentLine = currentLine.substr(currentLine.find("$")+1);
		cout << "BEFORE THE SWITCH : " + currentLine << endl;
		if(currentLine == "RUN_INFOS"){break;}
		switch(this->hashit_Header2(currentLine)){
			case STOP_SIMUL_AFTER:
				cout << "Entering case STOP_SIMUL_AFTER\n";
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					cout << currentLine << endl;
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$STOP_SIMUL_AFTER"){
						cout << "EXITING STOP_SIMUL_AFTER\n";
						break;
					}
					// If the string is empty, it was just a white space. Continue.
					if(currentLine == string()){continue;}
					// Find the position of the equal sign:
					std::size_t posEqual  = currentLine.find("=");
					// The property we want to set:
					std::string propName  = currentLine.substr(0,posEqual);
					// The property name the user gave:
					std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());

					cout << propName + "=" + propGiven << endl;
					cout << "To compare with " + currentLine << endl;

					if(propName == "stopTime"){
						this->stopTime = std::stod(propGiven);
						cout << "ADDED stopTime is " << this->stopTime << endl;
					}else if(propName != "stopTime"){
						printf("InputParser::readHeader_RUN_INFOS:: You didn't provide a ");
						printf("good member for $RUN_INFOS$STOP_SIMUL_AFTER.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;

			case TEMP_INIT:
				cout << "Entering case TEMP_INIT\n";
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					cout << currentLine << endl;
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$TEMP_INIT"){
						cout << "EXITING TEMP_INIT\n";
						break;
					}
					// If the string is empty, it was just a white space. Continue.
					if(currentLine == string()){continue;}
					// Find the position of the equal sign:
					std::size_t posEqual  = currentLine.find("=");
					// The property we want to set:
					std::string propName  = currentLine.substr(0,posEqual);
					// The property name the user gave:
					std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());

					cout << propName + "=" + propGiven << endl;
					cout << "To compare with " + currentLine << endl;

					// Create a substring with the material:
					std::string mat = propName.substr(
							propName.find("T_INIT")+sizeof("T_INIT"));
					cout << "SUBSTR IS " + mat << endl;

					if(mat != string()){
						double tempInit = std::stod(propGiven);
						this->GetInitTemp_FromMaterialName.insert(
							std::pair<std::string,double>(mat,tempInit)
						);
						cout << "For material " + mat;
						cout << ", we have initial temperature of ";
						cout << this->GetInitTemp_FromMaterialName[mat] << endl;
					}else{
						printf("InputParser::readHeader_RUN_INFOS:: You didn't provide a ");
						printf("good member for $RUN_INFOS$TEMP_INIT.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;

			default:
				printf("Should not end up here. Complain to Romin. Abort.");
				printf("(in file %s at %d)\n",__FILE__,__LINE__);
				std::cout << "Faulty line is " + currentLine << endl << endl;
				abort();
		}
	}
}

std::vector<double> InputParser::determineVectorFromStr(std::string str){
	std::stringstream stream(str);
	std::string word;
	std::vector<double> tempVec;
	cout << "STRING IS " + str << endl;
	while( getline(stream, word, ';') ){
		std::cout << word << "\n";
		tempVec.push_back(std::stod(word));
	}
		
	cout << "END\n";

	return tempVec;
}