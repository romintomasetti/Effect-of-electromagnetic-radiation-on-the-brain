#include "InputParser.h"


#include "JSON/json.hpp"
#include "rapidjson/document.h"
using namespace rapidjson;

#include "header_with_all_defines.hpp"

#include <algorithm>
#include <cctype>

#include <cstring>

#include <sys/stat.h>

#include <cmath>

#include <dirent.h>
#include <fstream>

#include <time.h>

#define DTTMFMT "%Y-%m-%d %H:%M:%S "
#define DTTMSZ 21

static char *getDtTm (char *buff) {
    time_t t = time (0);
    strftime (buff, DTTMSZ, DTTMFMT, localtime (&t));
    return buff;
}

bool directory_exists( const std::string &directory, bool createIt );

std::vector<size_t> findCharacterInsideString(std::string str,std::string charact){
	
	std::vector<size_t> positions;
	
	size_t pos = str.find(charact, 0);
	
	while(pos != string::npos)
	{
		positions.push_back(pos);
		pos = str.find(charact,pos+1);
	}
	return positions;
}

void removeFilesOfDirectory(std::string directoryOutputFiles, std::string extension){

	DIR *dir;
	struct dirent *ent;

	if(directoryOutputFiles != std::string()){
		if ((dir = opendir (directoryOutputFiles.c_str())) != NULL) {
			while ((ent = readdir (dir)) != NULL) {
				std::string fileName = ent->d_name;
				std::string fileExt  = fileName.substr(fileName.find('.')+1);
				if(fileExt == extension){
					/*printf ("%s has extension %s (will be deleted)\n", 
						ent->d_name,
						fileExt.c_str());*/
					std::remove( (directoryOutputFiles + "/" + fileName).c_str() );
				}else{
					/*printf ("%s has extension %s (will *NOT* be deleted)\n", 
						ent->d_name,
						fileExt.c_str());*/
				}
  			}
			closedir (dir);
		} else {
			DISPLAY_WARNING(
				"Could not open directory (%s) to delete %s files.",
				directoryOutputFiles.c_str(),extension.c_str()
			);
		}
	}
}

/**
 * @brief This function returns a folder name contained inside a string and the output file's name.
 */
std::vector<std::string> get_folder_from_name_parser(std::string outputName){
    
    // First element is the folder name, second is the file name.
    std::vector<std::string> returned_folder_name = {std::string(),outputName};

    if(outputName.find('/') != std::string::npos){
        /* The output name specifies an output folder name */

        // Folder name:
        returned_folder_name[0] = outputName.substr(0,outputName.find('/'));

        // File name:
        returned_folder_name[1] = outputName.substr(outputName.find('/')+1);
    }

    return returned_folder_name;
}

/// Check if we must delete the output files (vti and pvti):
void InputParser::deleteFiles(int MPI_RANK /* = 0 by default */){
	if(MPI_RANK == 0){
		std::vector<std::string> to_be_removed = {"remove_pvti","remove_vti"};
		for(size_t it = 0 ; it < to_be_removed.size() ; it++){
			if ( this->removeWhat_dico.find(to_be_removed[it]) == this->removeWhat_dico.end() ) {
				// Not found, meaning not specified, do nothing.
			} else {
				if(this->removeWhat_dico[to_be_removed[it]] == true){
					std::string extension = to_be_removed[it];
					extension = extension.substr(extension.find('_')+1);
					std::string directoryOutputFiles 
						= get_folder_from_name_parser(this->outputNames["output"])[0];
					removeFilesOfDirectory(directoryOutputFiles,extension);
				}
			}
		}
	}
}


// Small enums for the dollar strings (see InputParser::readHeader)
stringDollar_Header1 InputParser::hashit_Header1 (std::string const& inString) {
    if (inString == "INFOS") return INFOS;
    if (inString == "MESH") return MESH;
    if (inString == "RUN_INFOS") return RUN_INFOS;
	if (inString == "POST_PROCESSING") return POST_PROCESSING;
	else {
		DISPLAY_ERROR_ABORT(
			"Nothing corresponds to %s in headers 1.\n",
			inString.c_str()
		);
	}
}
stringDollar_Header2 InputParser::hashit_Header2 (std::string const& inString) {
    if (inString == "NAME")                  return NAME;
	if (inString == "REMOVE_EXISTING_FILES") return REMOVE_EXISTING_FILES;
    if (inString == "DELTAS")                return DELTAS;
    if (inString == "DOMAIN_SIZE") 			 return DOMAIN_SIZE;
	if (inString == "SOURCE")				 return SOURCE;
	if (inString == "STOP_SIMUL_AFTER") 	 return STOP_SIMUL_AFTER;
	if (inString == "TEMP_INIT") 			 return TEMP_INIT;
	if (inString == "TIME_STEP") 			 return TIME_STEP;
	if (inString == "BOUNDARY_CONDITIONS") 	 return BOUNDARY_CONDITIONS;
	if (inString == "OUTPUT_SAVING") 		 return OUTPUT_SAVING;
	if (inString == "MATERIALS") 			 return MATERIALS;
	if (inString == "ORIGINS")   			 return ORIGINS;
	if (inString == "PROBING_POINTS")        return PROBING_POINTS;
	if (inString == "ELECTRO_STEADY_STATE")  return ELECTRO_STEADY_STATE;
	else {
		printf("In file %s at %d. Complain to Romin. Abort().\n",__FILE__,__LINE__);
		cout << "Faulty string is ::" + inString + "::" << endl;
		abort();
	}
}

// Get lengths
double InputParser::get_length_WholeDomain(unsigned int direction,std::string type){
	if(strcmp(type.c_str(),"ELECTRO") == 0){
		// Direction = 0 gives along X, 1 along Y, 2 along Z:
		if(direction == 0){
			return this->lengthX_WholeDomain_Electro;
		}else if(direction == 1){
			return this->lengthY_WholeDomain_Electro;
		}else if(direction == 2){
			return this->lengthZ_WholeDomain_Electro;
		}else{
			printf("InputParser::get_length::ERROR:\n\tDirection should be between 0 and 2.");
			printf("Aborting (file %s at %d).\n",__FILE__,__LINE__);
			abort();
		}
	}else if(strcmp(type.c_str(),"THERMAL") == 0){
			// Direction = 0 gives along X, 1 along Y, 2 along Z:
		if(direction == 0){
			return this->lengthX_WholeDomain_Thermal;
		}else if(direction == 1){
			return this->lengthY_WholeDomain_Thermal;
		}else if(direction == 2){
			return this->lengthZ_WholeDomain_Thermal;
		}else{
			printf("InputParser::get_length::ERROR:\n\tDirection should be between 0 and 2.");
			printf("Aborting (file %s at %d).\n",__FILE__,__LINE__);
			abort();
		}
	}else{
		fprintf(stderr,"InputParser::get_length_WholeDomain::ERROR\n");
		fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
		abort();
	}
}

InputParser::InputParser(string file_name){

		this->filename = file_name;

}

void InputParser::defaultParsingFromFile(int MPI_RANK){
	if(this->filename == string()){
		fprintf(stderr,"In %s :: No file provided. Aborting.\n",
			__FUNCTION__);
		fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
		#ifdef MPI_COMM_WORLD
			MPI_Abort(MPI_COMM_WORLD,-1);
		#else
			abort();
		#endif
	}else{
		// Check that the file exist:
		if(this->is_file_exist(this->filename)){
			#if DEBUG > 1
			cout << "Input file exist !" << endl;
			#endif
		}else{
			cout << "InputParser::defaultParsingFromFile\n";
			cout << "Input file doesn't exist ! Aborting.\n";
			abort();
		}
		// Go to the parsing function:
		this->basicParsing(this->filename);
	}

	this->deleteFiles(MPI_RANK);
	
	if(this->source.number_of_sources.get() != 0 && this->conditionsInsideSources.empty()){
		fprintf(stderr,"In %s :: ERROR :: There is %d sources, but you specified"
						" no condition inside of it! Aborting.\n",__FUNCTION__,
						this->source.number_of_sources.get());
		fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
		#ifdef MPI_COMM_WORLD
			MPI_Abort(MPI_COMM_WORLD,-1);
		#else
			abort();
		#endif
	}
}

void InputParser::defaultParsingFromFile(std::string &filename, int MPI_RANK){
	// Check that the file exists:
	if(this->is_file_exist(filename)){
		/* Input file exists. */
	}else{
		fprintf(stderr,"%sIn %s :: No file provided/ not found / cannot open it (given file name is |%s|).\n"
					   " Have you tried '-inputfile my_file.input' ? Aborting.%s\n",
			ANSI_COLOR_RED,__FUNCTION__,filename.c_str(),ANSI_COLOR_RESET);
		fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
		#ifdef MPI_COMM_WORLD
			MPI_Abort(MPI_COMM_WORLD,-1);
		#else
			abort();
		#endif
	}
	// Go to the parsing function:
	this->basicParsing(filename);

	this->deleteFiles(MPI_RANK);
}

/**
 * @brief Checks if a file exists.
 */
bool InputParser::is_file_exist(const string &fileName){
	
    ifstream infile(fileName);
	
    return infile.good();
}

void InputParser::basicParsing(const string &filename){
	// Check the extension of the file:
	if(filename.substr(filename.find_last_of(".")+1) == "input"){
		// The extension is correct, proceed.
		// Open the file for reading:
		ifstream inputFile;
		inputFile.open(filename.c_str());
		if(inputFile.fail()){
			// Opening failed, aborting.
			fprintf(stderr,"In %s :: Cannot open input file %s. Aborting.\n",
				__FUNCTION__,filename.c_str());
			fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
			inputFile.clear();
			#ifdef MPI_COMM_WORLD
				MPI_Abort(MPI_COMM_WORLD,-1);
			#else
				abort();
			#endif
		}else if(inputFile.is_open()){
			// Contains the current read line of the input file:
			string currentLine;
			
			// Looping on the whole file:
			while(!inputFile.eof()){

				// Get line:
				getline(inputFile,currentLine);
				
				// Check that the line is not a comment:
				this->checkLineISNotComment(inputFile,currentLine);

				if(currentLine.find("$") != std::string::npos){
					// If there is a dollar, we begin a section:
					this->readHeader(inputFile,currentLine);
				}
			}
		}else{
			fprintf(stderr,"In %s :: Should not end up here. Aborting.\n",
				__FUNCTION__);
			fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
			inputFile.clear();
			#ifdef MPI_COMM_WORLD
				MPI_Abort(MPI_COMM_WORLD,-1);
			#else
				abort();
			#endif
		}
		inputFile.close();
	}else{
		fprintf(stderr,"In %s :: The input file is not under"
					" .input format (has %s). Please check your input file. Aborting.\n",
					__FUNCTION__,filename.c_str());
		fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
		#ifdef MPI_COMM_WORLD
			MPI_Abort(MPI_COMM_WORLD,-1);
		#else
			abort();
		#endif
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
	
	// Go with the switch:
	switch(hashit_Header1(strHeader1)){

		case INFOS    : 
			this->readHeader_INFOS(file);
			break;

		case MESH     : 
			this->readHeader_MESH (file);
			break;

		case RUN_INFOS: 
			this->readHeader_RUN_INFOS(file);
			break;

		case POST_PROCESSING:
			this->readHeader_POST_PROCESSING(file);
			break;

		default:
			fprintf(stderr,"In %s :: should not end up here. Aborting.\n",__FUNCTION__);
			fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
			#ifdef MPI_COMM_WORLD
				MPI_Abort(MPI_COMM_WORLD,-1);
			#else
				abort();
			#endif
	}
}

// Check that the line is not a comment:
bool InputParser::checkLineISNotComment(ifstream &file, string &currentLine){
	/* Check for line(s) being comments or blank */
	this->RemoveAnyBlankSpaceInStr(currentLine);
	while(currentLine == string() && !file.eof()){
		
		getline(file,currentLine);
		this->RemoveAnyBlankSpaceInStr(currentLine);

	}

	const std::string comment1_beg = "/*";
	const std::string comment1_end = "*/";
	const std::string comment2 = "//";

	if(currentLine.find(comment1_beg) != std::string::npos){
		// We have multiple lines comment ! Read all the comment before exiting.
		while(!file.eof()){
			if(currentLine.find(comment1_end) != std::string::npos){
				getline(file,currentLine);
				this->checkLineISNotComment(file,currentLine);
				break;
			}
			getline(file,currentLine);
		}
		return true;
	}else if(currentLine.find(comment2) != std::string::npos){
		// We have one line comment !
		while(!file.eof()){
			if(currentLine.find(comment2) == std::string::npos){
				this->checkLineISNotComment(file,currentLine);
				break;
			}
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

		if(currentLine == "INFOS"){break;}

		switch(this->hashit_Header2(currentLine)){
			case NAME:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					getline(file,currentLine);
					this->checkLineISNotComment(file,currentLine);
					this->RemoveAnyBlankSpaceInStr(currentLine);

					if(currentLine == "$NAME"){
						break;
					}

					if(currentLine == string()){continue;}

					std::size_t posEqual  = currentLine.find("=");
					std::string propName  = currentLine.substr(0,posEqual); 
					std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());
					
					if(propName != "output" && propName != "error" && propName != "profile"){
						printf("InputParser::readHeader_INFOS:: You didn't provide a ");
						printf("good member for $INFOS$NAME.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s:%d)\n",__FILE__,__LINE__);
						abort();
					}
					this->outputNames[propName] = propGiven;
				}
				break;

			case REMOVE_EXISTING_FILES:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					getline(file,currentLine);
					this->checkLineISNotComment(file,currentLine);
					this->RemoveAnyBlankSpaceInStr(currentLine);

					if(currentLine == "$REMOVE_EXISTING_FILES"){
						break;
					}
					if(currentLine == string()){continue;}
					std::size_t posEqual  = currentLine.find("=");
					std::string propName  = currentLine.substr(0,posEqual); 
					std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());
					
					if(propName == "remove_vti" || propName == "remove_pvti"){
						bool propGiven_bool = false;
						(propGiven == "true") ? propGiven_bool = true : propGiven_bool = false;

						this->removeWhat_dico.insert(
							 std::pair<std::string,bool>(propName,propGiven_bool)
						);

					}else{
						fprintf(stderr,"In %s :: REMOVE_EXISTING_FILES :: wrong property name. Aborting.\n",
							__FUNCTION__);
						fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
						#ifdef MPI_COMM_WORLD
						MPI_Abort(MPI_COMM_WORLD,-1);
						#else
						abort();
						#endif
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

		if(currentLine == "MESH"){break;}

		switch(this->hashit_Header2(currentLine)){
			case DELTAS:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$DELTAS"){
						
						// Check all electromagnetic deltas are equal !
						if(this->deltaX_Electro == this->deltaY_Electro
							&& this->deltaX_Electro == this->deltaZ_Electro
							&& this->deltaY_Electro == this->deltaZ_Electro)
						{
							// Check the ratio:
							double ratio = this->deltaX_Electro / this->delta_Thermal;

							// Use absolute value to keep rounding error away.

							if(abs(ratio - this->ratio_EM_TH_delta) > ratio*1E-8){
								fprintf(stderr,"InputParser::Wrong spatial steps\n");
								fprintf(stderr,"Spatial step EM is %f.",this->deltaX_Electro);
								fprintf(stderr,"Spatial step TH is %f.",this->delta_Thermal);
								fprintf(stderr,"Announced ratio is %f but has %f.\n",
									this->ratio_EM_TH_delta,ratio);
								fprintf(stderr,"Aborting.\nFile %s:%d\n",__FILE__,__LINE__);
								std::abort();
							}
						}else{
							fprintf(stderr,"InputParser::wrong electromagnetic time steps.\n");
							fprintf(stderr,"Spatial steps for EM grid should be equal.\nAborting.\n");
							fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
							std::abort();
						}
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

					// Spatial step delta along X for electromagnetic mesh:
					if(propName == "deltaX_Electro"){
						this->deltaX_Electro = std::stod(propGiven);

					// Spatial step delta along Y for electromagnetic mesh:
					}else if(propName == "deltaY_Electro"){
						this->deltaY_Electro = std::stod(propGiven);

					// Spatial step delta along Z for electromagnetic mesh:
					}else if(propName == "deltaZ_Electro"){
						this->deltaZ_Electro = std::stod(propGiven);

					// Spatial step delta along all three directions for the thermal mesh:
					}else if(propName == "delta_Thermal"){
						this->delta_Thermal = std::stod(propGiven);

					// Ratio between the spatial and thermal grids
					// If deltaElectro = 4 and deltaThermal = 2, the ratio is 0.5.
					}else if(propName == "ratio_EM_TH_delta"){
						this->ratio_EM_TH_delta = std::stod(propGiven);

					}else if(propName != "deltaX_electro" 
							&& propName != "deltaY_Electro" 
							&& propName != "deltaZ_Electro"
							&& propName != "delta_Thermal"){
						printf("InputParser::readHeader_MESH:: You didn't provide a ");
						printf("good member for $MESH$DELTAS. Has %s.\nAborting.\n",propName.c_str());
						cout << propName << endl;
						printf("(in file %s:%d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;
			
			
			case DOMAIN_SIZE:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);

					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$DOMAIN_SIZE"){
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

					if(propName == "L_X_ELECTRO"){
						this->lengthX_WholeDomain_Electro = std::stod(propGiven);

					}else if(propName == "L_Y_ELECTRO"){
						this->lengthY_WholeDomain_Electro = std::stod(propGiven);

					}else if(propName == "L_Z_ELECTRO"){
						this->lengthZ_WholeDomain_Electro = std::stod(propGiven);

					}else if(propName == "L_X_THERMAL"){
						this->lengthX_WholeDomain_Thermal = std::stod(propGiven);

					}else if(propName == "L_Y_THERMAL"){
						this->lengthY_WholeDomain_Thermal = std::stod(propGiven);

					}else if(propName == "L_Z_THERMAL"){
						this->lengthZ_WholeDomain_Thermal = std::stod(propGiven);
					
					}else if(propName != "L_X_ELECTRO" 
							&& propName != "L_Y_ELECTRO" 
							&& propName != "L_Z_ELECTRO"
							&& propName != "L_X_THERMAL"
							&& propName != "L_Y_THERMAL"
							&& propName != "L_Z_THERMAL"){
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
							if(this->conditionsInsideSources.empty()){
								fprintf(stderr,"In %s :: ERROR :: you have more than zero "
											"source but not imposed condition specified !"
											" Aborting.\n",
											__FUNCTION__);
								fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
								#ifdef MPI_COMM_WORLD
									MPI_Abort(MPI_COMM_WORLD,-1);
								#else
									abort();
								#endif
							}
							if(this->source_time.empty() || 
									this->source_time.size() != this->source.get_number_of_sources())
								DISPLAY_ERROR_ABORT(
									"You must specify source type ! (GAUSSIAN, SINE)"
								);
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

						if(propName != "NBR_SOURCES" && nbr_Sources_Defined == false){
							printf("InputParser::readHeader_MESH::CASE SOURCE\n");
							printf("You must first set the number of sources. Aborting.\n");
							std::abort();
						}

						if(propName == "NBR_SOURCES"){
							this->source.set_number_of_sources(std::stod(propGiven));
							nbr_Sources_Defined = true;

						}else if(propName == "L_X"){
							std::vector<double> temp = 
								this->determineVectorFromStr(
									propGiven,
									this->source.get_number_of_sources());
							this->source.setLengthAlongOneDir(0,temp);

						}else if(propName == "L_Y"){
							std::vector<double> temp = 
								this->determineVectorFromStr(
									propGiven,
									this->source.get_number_of_sources());
							this->source.setLengthAlongOneDir(1,temp);

						}else if(propName == "L_Z"){
							std::vector<double> temp = 
								this->determineVectorFromStr(
									propGiven,
									this->source.get_number_of_sources());
							this->source.setLengthAlongOneDir(2,temp);

						}else if(propName == "C_X"){
							std::vector<double> temp = 
								this->determineVectorFromStr(
									propGiven,
									this->source.get_number_of_sources());
							this->source.setCenterAlongOneDir(0,temp);

						}else if(propName == "C_Y"){
							std::vector<double> temp = 
								this->determineVectorFromStr(
									propGiven,
									this->source.get_number_of_sources());
							this->source.setCenterAlongOneDir(1,temp);

						}else if(propName == "C_Z"){
							std::vector<double> temp = 
								this->determineVectorFromStr(
									propGiven,
									this->source.get_number_of_sources());
							this->source.setCenterAlongOneDir(2,temp);

						}else if(propName == "FRQCY"){
							std::vector<double> temp = 
								this->determineVectorFromStr(
									propGiven,
									this->source.get_number_of_sources());
							this->source.setAllFrequencies(temp);
							
						}else if(propName == "IMPOSED"){
							/// Find the semi-colons:
							std::vector<size_t> pos_semi_col
								= findCharacterInsideString(propGiven, ";");
							if(pos_semi_col.size() != this->source.number_of_sources.get()){
								fprintf(stderr,"In %s :: ERROR :: IMPOSED :: "
												"You must have as many semi-colon(s) as "
												"you have sources. Aborting.\n",
												__FUNCTION__);
								fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
								#ifdef MPI_COMM_WORLD
									MPI_Abort(MPI_COMM_WORLD,-1);
								#else
									abort();
								#endif
							}
							
							for(unsigned int I = 0 ; I < pos_semi_col.size() ; I ++){
								std::string str;
								if(I == 0){
									str = propGiven.substr(
										0,pos_semi_col[I]);
								}else{
									str = propGiven.substr(
										pos_semi_col[I-1]+1,
										(pos_semi_col[I]-pos_semi_col[I-1])-1);
								}
								//printf("Has : %s \n",str.c_str());
								if(str != "DIPOLE" && str != "SIMPLE"){
									fprintf(stderr,"In %s :: ERROR :: Imposed condition on source"
												" should be either 'DIPOLE' or 'SIMPLE'. Aborting.\n",
												__FUNCTION__);
									fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
									#ifdef MPI_COMM_WORLD
										MPI_Abort(MPI_COMM_WORLD,-1);
									#else
										abort();
									#endif
								}
								this->conditionsInsideSources.push_back(str);
							}
							
						}else if(propName == "SOURCE_TIME"){
							/// Find the semi-colons:
							std::vector<size_t> pos_semi_col
								= findCharacterInsideString(propGiven, ";");
							if(pos_semi_col.size() != this->source.get_number_of_sources()){
								DISPLAY_ERROR_ABORT(
									"You must provided as many semi-colon(s)"
									" as there are source(s) (has %u source(s) and %s.",
									this->source.get_number_of_sources(),
									propGiven.c_str()
								);
							}
							for(size_t i = 0 ; i < pos_semi_col.size() ; i ++){
								size_t beg, length;
								if(i==0){beg = 0; length = pos_semi_col[0];}
								else{
									beg = pos_semi_col[i-1]+1;
									length = pos_semi_col[i]-pos_semi_col[i-1]-1;
								}
								std::string given = 
									propGiven.substr(beg,length);
								if(given != "GAUSSIAN" && given != "SINE"){
									DISPLAY_ERROR_ABORT(
										"The given property is different from GAUSSIAN or SINE (has %s).",
										given.c_str()
									);
								}
								this->source_time.push_back(given);
							}
							

						}else{
							printf("InputParser::readHeader_MESH:: You didn't provide a ");
							printf("good member for $MESH$SOURCE.\nAborting.\n");
							cout << propName << endl;
							printf("(in file %s at %d)\n",__FILE__,__LINE__);
							abort();
						}
					}
				}
				break;

			case MATERIALS:

				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$MATERIALS"){
						if(!is_directory(this->material_data_directory)){
							DISPLAY_ERROR_ABORT(
								"You mus provide a valid material data directory. Has %s.",
								this->material_data_directory.string().c_str()
							);
						}
						if(this->simulationType.get() == "USE_GEOMETRY_FILE"){
							if(this->file_containing_geometry == std::string()){
								DISPLAY_ERROR_ABORT(
									"You want to use a geometry file, but you don't provide it."
								);
							}
						}
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

					if(propName == "USE_AIR_EVERYWHERE"){
						if(propGiven == "true"){
							this->simulationType = "USE_AIR_EVERYWHERE";
						}else{
							printf("You set USE_AIR_EVERYWHERE");
							printf(" to false.Aborting.\n");
							std::abort();
						}
					
					}else if(propName == "MATERIAL_DIRECTORY"){
						boost::filesystem::path givenPath(propGiven);
						if(is_directory(givenPath)){
							this->material_data_directory = givenPath;
						}else{
							int ret = search_for_directory(
								propGiven,
								givenPath
							);
							if(ret != EXIT_SUCCESS){
								DISPLAY_ERROR_ABORT(
									"Cannot find directory %s.",
									propGiven.c_str()
								);
							}
							this->material_data_directory = givenPath;
						}
						//std::cout << "Material dirctory is " << this->material_data_directory << std::endl;

					}else if(propName == "USE_GEOMETRY_FILE"){
						// Parse the propGiven string:
						std::vector<size_t> pos_acc_open
							= findCharacterInsideString(propGiven,"{");
						std::vector<size_t> pos_acc_closed
							= findCharacterInsideString(propGiven,"}");
						std::vector<size_t> pos_comma
							= findCharacterInsideString(propGiven,",");

						if(    pos_acc_open.size()   != 1
							|| pos_acc_closed.size() != 1
							|| pos_comma.size()      != 1){
								std::string compl_info = std::string();
								if(pos_comma.size() != 1){
									compl_info = 
										"You possibly put a semi-colon separator"
										" instead of a comma separator.";
								}else{
									compl_info =
										"You probably miss a brace.";
								}
								DISPLAY_ERROR_ABORT(
									"You should provide something like"
									" USE_GEOMETRY_FILE={true,my_file.geomtry}"
									" (has %s).%s",
									propGiven.c_str(),
									compl_info.c_str()
								);
							}
						/// Check that the boolean is set to true:
						std::string boolVar = 
							propGiven.substr(
								pos_acc_open[0]+1,
								pos_comma[0]-pos_acc_open[0]-1
							);
						if(boolVar != "true"){
							DISPLAY_ERROR_ABORT(
								"You want USE_GEOMETRY_FILE but you set it to %s."
								" Please set it to true.",
								boolVar.c_str()
							);
						}else{
							this->simulationType = propName;
						}
						/// Get the geometry file name:
						std::string geometryFilename
							= propGiven.substr(
								pos_comma[0]+1,
								pos_acc_closed[0]-pos_comma[0]-1
							);
						/*if(geometryFilename.substr(geometryFilename.find('.')+1) != "geometry"){*/
						if(boost::filesystem::extension(geometryFilename) != ".geometry"){
							DISPLAY_ERROR_ABORT(
								"The geometry file should have the extension '.geometry'"
								" but it has '.%s'.",
								boost::filesystem::extension(geometryFilename).c_str()
							);
						}else{
							this->file_containing_geometry
								= geometryFilename;
						}
						
						

					}else if(propName == "TEST_PARAVIEW"){
						if(propGiven == "true"){
							this->simulationType = "TEST_PARAVIEW";
						}else{
							printf("You set TEST_PARAVIEW to false. Aborting.\n");
							std::abort();
						}

					}else if(propName == "TEST_PARAVIEW_MPI"){
						// Assign the simulation type:
						this->simulationType = "TEST_PARAVIEW_MPI";
						// Parse the propGiven:
						/**
						 * The 'propGiven' must be (A,B,C), i.e. comma-separated
						 * and surrounded with brackets.
						 */
						/// Remove the brackets:
						if(propGiven.c_str()[0] == '('){
							/// Check the last character is a bracket:
							if(propGiven.c_str()[propGiven.size()-1] == ')'){
								/// Erase the two parenthesis from the 'propGiven' string:
								propGiven = propGiven.substr(1, propGiven.size() -2);
							}else{
								fprintf(stderr,"In %s :: while parsing the arguments of 'TEST_PARAVIEW_MPI.\n",
									__FUNCTION__);
								fprintf(stderr,"You miss the closing parenthesis !\n");
								fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
								#ifdef MPI_COMM_WORLD
								MPI_Abort(MPI_COMM_WORLD,-1);
								#else
								abort();
								#endif
							}
						}else{
							fprintf(stderr,"In %s :: while parsing the arguments of 'TEST_PARAVIEW_MPI.\n",
									__FUNCTION__);
							fprintf(stderr,"You miss the opening parenthesis !\n");
							fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
							#ifdef MPI_COMM_WORLD
							MPI_Abort(MPI_COMM_WORLD,-1);
							#else
							abort();
							#endif
						}
						std::vector<std::string> args;
						std::stringstream ss(propGiven);
						while(ss.good()){
							std::string substr;
							getline(ss,substr,',');
							args.push_back(substr);
						}
						/// If there is more than 3 args, abort:
						if(args.size() > 3){
							fprintf(stderr,"In %s :: while parsing the arguments of 'TEST_PARAVIEW_MPI.\n",
									__FUNCTION__);
							fprintf(stderr,"You must given 3 args: TEMP=sthg,E=sthgElse,H=sthg !\n");
							fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
							#ifdef MPI_COMM_WORLD
							MPI_Abort(MPI_COMM_WORLD,-1);
							#else
							abort();
							#endif
						}
						std::vector<std::string> ARGS_ORDER = {"TEMP","E","H"};
						for(size_t i = 0 ; i < 3 ; i ++ ){
							if(args[i].find(ARGS_ORDER[i]) == std::string::npos){
								/// Cannot find the argument TEMP that should come first. Aborting.
								fprintf(stderr,"In %s :: while parsing the arguments of 'TEST_PARAVIEW_MPI.\n",
										__FUNCTION__);
								fprintf(stderr,"You must provide %s !\n",ARGS_ORDER[i].c_str());
								fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
								#ifdef MPI_COMM_WORLD
								MPI_Abort(MPI_COMM_WORLD,-1);
								#else
								abort();
								#endif
							}else{
								std::size_t posEqual   = args[i].find("=");
								std::string prop_name  = args[i].substr(0,posEqual);
								if(prop_name != ARGS_ORDER[i]){
									fprintf(stderr,"FATAL ERROR\n");
									fprintf(stderr,"FILE %s:%d\n",__FILE__,__LINE__);
									abort();
								}
								std::string prop_given = args[i].substr(posEqual+1,args[i].length());
								this->TEST_PARAVIEW_MPI_ARGS.insert(
									std::pair<std::string,std::string>(
										prop_name,
										prop_given
									));
							}
						}

						/* END OF TEST_PARAVIEW_MPI */

					}else if(propName == "MATERIAL_DATA_FILE"){
						// Determine the number of files:
						std::vector<size_t> pos_commas
							= findCharacterInsideString(propGiven,",");
						size_t nbr_files = pos_commas.size()+1;
						// Find equal signs:
						std::vector<size_t> pos_equal_signs
							= findCharacterInsideString(propGiven,"=");
						if(pos_equal_signs.size() != nbr_files){
							DISPLAY_ERROR_ABORT(
								"You have requested %zu material files but"
								" I only have %zu equal signs (should have %zu).",
								nbr_files,pos_equal_signs.size(),nbr_files
							);
						}

						size_t pos_acc_close = findCharacterInsideString(propGiven,"}")[0];
						
						for(size_t I = 0 ; I < nbr_files ; I++){
							size_t beg = 0;
							size_t size= 0;
							std::string type;
							if(I == 0){
								beg  = pos_equal_signs[0]+1;
								size = pos_commas[0]-beg;
								type = propGiven.substr(1,beg-2);
								//std::cout << type << std::endl;
							}else if( I == nbr_files - 1){
								beg = pos_equal_signs[nbr_files-1]+1;
								size= pos_acc_close-beg;
								type = propGiven.substr(
											pos_commas[nbr_files-2]+1,
											pos_equal_signs[nbr_files-1]-pos_commas[nbr_files-2]-1);
								//std::cout << type << std::endl;
							}else{
								beg  = pos_equal_signs[I-1]+1;
								size = pos_commas[I]-beg;
								type = propGiven.substr(
										pos_commas[I-1],
										pos_equal_signs[I-1]-pos_commas[I-1]-1);
								//std::cout << type << std::endl;
							}
							if(I == 0 && type != "OLD"){
								DISPLAY_ERROR_ABORT(
									"First provided file should be 'OLD'."
								);
							}
							if(I == 1 && type != "ELECTRO"){
								DISPLAY_ERROR_ABORT(
									"Second provided file should be 'ELECTRO'."
								);
							}
							std::string whichFile
								= propGiven.substr(beg,size);
							//std::cout << whichFile << std::endl;
							this->material_data_files.push_back(whichFile);
						}

					}else{
						if(propName == "MATERIAL_DATA_FILE")
							continue;
						printf("InputParser::readHeader_MESH:: You didn't provide a ");
						printf("good member for $MESH$MATERIALS (has %s).\nAborting.\n",propName.c_str());
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				/* Check that simulation type was effectively set: */
				if(this->simulationType.get_alreadySet() == false){
					printf("InputParser::readHeader_MESH::MATERIALS\n");
					printf("Exiting section MATERIALS without specifying");
					printf(" anything. Aborting().\n\n");
					abort();
				}
				break;

			case ORIGINS:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$ORIGINS"){						
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

					if(propName == "ORIGIN_ELECTRO_X"){
						this->origin_Electro_grid[0]= std::stod(propGiven);
						
					}else if(propName == "ORIGIN_ELECTRO_Y"){
						this->origin_Electro_grid[1] = std::stod(propGiven);
					
					}else if(propName == "ORIGIN_ELECTRO_Z"){
						this->origin_Electro_grid[2] = std::stod(propGiven);

					}else if(propName == "ORIGIN_THERMAL_X"){
						this->origin_Thermal_grid[0] = std::stod(propGiven);

					}else if(propName == "ORIGIN_THERMAL_Y"){
						this->origin_Thermal_grid[1] = std::stod(propGiven);

					}else if(propName == "ORIGIN_THERMAL_Z"){
						this->origin_Thermal_grid[2] = std::stod(propGiven);

					}else{
						fprintf(stderr,"In %s (%d) :: Cannot associate %s to any property. Aborting.\n",
							__FUNCTION__,__LINE__,propName.c_str());
						fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
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

		if(currentLine == "RUN_INFOS"){
			break;
		}
		switch(this->hashit_Header2(currentLine)){

			case STOP_SIMUL_AFTER:

				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$STOP_SIMUL_AFTER"){
						bool do_abort = false;
						/// Check number of steps have been set different from zero:
						if(this->maxStepsForOneCycleOfThermal == 0){
							fprintf(stderr,"In %s :: ERROR :: max number of steps for the thermal"
										" algorithm is equal to zero. Aborting.\n",
										__FUNCTION__);
							do_abort = true;
						}
						if(this->maxStepsForOneCycleOfElectro == 0){
							fprintf(stderr,"In %s :: ERROR :: max number of steps for the electro"
										" algorithm is equal to zero. Aborting.\n",
										__FUNCTION__);
							do_abort = true;
						}
						if(do_abort){
							fprintf(stderr,"In %s ::  Aborting (because of a number of steps set to zero).\n",
								__FUNCTION__);
							fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
							#ifdef MPI_COMM_WORLD
							MPI_Abort(MPI_COMM_WORLD,-1);
							#else
							abort();
							#endif
						}
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

					if(propName == "stopTime"){
						this->stopTime = std::stod(propGiven);

					}else if(propName == "maxStepsForOneCycleOfElectro"){
						/**
						 * This property is usefull to impose a maximum number of steps
						 * for the electromagnetic solver before stopping.
						 */
						/// Transform the string in a size_t with std::stold and a cast:
						this->maxStepsForOneCycleOfElectro = (size_t) std::stold(propGiven);;
						
						if(this->maxStepsForOneCycleOfElectro == 0){
							fprintf(stderr,"In %s ::maxStepsForOneCycleOfElectro has been set to zero. Aborting.\n",
								__FUNCTION__);
							fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
							#ifdef MPI_COMM_WORLD
							MPI_Abort(MPI_COMM_WORLD,-1);
							#else
							abort();
							#endif
						}
					}else if(propName == "maxStepsForOneCycleOfThermal"){
						/**
						 * This property is usefull to impose a maximum number of steps
						 * for the electromagnetic solver before stopping.
						 */
						/// Transform the string in a size_t with std::stold and a cast:
						this->maxStepsForOneCycleOfThermal = (size_t) std::stold(propGiven);
						/*printf("max time step thermal %zu\n\n",
							this->maxStepsForOneCycleOfThermal);*/
						
						if(this->maxStepsForOneCycleOfThermal == 0){
							fprintf(stderr,"In %s ::maxStepsForOneCycleOfThermal has been set to zero. Aborting.\n",
								__FUNCTION__);
							fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
							#ifdef MPI_COMM_WORLD
							MPI_Abort(MPI_COMM_WORLD,-1);
							#else
							abort();
							#endif
						}

					}else{
						printf("InputParser::readHeader_RUN_INFOS:: You didn't provide a ");
						printf("good member for $RUN_INFOS$STOP_SIMUL_AFTER.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;

			case TEMP_INIT:

				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$TEMP_INIT"){
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

					if(propName == "INIT_TEMP_FILE"){

						std::string init_temp_filename = propGiven;
						bool failed = true;

						if(this->is_file_exist(init_temp_filename)){
							failed = false;
							rapidjson::Document doc;
							read_json(init_temp_filename,doc);

							// Go through file:
							static const char* kTypeNames[] = 
    							{ "Null", "False", "True", "Object", "Array", "String", "Number" };
							for (Value::ConstMemberIterator itr = doc.MemberBegin();
								itr != doc.MemberEnd(); ++itr)
							{
								std::string name(itr->name.GetString());
								std::string value(kTypeNames[itr->value.GetType()]);
								// The file contains only .info and .initTemp fields:
								if(name.find(".initTemp") == std::string::npos
									&& name.find(".info") == std::string::npos){
										DISPLAY_ERROR_ABORT(
											"In file %s: the name %s doesn't contain '.initTemp'.",
											init_temp_filename.c_str(),
											name.c_str()
										);
								}
								if(name.find(".info") != std::string::npos)
									continue;
								if( value != "Number"){
									DISPLAY_ERROR_ABORT(
										"In file %s: the name %s has not a number value"
										" but a %s value.",
										init_temp_filename.c_str(),
										name.c_str(),
										kTypeNames[itr->value.GetType()]
									);
								}
								std::string mat = name.substr(0,name.find(".initTemp"));
								double tempInit = itr->value.GetDouble();
								//std::cout << "Adding material initial temperature:: " + mat << tempInit << std::endl;
								this->GetInitTemp_FromMaterialName.insert(
									std::pair<std::string,double>(mat,tempInit)
									); 

							}
							
						}
						if(failed){
							DISPLAY_ERROR_ABORT(
								"File %s is not found.",init_temp_filename.c_str()
							);
						}

					}else{
						printf("InputParser::readHeader_RUN_INFOS:: You didn't provide a ");
						printf("good member for $RUN_INFOS$TEMP_INIT.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;

			case TIME_STEP:
				// Read the given time steps:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$TIME_STEP"){
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

					if(propName == "THERMAL_TIME_STEP"){
						this->thermal_algo_time_step = std::stod(propGiven);
						
					}else{
						printf("InputParser::readHeader_RUN_INFOS:: You didn't provide a ");
						printf("good member for $RUN_INFOS$TIME_STEP.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;

			case ELECTRO_STEADY_STATE:
				// Read the given time steps:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$ELECTRO_STEADY_STATE"){
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

					if(propName == "CHECK_EVERY_POINT"){

					}else{
						DISPLAY_ERROR_ABORT(
							"In ELECTRO_STEADY_STATE :: nothing corresponds to %s.",
							propName.c_str()
						);
					}
				}
				break;

			case OUTPUT_SAVING:
				// Read the given time steps:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$OUTPUT_SAVING"){
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

					if(propName == "SAMPLING_FREQ_ELECTRO"){
						this->SAMPLING_FREQ_ELECTRO = std::stol(propGiven);

					}else if(propName == "SAMPLING_FREQ_THERMAL"){
						this->SAMPLING_FREQ_THERMAL = std::stol(propGiven);

					}else{
						printf("InputParser::readHeader_RUN_INFOS:: You didn't provide a ");
						printf("good member for $RUN_INFOS$TEMP_INIT.\nAborting.\n");
						cout << propName << endl;
						printf("(in file %s at %d)\n",__FILE__,__LINE__);
						abort();
					}
				}
				break;

			
			case BOUNDARY_CONDITIONS:
				// Read the given boundary conditions:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					// Read line:
					getline(file,currentLine);
					// Get rid of comments:
					this->checkLineISNotComment(file,currentLine);
					// Remove any blank in the string:
					this->RemoveAnyBlankSpaceInStr(currentLine);
					// If the string is "$DELTAS" it means the section ends.
					if(currentLine == "$BOUNDARY_CONDITIONS"){
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

					if(propName == "BC_FACE_0"
						|| propName == "BC_FACE_1"
						|| propName == "BC_FACE_2"
						|| propName == "BC_FACE_3"
						|| propName == "BC_FACE_4"
						|| propName == "BC_FACE_5")
					{
						// Find the position of the second underscore:
						size_t pos = propName.find('_', 0);
						pos        = propName.find('_',pos+1);
						size_t numberFace = std::stol(propName.substr(pos+1,propName.size()));
						
						// Parse the right hand side:
						
						// Detect first and last brackets:
						size_t pos_bra_in  = propGiven.find("{");
						size_t pos_bra_out = propGiven.find("}");
						size_t pos_comma   = propGiven.find(";");
						
						std::string condition = propGiven.substr(pos_bra_in+1,(pos_comma-pos_bra_in-1));
						std::string value     = propGiven.substr(pos_comma+1,(pos_bra_out-pos_comma-1));

						double value_dbl = std::stod(value);

						this->THERMAL_FACE_BC_TYPE.insert(std::pair<std::size_t,std::string>(numberFace,condition));

						this->THERMAL_FACE_BC_VALUE.insert(std::pair<std::size_t,double>(numberFace,value_dbl));
						
						/*printf("Face %d ! has |%s| => |%s| and |%s| (%zu,%zu,%zu) --> (%s,%lf)\n",(int)numberFace,
							propGiven.c_str(),condition.c_str(),value.c_str(),
							pos_bra_in,pos_comma,pos_bra_out,
							this->THERMAL_FACE_BC_TYPE[numberFace].c_str(),
							this->THERMAL_FACE_BC_VALUE[numberFace]);*/
						
						
											
					}else{
						printf("InputParser::readHeader_RUN_INFOS:: You didn't provide a ");
						printf("good member for $RUN_INFOS$BOUNDARY_CONDITIONS.\nAborting.\n");
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

std::vector<double> InputParser::determineVectorFromStr(
		std::string str,
		size_t size_to_verify_for /* = 0 */
	)
{
	std::stringstream stream(str);
	std::string word;
	std::vector<double> tempVec;
	while( getline(stream, word, ';') ){
		tempVec.push_back(std::stod(word));
	}

	if(tempVec.size() > size_to_verify_for){
		fprintf(stderr,"In %s :: The vector from the string contains more elements than annouced.\n",
			__FUNCTION__);
		fprintf(stderr,"Received %s, announced %zu elements but has %zu elements. Aborting.\n",
			str.c_str(),size_to_verify_for,tempVec.size());
		abort();
	}
		
	return tempVec;
}

void InputParser::readHeader_POST_PROCESSING(ifstream &file){
	
	std::string currentLine = string();
	
	while(currentLine != "POST_PROCESSING"){
		// Read line:
		getline(file,currentLine);
		// Get rid of comments:
		this->checkLineISNotComment(file,currentLine);
		// Remove any blank space:
		this->RemoveAnyBlankSpaceInStr(currentLine);
		// Remove Dollar sign:
		currentLine = currentLine.substr(currentLine.find("$")+1);

		if(currentLine == "POST_PROCESSING"){
			break;
		}

		switch(this->hashit_Header2(currentLine)){
			case PROBING_POINTS:
				while(!file.eof()){
					// Note: sections are ended by $the-section-name.
					getline(file,currentLine);
					this->checkLineISNotComment(file,currentLine);
					this->RemoveAnyBlankSpaceInStr(currentLine);

					if(currentLine == "$PROBING_POINTS"){
						break;
					}

					if(currentLine == string()){continue;}

					std::size_t posEqual  = currentLine.find("=");
					std::string propName  = currentLine.substr(0,posEqual); 
					std::string propGiven = currentLine.substr(posEqual+1,currentLine.length());
					
					if( propName == "probe_point"){
						// More than on point can be probed:
						std::vector<size_t> pos_commas
							= findCharacterInsideString(propGiven,",");
						std::vector<size_t> pos_accol_open
							= findCharacterInsideString(propGiven,"{");
						std::vector<size_t> pos_accol_close
							= findCharacterInsideString(propGiven,"}");

						if(    pos_commas.size()      != 4 
							|| pos_accol_open.size()  != 1
							|| pos_accol_close.size() != 1)
						{
							DISPLAY_ERROR_ABORT(
								"probe_point :: Wrong input (has %s)"
								" but expected is something like"
								" {Ex,0.2,0.2,0.2,ALL}.",
								propGiven.c_str()
							);
						}

						std::string type_field = 
							propGiven.substr(pos_accol_open[0]+1,
								(pos_commas[0]-pos_accol_open[0])-1);
						std::vector<double> coord(3);
						coord[0] = std::stod(propGiven.substr(pos_commas[0]+1,
										(pos_commas[1]-pos_commas[0])+1));
						coord[1] = std::stod(propGiven.substr(pos_commas[1]+1,
										(pos_commas[2]-pos_commas[1])+1));
						coord[2] = std::stod(propGiven.substr(pos_commas[2]+1,
										(pos_commas[3]-pos_commas[2])+1));
						std::string at_which_time
							= propGiven.substr(pos_commas[3]+1,
										(pos_accol_close[0]-pos_commas[3])-1);



						std::string filename_ = "probe_point/";
						filename_.append(type_field);
						filename_.append("_");
						filename_.append(to_string(coord[0]));
						filename_.append("_");
						filename_.append(to_string(coord[1]));
						filename_.append("_");
						filename_.append(to_string(coord[2]));
						filename_.append("_");
						filename_.append(at_which_time);
						filename_.append(".txt");

						/*printf("{|%s|,%lf,%lf,%lf,|%s|} --> %s\n",
							type_field.c_str(),
							coord[0],coord[1],coord[2],
							at_which_time.c_str(),
							filename_.c_str());*/

						probed_point temp = {
							type_field,   //.type_field    
							coord,        //.coordinates 
							at_which_time,//.at_which_time
							filename_     //.filename
						};
						this->points_to_be_probed.push_back(temp);
						if(this->MPI_rank == 0){
							const std::string dir = "probe_point";
							directory_exists(dir,true);

							/// Check that no file with the same name exists. If there is one, delete it.
							if(is_file_exist(filename_)){
								remove(filename_.c_str());
							}
							/// Create the file:
							std::ofstream outfile (filename_,std::ofstream::out);
							if(!outfile.is_open()){
								DISPLAY_ERROR_ABORT(
									"Cannot create file %s.",filename_.c_str()
								);
							}
							char buff[DTTMSZ];
							outfile << "Created on " << getDtTm (buff);
							outfile << " | contains the field " + type_field;
							outfile << " at time(" + at_which_time + ")";
							outfile << " and at point (" << coord[0];
							outfile << "," << coord[1] << "," << coord[2] << ")" << std::endl;
							outfile.close();
						}
					}else if(propName == "probe_line"){
						/**
						 * @brief Parse the input file when the user wants to probe lines.
						 */
						// Example: probe_line={Ex,x=1,y=ALL,z=1}
						// This will probe Ex on the line (x=1,z=1) and y varying.
						DISPLAY_ERROR_ABORT("Not yet implemented.");
					}else{
						DISPLAY_ERROR_ABORT(
							"In $PROBING_POINTS :: no property corresponds to %s.",
							propName.c_str()
						);
					}
				}
				break;

			default:
				DISPLAY_ERROR_ABORT(
					"Should not end up here. Faulty line is %s.",
					currentLine.c_str()
				);
		}
	}
}


bool directory_exists( const std::string &directory, bool createIt )
{
    DIR* dir;

	if(NULL == (dir = opendir(directory.c_str()))){
		printf("Directory %s doesn't exist.\n",directory.c_str());
		if(createIt){
			/// Create the directory.
			int success;
			#ifdef __linux__
				success = mkdir(directory.c_str(), 0777); 
			#else
				success = _mkdir(directory.c_str());
			#endif
			if(success != 0){
				/// Directory creation failed:
				DISPLAY_ERROR_ABORT(
					"MKDIR(%s) failed.",directory.c_str()
				);
			}else{
				printf("Directory %s successfully created !\n",directory.c_str());
			}
		}
		return false;
	}

    return true;
}