#include "Materials.h"
#include <fstream>
#include <cstring>
#include <algorithm>

#include "CSV_parser.hpp"
#include "header_with_all_defines.hpp"

void Materials::get_properties_from_file_ELECTRO(std::string const &filename_ELECTRO_PROPS){
	/// Check that the file is a CSV file:
	std::string filename = filename_ELECTRO_PROPS;
	if(filename.substr(filename.find(".")+1) != "csv"){
		DISPLAY_ERROR_ABORT("The property file %s has not a '.csv' extension.",filename.c_str());
	}

	/// Read the whole CSV file:
	std::vector<std::vector<std::string> > data = parse2DCsvFile(filename);
 
	/// Count the number of different materials:
	size_t nbr_mat = 0;
	std::string current_mat = data[1][0];
	std::vector<std::string> list_mat;
    for (size_t I = 1 ; I < data.size() ;  I++) {
		std::vector<std::string> l = data[I];
        if(l[0] != current_mat){
			/// Check that the material was not already encountered:
			std::vector<std::string>::iterator it;
			it = std::find (list_mat.begin(), list_mat.end(), l[0]);
			if (it == list_mat.end()){
				//std::cout << "New mat : " << l[0] << std::endl;
				this->materialID_FromMaterialName_ELECTRO[l[0]]    = nbr_mat;
				this->materialName_FromMaterialID_ELECTRO[nbr_mat] = l[0];
				if(nbr_mat == 0 && l[0] != "Air"){
					DISPLAY_ERROR_ABORT(
						"The material with ID 0 MUST be air but as %s instead.",
						this->materialName_FromMaterialID_ELECTRO[nbr_mat].c_str()
					);
				}
				
				nbr_mat++;
				list_mat.push_back(l[0]);
				current_mat = l[0];
			}
		}
    }

	/// Count the number of properties:
	size_t nbr_prop = data[1].size()-1;

	/// Get the properties names (start for loop at 1 becase first column is material name):
	std::vector<std::string> propNames(nbr_prop);
	for(size_t I = 1 ; I <= nbr_prop ; I ++){
		propNames[I-1] = data[1][I];
		//std::cout << "New prop " + propNames[I-1] << std::endl;
	} 

	this->list_of_properties_of_list_of_materials_ELECTRO
		= propNames;

	/// Fill in this->list_of_materials:
	this->list_of_materials_ELECTRO = std::vector<material_struct>(nbr_mat);
	for(size_t I = 0 ; I < nbr_mat ; I ++){
		this->list_of_materials_ELECTRO[I].name = list_mat[I];
		this->list_of_materials_ELECTRO[I].ID   = I;
		for(size_t J = 0 ; J < nbr_prop ; J ++){
			this->list_of_materials_ELECTRO[I].properties.insert(
				std::pair<std::string,double>(propNames[J],0.0)
			);
		}
	}

	for(size_t I = 1 ; I < data.size() ; I ++){
		std::vector<std::string> l = data[I];

		std::vector<std::string>::iterator it;
		it = std::find (list_mat.begin(), list_mat.end(), l[0]);

		if (it != list_mat.end()){
			std::string mat = l[0];
			size_t ID = this->materialID_FromMaterialName_ELECTRO[mat];

			//printf("Current mat %s :: Id %zu\n",mat.c_str(),ID);
			for(size_t J = 1 ; J <= nbr_prop ; J ++){
				/// Identify the property:
				std::string propName = propNames[J-1];
				//std::cout << propName + "=" + l[J] << std::endl;

				double propVal = 0.0;

				if(l[J] == "\"N/A\""){
					propVal = nan("");
					//printf(">>>> Has %lf...\n",propVal);
				}else{
					propVal = std::stod(l[J]);
				}
				
				this->list_of_materials_ELECTRO[ID].properties[propName] = propVal;
			}
		}
	}

	#ifndef NDEBUG
	std::vector<material_struct> list_mat_prop = this->list_of_materials_ELECTRO;
	for(size_t I = 0 ; I < list_mat_prop.size() ; I ++){
		printf(">>> EM properties of material %s:\n",list_mat_prop[I].name.c_str());
		map<std::string, double>::iterator it;
		for(it = list_mat_prop[I].properties.begin(); it != list_mat_prop[I].properties.end() ; it ++){
			printf("\t> %s = %lf\n",
				it->first.c_str(),
				it->second);
		}
	}
	#endif
}

/****************************************/
/* Conventions for the properties file  */
/* This is an example:
 * Material,Temperature,Conductivity,Permittivity
 * AIR,25,2,5
 * AIR,35,56,8
 * WATER,56,5,9
 * WATER,85,4,9
*/
/* PLEASE NOTE THAT TEMPERATURE MUST BE GIVEN IN KELVIN */

/* You cannot mix water and air properties in the above example. You must first provide all air
 * properties and then all water properties, or the contrary. No mixing.*/

void Materials::getPropertiesFromFile(string filename,string filename_ELECTRO_PROPS){
	this->get_properties_from_file_ELECTRO(filename_ELECTRO_PROPS);
	/* Opening the file */
	ifstream file;
	file.open(filename,fstream::in);
	if(file.fail())
	{
		/* Opening failed - aborting ... */
		printf("There was an error while opening the file containing the properties '%s' (line %d).\n",
						filename.c_str(),__LINE__);
		printf("Check the access path.\n");
		abort();
	}
	////////////////////
	/* Acquiring data */
	////////////////////
	
	/* Will contain the current line read in the file */
	string currentLine;
	/* Will contain the current and the previous material - useful for counting the number of materials
		and building 3D array containing all the properties */
	string currentMaterial, previousMaterial;
	/* The first line of the file is not useful */
	getline(file,currentLine);

	/* Determine the number of information per material and the number of different materials */
	// Also, determine the number of properties.
	unsigned int maxNumberOfTemp = 0, counter = 1, numberOfMaterials = 0;
	// Note: numberOfMaterials starts at -1 because at the first comparaison between 
	// 	 current- and previousMaterial, it will always trigger "numberOfMaterials++"
	//	 even if it is not necessary.
	bool numberOfPropertiesIsKnown = false;
	while(!file.eof())
	{
		// Get the current line:
		getline(file,currentLine);
		#ifndef NDEBUG
			std::cout << currentLine << std::endl;
		#endif
		if(currentLine.empty()){
			// Do nothing, the line is empty !
		}else{
			// Get the position of the first comma:
			bool materialName = false;
			for(unsigned int i = 0 ; i < currentLine.length() ; i++){
					if(currentLine[i] == ','){
						// The current material is for example AIR,
						if(materialName == false){
							currentMaterial.assign(currentLine,0,i);
							materialName = true;
						}
						if(numberOfPropertiesIsKnown == true)
							break;
						else
							this->numberOfProperties ++;
					}
			}
			numberOfPropertiesIsKnown = true;
			// Compare current and previous materials to see if it has changed:
			if(currentMaterial.compare(previousMaterial) != 0){
				// Increment the number of materials
				numberOfMaterials++;

				//cout << "C::" + currentMaterial << "/P::" + previousMaterial << endl;

				// We want to track the material for which we have the highest number of 
				// properties specifications - it will be useful to initialize the table with
				// all the properties.
				counter ++;
				if(maxNumberOfTemp < counter){
					maxNumberOfTemp = counter;
				}
				counter = 1;
			}else{
				counter++;
				//cout << "Counter is " << counter << endl;
			}
			previousMaterial = currentMaterial;
		}
	}
	if(maxNumberOfTemp < counter){
		maxNumberOfTemp = counter;
	}
	// Set the oject fields:
	this->numberOfMaterials  = numberOfMaterials;
	this->numberOfProperties = numberOfProperties;
	this->maxNumberOfTemp    = maxNumberOfTemp; 	

	//cout << "maxNumberOfTemp = " << this->maxNumberOfTemp << endl;

	// Allocate the this->properties 3d array.
	/*this->properties = new double**[numberOfMaterials];
	for(unsigned int I = 0 ; I < numberOfMaterials ; I ++){
		this->properties[I] = new double*[numberOfProperties];
		for(unsigned int J = 0 ; J < numberOfProperties ; J ++){
			this->properties[I][J] = new double[maxNumberOfTemp];
			for(unsigned int K = 0 ; K < maxNumberOfTemp ; K ++){
				this->properties[I][J][K] = -1.;
			}
		}
	}*/
	
	this->properties.set_size_data(this->numberOfMaterials,
								   this->numberOfProperties,
								   this->maxNumberOfTemp);
	this->properties.fillIn(-1.0);
	//this->printAllProperties();
	
	// Empty the previous material variable for later.
	previousMaterial = string();
	/*
		cout << "The highest number of specification is " << maxNumberOfTemp << "." << endl;
		cout << "Number of different materials is " << numberOfMaterials << "." << endl;
		cout << "Number of properties is " << this->numberOfProperties << endl;
	*/
	// Rewind the position in the file for later:	
	file.clear();
	file.seekg(0);
	// The first line is not useful:	
	getline(file,currentLine);
	// Useful for later:
	string valueStr = string();
	// Go over the file again to collect the properties and store them:
	unsigned int counterMaterial = 0, counterProperties = 0, counterTemp = 0;	
	bool firstTime = true;
	while(!file.eof())
        {
		// Get line:
        getline(file,currentLine);
		bool MATERIAL = true;
		// Go over the whole line:
		for(unsigned int i = 0 ; i < currentLine.length() ; i++){
			// Comma detected, it was the name of the material before:
			if(currentLine[i] == ',' && MATERIAL){
				currentMaterial.assign(currentLine,0,i);
				if(firstTime){
					firstTime = false;
					this->materialID_FromMaterialName[currentMaterial] = (unsigned char)counterMaterial; 
					this->materialName_FromMaterialID[(unsigned char)counterMaterial] = currentMaterial;
					previousMaterial = currentMaterial;
				}
				MATERIAL = false;
				if(currentMaterial.compare(previousMaterial) != 0){
					counterMaterial ++;
					// Fill in the dictionnary of material and corresponding ID:
					this->materialID_FromMaterialName[currentMaterial] = (unsigned char)counterMaterial; 
					this->materialName_FromMaterialID[(unsigned char)counterMaterial] = currentMaterial;
					this->numberOFTempForTheMaterial.push_back(counterTemp);
					//cout << "COUNTER TEMP : " << counterTemp << endl;
					counterTemp     = 0;
					
					//cout << "On change de matériau.\n";
					
				}
				continue;
			}else if(MATERIAL){
				continue;
			}
			if(currentLine[i] == ',' && MATERIAL == false){
				// Get double from string:	
				size_t offset = 0;			
				double a = stod(valueStr,&offset);
				// Push this double in the properties 3D table:
				this->properties(counterMaterial,counterProperties,counterTemp) = a;
				
				//cout << "Added prop. : " << a << "--" << valueStr << "--";
				//cout <<  "at (" << counterMaterial << "," << counterProperties << "," << counterTemp << ")" << endl;
				
				valueStr = string();
				counterProperties++;
			}else{
				valueStr.push_back(currentLine[i]);
				if(i == currentLine.length()-1){
					size_t offset = 0;
					double a = stod(valueStr,&offset);
					this->properties(counterMaterial,counterProperties,counterTemp) = a;
					
					//cout << "Added prop. : " << a << "--" << valueStr << "--";
					//cout <<  "at (" << counterMaterial << "," << counterProperties << "," << counterTemp << ")" << endl;
					
					valueStr = string();
				}
			}
		}
		counterTemp++;
		counterProperties = 0;
		valueStr = string();
		previousMaterial = currentMaterial;
        }
	this->numberOFTempForTheMaterial.push_back(counterTemp-1);
	//cout << "COUNTER TEMP : " << counterTemp << endl;
	/* Closing file */
	file.close();
	/*
	for(auto& x : this->materialID_FromMaterialName)
	{
		cout << x.first << "," << (int)x.second << endl;
	}
	*/
	
}

/*
 * Materials::getProperty
 * 	Inputs:	1) temperature (we will take the nearest)
 *			2) material
 *			3) property
 *			4) if interpolation=true, interpolate property
 */
double Materials::getProperty(double temperature,unsigned char material, unsigned char property,
							 bool interpolation /*= false*/){
	if(interpolation == true){
		cout << "Sorry, not yet implemented!\n";
		abort();
	}else{
		/* WE DON'T INTERPOLATE */

		if(this->numberOfMaterials < material){
			printf("Materials::getProperty::ERROR\n");
			printf("\tAsking for material %d but has only %d materials in the database.",material,
									this->numberOfMaterials);
			printf("\nAborting.\n\n");
			abort();
		}
		/* Next : smaller or equal because one of the column is for temperature !*/
		if(this->numberOfProperties <= property){
			fprintf(stderr,"Materials::getProperty::ERROR\n");
			fprintf(stderr,"Given property(%d) is higher than the number of properties(%d)",
				property,this->numberOfProperties-1);
			// UP: the -1 is to account for the column with temperature, whch is
			// not a property !
			fprintf(stderr,"Aborting. In file %s:%d\n",__FILE__,__LINE__);
			abort();
		}

		// First, we retrieve the number of temperature info for this material:
		unsigned char numberOfTempLines = this->numberOFTempForTheMaterial[material];
		unsigned int temperatureIndex = 0;
		double tempDiff = 1E10;
		double temp;

		//printf("Materials::getProperty: looking for material %d.\n",material);

		for(unsigned char I = 0 ; I < numberOfTempLines ; I ++){
			
			temp = abs(this->properties(material,0,I)-temperature);
			if( temp < tempDiff){
				tempDiff = temp;
				temperatureIndex = I;
			}else{
				temperatureIndex = I;
				break;
			}
		}
		
		// Find the temperature that is the nearest of the input temperature:
		//for(unsigned int I = 0 ; I < 
		// Accès au tableau: this->properties[material][sigma, mu, autre][Temperature]
		return this->properties(material,property,temperatureIndex);
	}
	return 0.0;
}

/* Destructor */
Materials::~Materials(void){
	/*if(this->properties != NULL)
		this->freeProperties();*/
	#if DEBUG > 2
	cout << "Materials::destructor::out\n";
	#endif
}

/* Constructor */
Materials::Materials(void){
}

void Materials::printAllProperties(void){
	cout << "---------------------PRINT PROPERTIES-----------------------\n";
	cout << "nbrMat = " << this->numberOfMaterials << " - maxNbrTemp = ";
	cout << this->maxNumberOfTemp << " - nbrProp = " << this->numberOfProperties;
	cout << endl << "***---" << endl;
	// Going through the materials:
	for(size_t I = 0 ; I < this->numberOfMaterials ; I ++){
		// Going through each temperature specification:
		for(size_t K = 0 ; K < this->maxNumberOfTemp ; K ++){
			// Going through each property:
			for(size_t J = 0 ; J < this->numberOfProperties ; J ++){
				printf("%f", this->properties(I,J,K));
				cout <<  "(" << I << "," << J << "," << K << ") |";
			}
		cout << endl;
		}
		cout << "----" << endl;
	}
	cout << "---------------------------------------------\n";
}

/*void Materials::freeProperties(void){
	for(unsigned int I = 0 ; I < this->numberOfMaterials ; I ++){
		for(unsigned int J = 0 ; J < this->numberOfProperties ; J ++){
			delete[] this->properties[I][J];
		}
		delete[] this->properties[I];
	}
	delete[] this->properties;
	if(this->properties != NULL)
		this->properties = NULL;
}*/

// Get Dictionnary with the materials and the chosen unsigned char assigned to it:
map<string,unsigned char> Materials::get_dictionnary_MaterialToID(void){
	return this->materialID_FromMaterialName;	
}

map<unsigned char,string> Materials::get_dictionnary_IDToMaterial(void){
	return this->materialName_FromMaterialID;
}

void Materials::printNumberOfTempLinePerMat(void){
	// Print the field numberOFTempForTheMaterial
	for(unsigned int I = 0 ; I < this->numberOFTempForTheMaterial.size() ; I++){
		cout << "For material " + this->materialName_FromMaterialID[I];
		cout << ", there are " << this->numberOFTempForTheMaterial[I];
		cout << " number of temperature infos." << endl;
	}
}
