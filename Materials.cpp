#include "Materials.h"
#include <fstream>
#include <cstring>

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

void Materials::getPropertiesFromFile(string filename){
	/* Opening the file */
	ifstream file;
	file.open(filename,fstream::in);
	if(file.fail())
	{
		/* Opening failed - aborting ... */
		printf("There was an error while opening the file containing the properties (line %d).\n",__LINE__);
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
	unsigned int maxNumberOfTemp = 0, counter = 1, numberOfMaterials = -1;
	// Note: numberOfMaterials starts at -1 because at the first comparaison between 
	// 	 current- and previousMaterial, it will always trigger "numberOfMaterials++"
	//	 even if it is not necessary.
	bool numberOfPropertiesIsKnown = false;
	while(!file.eof())
	{
		// Get the current line:
		getline(file,currentLine);
		if(currentLine.empty()){
			// Do nothing, the line is empty !
		}else{
			// Get the position of the first comma:
			for(unsigned int i = 0 ; i < currentLine.length() ; i++){
					if(currentLine[i] == ','){
						// The current material is for example AIR,
						currentMaterial.assign(currentLine,0,i);
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
				// We want to track the material for which we have the highest number of 
				// properties specifications - it will be useful to initialize the table with
				// all the properties.
				counter ++;
				if(maxNumberOfTemp < counter)
					maxNumberOfTemp = counter;
				counter = 1;
			}else{
				counter++;
			}
			previousMaterial = currentMaterial;
		}
	}
	// Set the oject fields:
	this->numberOfMaterials  = numberOfMaterials;
	this->numberOfProperties = numberOfProperties;
	this->maxNumberOfTemp    = maxNumberOfTemp; 	

	// Allocate the this->properties 3d array.
	this->properties = new double**[numberOfMaterials];
	for(unsigned int I = 0 ; I < numberOfMaterials ; I ++){
		this->properties[I] = new double*[numberOfProperties];
		for(unsigned int J = 0 ; J < numberOfProperties ; J ++){
			this->properties[I][J] = new double[maxNumberOfTemp];
			for(unsigned int K = 0 ; K < maxNumberOfTemp ; K ++){
				this->properties[I][J][K] = -1.;
			}
		}
	}
	
	// Empty the previous material variable for later.
	previousMaterial = string();
	#if DEBUG > 2
		cout << "The highest number of specification is " << maxNumberOfTemp << "." << endl;
		cout << "Number of different materials is " << numberOfMaterials << "." << endl;
		cout << "Number of properties is " << this->numberOfProperties << endl;
	#endif
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
					this->materialNameForMaterialID[currentMaterial] = (unsigned char)counterMaterial;
					previousMaterial = currentMaterial;
				}
				MATERIAL = false;
				if(currentMaterial.compare(previousMaterial) != 0){
					counterMaterial ++;
					// Fill in the dictionnary of material and corresponding ID:
					this->materialNameForMaterialID[currentMaterial] = (unsigned char)counterMaterial; 
					counterTemp     = 0;
					#if DEBUG > 2
					cout << "On change de matériau.\n";
					#endif
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
				this->properties[counterMaterial][counterProperties][counterTemp] = a;
				#if DEBUG > 2
				cout << "Added prop. : " << a << "--" << valueStr << "--";
				cout <<  "at (" << counterMaterial << "," << counterProperties << "," << counterTemp << ")" << endl;
				#endif
				valueStr = string();
				counterProperties++;
			}else{
				valueStr.push_back(currentLine[i]);
				if(i == currentLine.length()-1){
					size_t offset = 0;
					double a = stod(valueStr,&offset);
					this->properties[counterMaterial][counterProperties][counterTemp] = a;
					valueStr = string();
				}
			}
		}
		counterTemp++;
		counterProperties = 0;
		valueStr = string();
		previousMaterial = currentMaterial;
        }
	/* Closing file */
	file.close();
	#if DEBUG > 0
	for(auto& x : this->materialNameForMaterialID)
	{
		cout << x.first << "," << (int)x.second << endl;
	}
	#endif
}


double Materials::getProperty(double,unsigned char){
	// Faire interpolation.
	// Dire que sigma = 0 , mu = 1
	// Accès au tableau: this->properties[material][sigma, mu, autre][Temperature]
	return 0.0;
}

/* Destructor */
Materials::~Materials(void){
	if(this->properties != NULL)
		this->freeProperties();
}

/* Constructor */
Materials::Materials(void){
}

void Materials::printAllProperties(void){
	cout << "---------------------PRINT PROPERTIES-----------------------\n";
	// Going through the materials:
	for(unsigned int I = 0 ; I < this->numberOfMaterials ; I ++){
		// Going through each temperature specification:
		for(unsigned int K = 0 ; K < this->maxNumberOfTemp ; K ++){
			// Going through each property:
			for(unsigned int J = 0 ; J < this->numberOfProperties ; J ++){
				cout << this->properties[I][J][K];
				cout <<  "(" << I << "," << J << "," << K << ") |";
			}
		cout << endl;
		}
	}
	cout << "---------------------------------------------\n";
}

void Materials::freeProperties(void){
	for(unsigned int I = 0 ; I < this->numberOfMaterials ; I ++){
		for(unsigned int J = 0 ; J < this->numberOfProperties ; J ++){
			delete[] this->properties[I][J];
		}
		delete[] this->properties[I];
	}
	delete[] this->properties;
	if(this->properties != NULL)
		this->properties = NULL;
}

// Get Dictionnary with the materials and the chosen unsigned char assigned to it:
map<string,unsigned char> Materials::get_dictionnary_MaterialToID(void){
	return this->materialNameForMaterialID;	
}
