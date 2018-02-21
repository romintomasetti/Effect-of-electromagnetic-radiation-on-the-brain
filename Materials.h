#ifndef MATERIALS_H
#define MATERIALS_H

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <limits.h>
#include <cmath>

#include "Array_3D_Template.h"

using namespace std;



class Materials{
	private:
		// Contains all the properties of all materials:
		Array_3D_Template<double> properties;
		// Number of properties:
		unsigned int numberOfProperties = 0;
		
		// Maximum number of temperature specifications:
		unsigned int maxNumberOfTemp    = 0;
		// Dictionnary with the materials and the chosen unsigned char assigned to it:
		map<string,unsigned char> materialID_FromMaterialName;
		map<unsigned char,string> materialName_FromMaterialID;
		
		vector<unsigned int> numberOFTempForTheMaterial;
	
		// Free the properties array (called in the destructor):
		//void   freeProperties(void);
	public:
		// Get all the properties specified in a file, and put them in a 3D array:
		void   getPropertiesFromFile(string);
		// Get a property for a given material at a given temperature:
		double getProperty(double, unsigned char, unsigned char,bool interpolation = false);
		// Print all the properties:
		void   printAllProperties(void);
		// Print the number of temperature lines per material:
		void   printNumberOfTempLinePerMat(void);
		// Number of materials (pour AlgoElectro.cpp)
		unsigned int numberOfMaterials  = 0;
		// Constructor:
		Materials();
		// Destructor:
		~Materials();
		// Get Dictionnary with the materials and the chosen unsigned char assigned to it:
		map<string,unsigned char> get_dictionnary_MaterialToID(void);
		map<unsigned char,string> get_dictionnary_IDToMaterial(void);
};

#endif
