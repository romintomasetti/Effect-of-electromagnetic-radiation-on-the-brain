#include "Materials.h"
#include <fstream>
#include <cstring>

void Materials::getPropertiesFromFile(string filename){
	/* Opening the file */
	ifstream file;
	file.open(filename,fstream::in);
	if(file.fail())
	{
		printf("There was an error while opening the file containing the properties (line %d).\n",__LINE__);
		abort();
	}
	/* Acquiring data */
	string currentLine;
	string currentMaterial, previousMaterial;
	getline(file,currentLine);

	/* Déterminer le nombre de lignes par matériau et prendre le plus grand.*/
	unsigned int maxNumberOfTemp = 0, counter = 0, numberOfMaterials = 0;
	while(!file.eof())
	{
		getline(file,currentLine);
		for(unsigned int i = 0 ; i < currentLine.length() ; i++){
				if(currentLine[i] == ','){
					currentMaterial.assign(currentLine,0,i);
					break;
				}
		}
		if(currentMaterial.compare(previousMaterial) != 0){
			cout << "On change de matériau.\n";
			numberOfMaterials++;
			if(maxNumberOfTemp < counter)
				maxNumberOfTemp = counter;
			counter = 0;
		}else{
			counter++;
		}
		previousMaterial = currentMaterial;
	}
	previousMaterial = string();
	cout << "Max number of temp line is " << maxNumberOfTemp << endl;
	cout << "Number of materials is " << numberOfMaterials << endl;
	file.clear();
	file.seekg(0);
	getline(file,currentLine);
	while(!file.eof())
        {
		// Create a vector for the properties at a given T:
		vector<double> propForThisTemp;
		string valueStr;
        	getline(file,currentLine);
		cout << currentLine << endl;
		bool MATERIAL = true;
		for(unsigned int i = 0 ; i < currentLine.length() ; i++){
			if(currentLine[i] == ',' && MATERIAL){
				currentMaterial.assign(currentLine,0,i);
				MATERIAL = false;
				continue;
			}else if(MATERIAL){
				continue;
			}
			if(currentLine[i] == ',' && MATERIAL == false){
				size_t offset = 0;
				double a = stod(valueStr,&offset);
				propForThisTemp.insert(propForThisTemp.end(),a);
				for (int i = 0; i < (int)propForThisTemp.size(); i++)
        				cout << propForThisTemp.at(i) << ' ';
				cout << endl;
				
				valueStr = string();
			}else{
				valueStr.push_back(currentLine[i]);
			}
		}
		if(currentMaterial.compare(previousMaterial) != 0){
			cout << "On change de matériau.\n";
		}
        	cout<< currentLine << " (prev.)" << previousMaterial << "\n"; 
		previousMaterial = currentMaterial;
        }
	/* Closing file */
	file.close();
	// Lire le fichier et extraire les paramètres.
}


double Materials::getProperty(double,unsigned char){
	// Faire interpolation.
	// Dire que sigma = 0 , mu = 1
	// Accès au tableau: this->properties[material][sigma, mu, autre][Temperature]
	return 0.0;
}
