#include <vector>
#include <string>
#include <iostream>

using namespace std;



class Materials{
	private:
		// Contains all the properties of all materials:
		double ***properties = NULL;
		// Number of properties:
		unsigned int numberOfProperties = 0;
		// Number of materials:
		unsigned int numberOfMaterials  = 0;
		// Maximum number of temperature specifications:
		unsigned int maxNumberOfTemp    = 0;
		
		// Free the properties array (called in the destructor):
		void   freeProperties(void);
	public:
		// Get all the properties specified in a file, and put them in a 3D array:
		void   getPropertiesFromFile(string);
		// Get a property for a given material at a given temperature:
		double getProperty(double,unsigned char);
		// Print all the properties:
		void   printAllProperties(void);
		// Constructor:
		Materials();
		// Destructor:
		~Materials();
};
