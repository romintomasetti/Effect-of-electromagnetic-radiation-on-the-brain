#include <vector>
#include <string>
#include <iostream>

using namespace std;

class Materials{
	private:
		// Contains all the properties of all materials:
		vector<vector<vector<double>>> properties;	
	public:
		void   getPropertiesFromFile(string);
		double getProperty(double,unsigned char);
};
