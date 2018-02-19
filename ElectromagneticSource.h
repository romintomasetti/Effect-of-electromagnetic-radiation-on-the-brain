/* This class defines the behavior and properties of an electromagnetic source */
#ifndef ELECTROMAGNETICSOURCE_H
#define ELECTROMAGNETICSOURCE_H

#include <vector>

using namespace std;

class ElectromagneticSource{
	private:
		double lengthX   = 0.0;
		double lengthY   = 0.0;
		double lengthZ   = 0.0;
		double centerX   = 0.0;
		double centerY   = 0.0;
		double centerZ   = 0.0;
		double frequency = 0.0;
	public:
		// Constructor:
		ElectromagneticSource(){};
		// Destructor:
		~ElectromagneticSource(){};
		// Set length in each direction:
		void setLengths(const double, const double, const double);
		// Set center in each direction:
		void setCenter (const double, const double, const double);
		// Set frequency:
		void setFrequency(const double);
		// Get length in each direction:
		vector<double> getLengths(void){
			vector<double> lengths = {this->lengthX,this->lengthY,this->lengthZ};
			return lengths;
		}
		// Get center in each direction:
		vector<double> getCenter(void){
			vector<double> center = {this->centerX,this->centerY,this->centerZ};
			return center;
		}
		// Get frequency:
		double getFrequency(void){
			return this->frequency;
		}
};

#endif