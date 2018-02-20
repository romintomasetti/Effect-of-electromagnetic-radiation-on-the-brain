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
		//
		size_t nbrNodeCorner1_X;
		size_t nbrNodeCorner1_Y;
		size_t nbrNodeCorner1_Z;
		size_t nodesInsideAlong_X;
		size_t nodesInsideAlong_Y;
		size_t nodesInsideAlong_Z;
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
		// From deltaX, deltaY and deltaX, from centerX, centerY, centerZ,
		// determine the nodes inside the antenna:
		void computeNodesInsideSource(const double,const double, const double,
									 const double, const double, const double);
		// Check that a node is inside the source:
		bool isInsideSource(const size_t, const size_t, const size_t);
};

#endif