/* This class defines the behavior and properties of an electromagnetic source */
#ifndef ELECTROMAGNETICSOURCE_H
#define ELECTROMAGNETICSOURCE_H

#include <vector>
#include "SetOnceVariable_Template.h"
#include "GridCreator.h"

using namespace std;

class GridCreator;

class ElectromagneticSource{
	private:
		SetOnceVariable_Template<unsigned int> number_of_sources;
		

		std::vector<bool> lengthsAlreadySet;
		std::vector<bool> centersAlreadySet;
		std::vector<bool> nodesInsideAlreadySet;

		//
		std::vector<size_t> nbrNodeCorner1_X;
		std::vector<size_t> nbrNodeCorner1_Y;
		std::vector<size_t> nbrNodeCorner1_Z;
		std::vector<size_t> nodesInsideAlong_X;
		std::vector<size_t> nodesInsideAlong_Y;
		std::vector<size_t> nodesInsideAlong_Z;
	public:
		// Constructor:
		ElectromagneticSource(){};
		// Destructor:
		~ElectromagneticSource(){};
		// Set the number of sources:
		// Lengths and centers in each direction of each source:
		std::vector<double> lengthX;
		std::vector<double> lengthY;
		std::vector<double> lengthZ;
		std::vector<double> centerX;
		std::vector<double> centerY;
		std::vector<double> centerZ;
		std::vector<double> frequency;
		void set_number_of_sources(const unsigned int);
		unsigned int get_number_of_sources(void){
			return this->number_of_sources.get();
		}
		// Set length in each direction:
		void setLengths(const double, const double, const double, const unsigned int);
		// Set lengths along one direction:
		void setLengthAlongOneDir(const unsigned int direction,
									vector<double> values);

		// Set center in each direction:
		void setCenter (const double, const double, const double, const unsigned int);
		// Set centers along one direction:
		void setCenterAlongOneDir(const unsigned int direction,
									vector<double> values);
		// Set frequency:
		void setFrequency(const double, const unsigned int);
		void setAllFrequencies(vector<double> freqs);
		
		// Get length in each direction:
		vector<double> getLengths(const unsigned int i){
			vector<double> lengths = {this->lengthX[i],
										this->lengthY[i],
										this->lengthZ[i]};
			return lengths;
		}
		// Get center in each direction:
		vector<double> getCenter(const unsigned int i){
			vector<double> center = {this->centerX[i],
										this->centerY[i],
										this->centerZ[i]};
			return center;
		}
		// Get frequency:
		double getFrequency(const unsigned int i){
			return this->frequency[i];
		}
		// From deltaX, deltaY and deltaX, from centerX, centerY, centerZ,
		// determine the nodes inside the antenna:
		void computeNodesInsideSource(const double,const double, const double,
									 const double, const double, const double,
									 const unsigned int);
		// Check that a node is inside the source:
		bool isInsideSource(const size_t, const size_t, const size_t,
							const unsigned int);
		//Get value source  
		//From mesh, t_current, i,j,k
        void computeSourceValue(GridCreator, double,int ,int,int,unsigned int);
};

#endif