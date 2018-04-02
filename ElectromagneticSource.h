/* This class defines the behavior and properties of an electromagnetic source */
#ifndef ELECTROMAGNETICSOURCE_H
#define ELECTROMAGNETICSOURCE_H

#include <vector>
#include "SetOnceVariable_Template.h"

#include <climits>

//#include "GridCreator.h"

using namespace std;

class GridCreator;

class ElectromagneticSource{
	private:
		

		//

	public:
		// Constructor:
		ElectromagneticSource(){};
		// Destructor:
		~ElectromagneticSource(){};
		
		SetOnceVariable_Template<unsigned int> number_of_sources;
		

		std::vector<bool> lengthsAlreadySet;
		std::vector<bool> centersAlreadySet;
		std::vector<bool> nodesInsideAlreadySet;
		
		
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
		
		void set_airGaps(const std::vector<double> airGaps);

		// Get length in each direction:
		vector<double> getLengths(const unsigned int i){
			vector<double> lengths = {this->lengthX[i],
										this->lengthY[i],
										this->lengthZ[i]};
			return lengths;
		}
		// Get center in each direction:
		void getCenter(const unsigned int i, double *vec){
			/*vector<double> center = {this->centerX[i],
										this->centerY[i],
										this->centerZ[i]};
			return center;*/
			vec[0] = this->centerX[i];
			vec[1] = this->centerY[i];
			vec[2] = this->centerZ[i];
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
		// Check that a node is inside the source, you give the source:
		bool isInsideSource(const size_t, const size_t, const size_t,const unsigned int);
		// Check that a node is inside a source, don't specify which one:
		bool isInsideSource(const size_t, const size_t, const size_t);
		//Get value source  
		//From mesh, t_current, i,j,k
        void computeSourceValue(GridCreator&,const size_t,const size_t,const size_t,double,char);

		int DetermineInWhichSourceWeAre(const size_t,const size_t,const size_t);

		bool isInsideAirGap(const size_t i_global, 
							const size_t j_global, 
							const size_t k_global,
							const unsigned int i);

		bool isInsideAirGap(const size_t i_global, 
											const size_t j_global, 
											const size_t k_global);

		/**
		 * @brief Returns true if the node is located inside or on the source.
		 * 
		 * Argments are: 1) Global indices of the node.
		 * 				 2) Spatial step along each direction.
		 *               3) Type of the node. Accepted types are:
		 * 						Ex, Ey, Ez, Hx, Hy, Hz
		 * 					Default is "Not_given", which leads to aborting.
		 *               4) Source ID. Default is UCHAR_MAX (~255), which leads to aborting.
		 *               5) Origin of the whole grid (EM), in coordinates.
		 *                  By default, the origin is (0,0,0).
		 * 
		 */
		std::string is_inside_source_Romin(
			const size_t I_gl, 
			const size_t J_gl, 
			const size_t K_gl,
			const std::vector<double> &deltas_Electro,
			const std::string &type = "Not_given",
			const std::string &source_type = "NOT_GIVEN",
			const unsigned char ID_Source = UCHAR_MAX,
			const std::vector<double> &origin_whole_grid = {0.0,0.0,0.0});
};

#endif
