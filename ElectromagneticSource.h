/* This class defines the behavior and properties of an electromagnetic source */
#ifndef ELECTROMAGNETICSOURCE_H
#define ELECTROMAGNETICSOURCE_H

#include <vector>
#include <algorithm>
#include "SetOnceVariable_Template.h"

#include "header_with_all_defines.hpp"

#include <climits>

#include "UTILS/vector_utilities.hpp"

//#include "GridCreator.h"

using namespace std;

class GridCreator;

class ElectromagneticSource{
	private:
		

		//

	public:
		std::vector<bool> there_is_at_least_one_element_non_zero_in_source;
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
		inline std::string is_inside_source_Romin(
					const size_t I_gl, 
					const size_t J_gl, 
					const size_t K_gl,
					const std::vector<double> &deltas_Electro,
					const std::string &type ,//= "Not_given"
					const std::string &source_type, // = "NOT_GIVEN"
					const unsigned char ID_Source, // = UCHAR_MAX,
					const std::vector<double> &origin_whole_grid)
		{
			/// Verify arguments:
			if(ID_Source == UCHAR_MAX){
				fprintf(stderr,"%s::ID_Source is equal to UCHAR_MAX. Aborting\n",__FUNCTION__);
				fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
				abort();
			}
			if(type == "Not_given"){
				fprintf(stderr,"%s:: you didn't specify the type of the node. Aborting.\n",__FUNCTION__);
				fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
				abort();
			}
			std::vector<std::string> avail_source_types 
				= {
						"DIPOLE",
						"SIMPLE",
						"FACE_EX",
						"FACE_EY",
						"FACE_EZ",
						"FACE_Minus_EX",
						"FACE_Minus_EY",
						"FACE_Minus_EZ"
				};
			bool source_types_is_ok = false;
			if(std::find(avail_source_types.begin(), avail_source_types.end(), source_type) 
					!= avail_source_types.end()) {
				source_types_is_ok = true;
			} else {
				source_types_is_ok = false;
			}
			if(source_type == "NOT_GIVEN" || !source_types_is_ok){
				DISPLAY_ERROR_ABORT(
					"Wrong source type. Has %s but accepted types are:\n%s.",
					source_type.c_str(),
					vector_string_to_one_string(avail_source_types,"\t>","\n").c_str()
				);
			}
			if(source_type == "SIMPLE"){

				/// If not electric field type, return 0:
				if( type != "Ex" && type != "Ey" && type != "Ez"){
					return "false";
				}

				/// Compute the coordinates of the node w.r.t. the origin of the whole grid:
				double X_coord = origin_whole_grid[0] + I_gl * deltas_Electro[0];
				double Y_coord = origin_whole_grid[1] + J_gl * deltas_Electro[1];
				double Z_coord = origin_whole_grid[2] + K_gl * deltas_Electro[2];

				/// Determine if it is inside the source:
				double EPS = deltas_Electro[0]*1E-5;
				if(    X_coord >= (this->centerX[ID_Source] - this->lengthX[ID_Source]/2.)-EPS
					&& X_coord <= (this->centerX[ID_Source] + this->lengthX[ID_Source]/2.)+EPS
					&& Y_coord >= (this->centerY[ID_Source] - this->lengthY[ID_Source]/2.)-EPS
					&& Y_coord <= (this->centerY[ID_Source] + this->lengthY[ID_Source]/2.)+EPS
					&& Z_coord >= (this->centerZ[ID_Source] - this->lengthZ[ID_Source]/2.)-EPS
					&& Z_coord <= (this->centerZ[ID_Source] + this->lengthZ[ID_Source]/2.)+EPS)
				{
					/// The node is inside the source.
					bool isOnFace_e_x = false;
					bool isOnFace_e_y = false;

					if(X_coord + deltas_Electro[0] >= (this->centerX[ID_Source] + this->lengthX[ID_Source]/2)+EPS){
						isOnFace_e_x = true;
					}

					if(Y_coord + deltas_Electro[1] >= (this->centerY[ID_Source] + this->lengthY[ID_Source]/2)+EPS){
						isOnFace_e_y = true;
					}

					if(isOnFace_e_x == true && isOnFace_e_y == true){
						// Impose none of Ex and Ey, return false:
						if( type == "Ex" || type == "Ey" ){
							return "false";
						}
					}else if(isOnFace_e_x == true){
						if(type == "Ex"){
							return "false";
						}else if(type == "Ey"){
							return "true";
						}
					}else if(isOnFace_e_y == true){
						if(type == "Ex"){
							return "true";
						}else if(type == "Ey"){
							return "false";
						}
					}
					this->there_is_at_least_one_element_non_zero_in_source[ID_Source] = true;
					return "true";
				}
			}else if(source_type == "DIPOLE"){

				/// If not electric field type, return 0:
				if( type != "Ex" && type != "Ey" && type != "Ez"){
					return "false";
				}


				/// Retrieve the frequency of the source:
				double speedOfLight = 3E8;
				double freq   = this->frequency[ID_Source];
				double lambda = speedOfLight / freq;

				double length_X = lambda/4;
				double length_Y = lambda/4;
				double length_Z = 2 * lambda/4 + deltas_Electro[2];

				/// Compute the coordinates of the node w.r.t. the origin of the whole grid:
				double X_coord = origin_whole_grid[0] + I_gl * deltas_Electro[0];
				double Y_coord = origin_whole_grid[1] + J_gl * deltas_Electro[1];
				double Z_coord = origin_whole_grid[2] + K_gl * deltas_Electro[2];

				double EPS = deltas_Electro[0]*1E-5;

				if(    X_coord >= (this->centerX[ID_Source] - length_X/2.)-EPS
					&& X_coord <= (this->centerX[ID_Source] + length_X/2.)+EPS
					&& Y_coord >= (this->centerY[ID_Source] - length_Y/2.)-EPS
					&& Y_coord <= (this->centerY[ID_Source] + length_Y/2.)+EPS
					&& Z_coord >= (this->centerZ[ID_Source] - length_Z/2.)-EPS
					&& Z_coord <= (this->centerZ[ID_Source] + length_Z/2.)+EPS)
				{
					/// Check if Ez is inside the airgap:
					if( abs(Z_coord - this->centerZ[ID_Source] ) <= 2*EPS ){
						if( type == "Ez"){
							this->there_is_at_least_one_element_non_zero_in_source[ID_Source] = true;
							return "true";
						}else if(type == "Ex" || type == "Ey"){
							return "false";
						}
					}
					/// The node is inside the source.
					bool isOnFace_e_PlusX  = false;
					bool isOnFace_e_MinusX = false;
					bool isOnFace_e_PlusY  = false;
					bool isOnFace_e_MinusY = false;
					bool isOnFace_e_PlusZ  = false;
					bool isOnFace_e_MinusZ = false;

					if(X_coord + deltas_Electro[0] >= (this->centerX[ID_Source] + length_X/2.)+EPS){
						isOnFace_e_PlusX  = true;
					}
					if(X_coord - 0*deltas_Electro[0] <= (this->centerX[ID_Source] - length_X/2.)-EPS){
						isOnFace_e_MinusX = true;
					}

					if(Y_coord + deltas_Electro[1] >= (this->centerY[ID_Source] + length_Y/2.)+EPS){
						isOnFace_e_PlusY = true;
					}
					if(Y_coord - 0*deltas_Electro[1] <= (this->centerY[ID_Source] - length_Y/2.)-EPS){
						isOnFace_e_MinusY = true;
					}

					if(Z_coord + deltas_Electro[2] >= (this->centerZ[ID_Source] + length_Z/2.)+EPS){
						isOnFace_e_PlusZ = true;
					}
					if(Z_coord + deltas_Electro[2] <= (this->centerZ[ID_Source] - length_Z/2.)-EPS){
						isOnFace_e_MinusZ = true;
					}

					if(    (isOnFace_e_PlusX  == true && isOnFace_e_PlusY  == true)
						|| (isOnFace_e_PlusX  == true && isOnFace_e_MinusY == true)
						|| (isOnFace_e_MinusX == true && isOnFace_e_MinusY == true)
						|| (isOnFace_e_MinusX == true && isOnFace_e_PlusY  == true)){
						// Edge. Impose both Ex, Ey and Ez to zero.
						if( type == "Ex" || type == "Ey"){
							return "false";
						}
					}else if(isOnFace_e_PlusX == true || isOnFace_e_MinusX == true){
						// Face with normal (+x) or (-x). Impose Ey and Ez to zero.
						if( type == "Ey" || type == "Ez" ){
							return "0";
						}
						if( type == "Ex"){
							return "false";
						}
						
					}else if(isOnFace_e_PlusY == true || isOnFace_e_MinusY == true){
						// Face with normal (+y) or (-y). Impose both Ex and Ez to zero.
						if(type == "Ex" || type == "Ez"){
							return "0";
						}else if(type == "Ey"){
							return "false";
						}
					}else if(isOnFace_e_MinusZ == true || isOnFace_e_PlusZ == true){
						// Face with normal (+z) or (-z). Impose Ex=Ey=0.
						if( type == "Ex" || type == "Ey" || type == "Ez" ){
							return "0";
						}else{
							return "false";
						}
					}

					// In the bulk, impose Ex, Ey and Ez to zero:
					return "0";
				}
				
			}else if(
				   source_type == "FACE_EX"
				|| source_type == "FACE_EY"
				|| source_type == "FACE_EZ"
				|| source_type == "FACE_Minus_EX"
				|| source_type == "FACE_Minus_EY"
				|| source_type == "FACE_Minus_EZ"){
				/**
				 * @brief Put a source on the face with normal.
				 * 	This kind of source has been hardcoded !
				 */
				this->there_is_at_least_one_element_non_zero_in_source[ID_Source] = true;
				return "false";
			}else{
				DISPLAY_ERROR_ABORT(
					"The source type %s doesn't match any source implementation."
					"Available source types are:\n%s",
					source_type.c_str(),
					vector_string_to_one_string(avail_source_types,"\t>","\n").c_str()
				);
			}
			
			/// By default, return false.
			return "false";
		}
};

#endif
