#include "ElectromagneticSource.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

#include <cmath>

#ifndef M_PI
  #define M_PI 3.14
#endif 
// Set the number of sources:
void ElectromagneticSource::set_number_of_sources(const unsigned int nbrSources){
	if(!this->number_of_sources.get_alreadySet()){
		
		this->number_of_sources = nbrSources;

		this->lengthsAlreadySet.reserve(this->number_of_sources.get());
		this->centersAlreadySet.reserve(this->number_of_sources.get());

		this->lengthX.reserve(this->number_of_sources.get());
		this->lengthY.reserve(this->number_of_sources.get());
		this->lengthZ.reserve(this->number_of_sources.get());
		this->centerX.reserve(this->number_of_sources.get());
		this->centerY.reserve(this->number_of_sources.get());
		this->centerZ.reserve(this->number_of_sources.get());
	}else{
		printf("ElectromagneticSource::set_number_of_sources::ERROR\n");
		printf("\tThe property of the field 'number_of_sources' was already set.");
		printf("\n\tAborting (at file %s at line %d)\n\n",__FILE__,__LINE__);
		abort();
	}
}

void ElectromagneticSource::setLengthAlongOneDir(
									const unsigned int direction,
									vector<double> values){
	if(values.size() != this->number_of_sources.get()){
		printf("You should give as many lengths in the specified direction");
		printf(" as there are sources. You gave %ld",values.size());
		printf(" values for %d sources. Aborting.\n",
										this->number_of_sources.get());
		std::abort();
	}

	if(direction == 0){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->lengthX[I] = values[I];
			this->lengthsAlreadySet[I] = true;
		}
	}else if(direction == 1){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->lengthY[I] = values[I];
			this->lengthsAlreadySet[I] = true;
		}
	}else if(direction == 2){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->lengthZ[I] = values[I];
			this->lengthsAlreadySet[I] = true;
		}
	}else{
		printf("Direction should be between 0 and 2. Aborting.\n");
		std::abort();
	}
}

void ElectromagneticSource::setAllFrequencies(vector<double> freqs){
	if(freqs.size() != this->number_of_sources.get()){
		printf("You should give as many frequencies");
		printf(" as there are sources. You gave %ld",freqs.size());
		printf(" values for %d sources. Aborting.\n",
										this->number_of_sources.get());
		std::abort();
	}
	this->frequency = freqs;
}

// Set centers along one direction:
void ElectromagneticSource::setCenterAlongOneDir(
								const unsigned int direction,
								vector<double> values){
	if(values.size() != this->number_of_sources.get()){
		printf("You should give as many centers in the specified direction");
		printf(" as there are sources. You gave %ld",values.size());
		printf(" values for %d sources. Aborting.\n",
										this->number_of_sources.get());
		std::abort();
	}

	if(direction == 0){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->centerX[I] = values[I];
			this->centersAlreadySet[I] = true;
		}
	}else if(direction == 1){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->centerY[I] = values[I];
			this->centersAlreadySet[I] = true;
		}
	}else if(direction == 2){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->centerZ[I] = values[I];
			this->centersAlreadySet[I] = true;
		}
	}else{
		printf("ElectromagneticSource::setCenterAlongOneDir\n");
		printf("Direction should be between 0 and 2. Aborting.\n");
		std::abort();
	}
}

void ElectromagneticSource::setLengths(const double L_X,
										const double L_Y, 
										const double L_Z,
										const unsigned int i){
	// i represents the desired source, from 0 to number_of_source.
	if(this->number_of_sources.get_alreadySet() && 
			this->lengthsAlreadySet[i] != true){
		this->lengthX[i] = L_X;
		this->lengthY[i] = L_Y;
		this->lengthZ[i] = L_Z;
		this->lengthsAlreadySet[i] = true;
	}
}

void ElectromagneticSource::setCenter (const double C_X,
																			 const double C_Y, 
																			 const double C_Z,
																			 const unsigned int i){
	// i represents the desired source, from 0 to number_of_source.
	if(this->number_of_sources.get_alreadySet() &&
			this->centersAlreadySet[i] != true){
		this->centerX[i] = C_X;
		this->centerY[i] = C_Y;
		this->centerZ[i] = C_Z;
		this->centersAlreadySet[i] = true;
	}
}

/////////////////////////////
std::string ElectromagneticSource::is_inside_source_Romin(
			const size_t I_gl, 
			const size_t J_gl, 
			const size_t K_gl,
			const std::vector<double> &deltas_Electro,
			const std::string &type /*= "Not_given"*/,
			const std::string &source_type /* = "NOT_GIVEN" */,
			const unsigned char ID_Source/* = UCHAR_MAX*/,
			const std::vector<double> &origin_whole_grid/* = {0.0,0.0,0.0}*/)
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
	std::vector<std::string> avail_source_types = {"DIPOLE","SIMPLE"};
	bool source_types_is_ok = false;
	if(std::find(avail_source_types.begin(), avail_source_types.end(), source_type) 
			!= avail_source_types.end()) {
		source_types_is_ok = true;
	} else {
		source_types_is_ok = false;
	}
	if(source_type == "NOT_GIVEN" || !source_types_is_ok){
		fprintf(stderr,"In %s :: ERROR :: wrong source type ! Aborting.\n",
						__FUNCTION__);
		fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
		#ifdef MPI_COMM_WORLD
			MPI_Abort(MPI_COMM_WORLD,-1);
		#else
			abort();
		#endif
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
					//printf("Sur l'arête : (%zu,%zu,%zu)\n",I_gl,J_gl,K_gl);
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

		double shift_X = 0.0;
		double shift_Y = 0.0;

		//if(type == "Ex")
			//shift_X = -deltas_Electro[0];

		//if(type == "Ey")
			//shift_Y = -deltas_Electro[1];

		/// Compute the coordinates of the node w.r.t. the origin of the whole grid:
		double X_coord = origin_whole_grid[0] + I_gl * deltas_Electro[0] + shift_X;
		double Y_coord = origin_whole_grid[1] + J_gl * deltas_Electro[1] + shift_Y;
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
			if( abs(Z_coord - this->centerZ[ID_Source] ) < EPS ){
				if( type == "Ez"){
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
			if(X_coord - deltas_Electro[0] <= (this->centerX[ID_Source] - length_X/2.)-EPS){
				isOnFace_e_MinusX = true;
			}

			if(Y_coord + deltas_Electro[1] >= (this->centerY[ID_Source] + length_Y/2.)+EPS){
				isOnFace_e_PlusY = true;
			}
			if(Y_coord - deltas_Electro[1] <= (this->centerY[ID_Source] - length_Y/2.)+EPS){
				isOnFace_e_MinusY = true;
			}

			if(    (isOnFace_e_PlusX  == true && isOnFace_e_PlusY  == true)
				|| (isOnFace_e_PlusX  == true && isOnFace_e_MinusY == true)
				|| (isOnFace_e_MinusX == true && isOnFace_e_MinusY == true)
				|| (isOnFace_e_MinusX == true && isOnFace_e_PlusY  == true)){
				// Edge. Impose both Ex, Ey and Ez to zero.
				if( type == "Ex" || type == "Ey"){
					return "false";
				}
			}else if(isOnFace_e_PlusX == true || isOnFace_e_MinusX){
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
				// Face with normal (+z) or (-z). Do nothing.
				return "false"; 
			}

			// In the bulk, impose Ex, Ey and Ez to zero:
			return "0";
		}
		
	}
	
	/// By default, return false.
	return "false";
}
