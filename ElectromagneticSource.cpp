#include "ElectromagneticSource.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>


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
bool ElectromagneticSource::is_inside_source_Romin(
			const size_t I_gl, 
			const size_t J_gl, 
			const size_t K_gl,
			const std::vector<double> &deltas_Electro,
			const std::string &type /*= "Not_given"*/,
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
		fprintf(stderr,"%s:: you didn't specify the type of the node. Abortin.\n",__FUNCTION__);
		fprintf(stderr,"In %s:%d\n",__FILE__,__LINE__);
		abort();
	}
	/// As a function of the type of node, the spatial shift is different.

	double shift_X = 0.0;
	double shift_Y = 0.0;
	double shift_Z = 0.0;

	double X_coord = -1.0;
	double Y_coord = -1.0;
	double Z_coord = -1.0;

	if(type == "Ex"){
		// The node is of part of the electric field's X component.
		// Shift for this type of node is deltaX along X.
		//shift_X = deltas_Electro[0];
		//I_shift++;


	}else if (type == "Ey"){
		// The node is of part of the electric field's Y component.
		// Shift for this type of node is deltaY along Y.
		//shift_Y = deltas_Electro[1];
		//J_shift++;

	}else if (type == "Ez"){
		// The node is of part of the electric field's Z component.
		// Shift for this type of node is deltaZ along Z.
		//shift_Z = deltas_Electro[2];

	}else if (type == "Hx"){
		// The node is of part of the magnetic field's X component.
		// Shift for this type of node is deltaY along Y and deltaZ along Z.
		//shift_Y = deltas_Electro[1];
		//shift_Z = deltas_Electro[2];

	}else if (type == "Hy"){
		// The node is of part of the magnetic field's Y component.
		// Shift for this type of node is deltaX along X and deltaZ along Z.
		//shift_X = deltas_Electro[0];
		//shift_Z = deltas_Electro[2];

	}else if (type == "Hz"){
		// The node is of part of the magnetic field's Z component.
		// Shift for this type of node is deltaX along X and deltaY along Y.
		//shift_X = deltas_Electro[0];
		//shift_Y = deltas_Electro[1];

	}else{
		fprintf(stderr,"In %s :: invalid type (has %s). Aborting.\n",
			__FUNCTION__,type.c_str());
		fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
	}

	/// Compute the coordinates of the node w.r.t. the origin of the whole grid:
	X_coord = origin_whole_grid[0] + I_gl * deltas_Electro[0] + shift_X;
	Y_coord = origin_whole_grid[1] + J_gl * deltas_Electro[1] + shift_Y;
	Z_coord = origin_whole_grid[2] + K_gl * deltas_Electro[2] + shift_Z;

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
				printf("Sur l'arête : (%zu,%zu,%zu)\n",I_gl,J_gl,K_gl);
				return false;
			}
		}else if(isOnFace_e_x == true){
			if(type == "Ex"){
				return false;
			}else if(type == "Ey"){
				return true;
			}
		}else if(isOnFace_e_y == true){
			if(type == "Ex"){
				return true;
			}else if(type == "Ey"){
				return false;
			}
		}

		return true;
	}
	
	/// By default, return false.
	return false;
}