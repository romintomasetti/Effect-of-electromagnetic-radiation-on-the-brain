#include "ElectromagneticSource.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>

// Set the number of sources:
void ElectromagneticSource::set_number_of_sources(const unsigned int nbrSources){
	if(!this->number_of_sources.get_alreadySet()){
		this->number_of_sources = nbrSources;

		this->lengthsAlreadySet.reserve(this->number_of_sources.get());
		this->centersAlreadySet.reserve(this->number_of_sources.get());
		this->nodesInsideAlreadySet.reserve(this->number_of_sources.get());

		this->lengthX.reserve(this->number_of_sources.get());
		this->lengthY.reserve(this->number_of_sources.get());
		this->lengthZ.reserve(this->number_of_sources.get());
		this->centerX.reserve(this->number_of_sources.get());
		this->centerY.reserve(this->number_of_sources.get());
		this->centerZ.reserve(this->number_of_sources.get());

		this->nbrNodeCorner1_X.reserve(this->number_of_sources.get());
		this->nbrNodeCorner1_Y.reserve(this->number_of_sources.get());
		this->nbrNodeCorner1_Z.reserve(this->number_of_sources.get());

		this->nodesInsideAlong_X.reserve(this->number_of_sources.get());
		this->nodesInsideAlong_Y.reserve(this->number_of_sources.get());
		this->nodesInsideAlong_Z.reserve(this->number_of_sources.get());
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
		}
	}else if(direction == 1){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->lengthY[I] = values[I];
		}
	}else if(direction == 2){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->lengthZ[I] = values[I];
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
		}
	}else if(direction == 1){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->centerY[I] = values[I];
		}
	}else if(direction == 2){
		for(unsigned int I = 0 ; I < values.size() ; I ++){
			this->centerZ[I] = values[I];
		}
	}else{
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

void ElectromagneticSource::computeNodesInsideSource(const double L_dom_X,
													 const double L_dom_Y,
													 const double L_dom_Z,
													 const double deltaX,
													 const double deltaY,
													 const double deltaZ,
													 const unsigned int i){
	// Number of the node of corner 1
	if(this->number_of_sources.get_alreadySet() == false ||
			this->centersAlreadySet[i] != true ||
			this->lengthsAlreadySet[i] != true){
		printf("ElectromagneticSource::computeNodesInsideSource::ERROR\n");
		printf("The number of sources hasn't been set, aborting.\n");
		std::abort();
	}
	this->nbrNodeCorner1_X[i] = (size_t) (this->centerX[i] - this->lengthX[i] / 2) / deltaX;
	this->nbrNodeCorner1_Y[i] = (size_t) (this->centerY[i] - this->lengthY[i] / 2) / deltaY;
	this->nbrNodeCorner1_Z[i] = (size_t) (this->centerZ[i] - this->lengthZ[i] / 2) / deltaZ;
	
	this->nodesInsideAlong_X[i] = (size_t) this->lengthX[i] / deltaX;
	this->nodesInsideAlong_Y[i] = (size_t) this->lengthY[i] / deltaY;
	this->nodesInsideAlong_Z[i] = (size_t) this->lengthZ[i] / deltaZ;

	this->nodesInsideAlreadySet[i] = true;
}

bool ElectromagneticSource::isInsideSource(const size_t x, 
																					 const size_t y, 
																					 const size_t z,
																					 const unsigned int i){
	// i represents the desired source.
	if(this->nodesInsideAlreadySet[i] != true){
		printf("ElectromagneticSource::isInsideSource::ERROR.\n");
		printf("Please call ElectromagneticSource::computeNodesInsideSource before.\n");
		printf("Aborting.\n\n");
		std::abort();
	}
	if( ( x >= this->nbrNodeCorner1_X[i] &&
		  x <= (this->nbrNodeCorner1_X[i] + this->nodesInsideAlong_X[i]))
		&&
		( y >= this->nbrNodeCorner1_Y[i] &&
		  y <= (this->nbrNodeCorner1_Y[i] + this->nodesInsideAlong_Y[i]))
		&&
		( z >= this->nbrNodeCorner1_Z[i] &&
		  z <= (this->nbrNodeCorner1_Z[i] + this->nodesInsideAlong_Z[i]))){
		return true;
	}
	// By default, return false:
	return false;
}