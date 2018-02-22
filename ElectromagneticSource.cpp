#include "ElectromagneticSource.h"
#include "GridCreator.h"

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

		this->airgap.reserve(this->number_of_sources.get());

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
											const size_t z){

	/* DETERMINE IN WHICH SOURCE WE ARE */
	
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


void ElectromagneticSource::set_airGaps(const std::vector<double> airGaps){
	if(this->number_of_sources.get_alreadySet()){
		printf("ElectromagneticSource::computeNodesInsideSource::ERROR\n");
		printf("The number of sources hasn't been set, aborting.\n");
		std::abort();
	}

	for(unsigned int I = 0 ; I < this->number_of_sources.get() ; I ++ ){
		this->airgap[I] = airGaps[I];
	}
}

/* --------------------------------------------------------------------------------------------------------------------- */
/* Here, i,j,k are local indices */
void ElectromagneticSource::computeSourceValue(GridCreator &mesh,
				 double tCurrent, int i_global, int j_global,
				 int k_global,char CHAMP)
{
	//double AirGap = 1;


	/* FIND IN WHICH SOURCE WE ARE */

	/* Size d'une antenne du dipole */
	double LengthDipoleX = this->lengthX[ID_Source];//mesh.elec_source.LengthX(ID_Source);
	double LengthDipoleY = this->lengthY[ID_Source];//mesh.elec_source.LengthY(ID_Source);
	double LengthDipoleZ = this->lengthZ[ID_Source];//(mesh.elec_source.LengthZ(ID_Source) - AirGap)/2;

	/* 	Here, GlobalIndices will contain the global indices corresponding to the local indices i,j,k, which are in the source */
	/*int GlobalIndices[3];
	GlobalIndices[0] = mesh.Transformation(i, j, k, mesh.myrank, mesh.numberofprocess);
	GlobalIndices[1] = mesh.Transformation(i, j, k, mesh.myrank, mesh.numberofprocess);
	GlobalIndices[2] = mesh.Transformation(i, j, k, mesh.myrank, mesh.numberofprocess);
	*/
	/* We know that for a dipole antenna E_x, E_y, H_x, H_y and H_z are all equal to 0, whereas E_z is different if we are in the air gap or not */
	if(CHAMP == 'H'){
		mesh.nodeMagn(i,j,k).field[0] = 0.0;
		mesh.nodeMagn(i,j,k).field[1] = 0.0;
		mesh.nodeMagn(i,j,k).field[2] = 0.0;
	}else if(CHAMP == 'E'){
		mesh.nodeElec(i,j,k).field[0] = 0.0;
		mesh.nodeElec(i,j,k).field[1] = 0.0;
		/* Ask Romin if the function returns the indices or the physical coordinates */
		double CenterAntenna[3] = mesh.elec_source.getCenter(ID_Source);

		/* If we are in the antenna */
		if(CenterAntenna[0]-(LengthDipoleX/2)/mesh.deltaX <= GlobalIndices[0]  && GlobalIndices[0] <= CenterAntenna[0]-(LengthDipoleX/2)/mesh.deltaX )
		{
			if(CenterAntenna[1]-(LengthDipoleY/2)/mesh.deltaY <= GlobalIndices[1] && GlobalIndices[1] <= CenterAntenna[1]-(LengthDipoleY/2)/mesh.deltaY)
			{
				if(CenterAntenna[2]-(LengthDipoleZ/2)/mesh.deltaZ <= GlobalIndices[2] && GlobalIndices[2] <= CenterAntenna[2]-(LengthDipoleY/2)/mesh.deltaZ)
					mesh.nodesElec(i,j,k).field[2] = sin(2*M_PI*mesh.elec_source.frequency[ID_Source]*tCurrent);
			}
		}
		else
		{
			mesh.nodesElec(i,j,k).field[2] = 0.0;	
		}
	}else{
		printf("ElectromagneticSource::computeSourceValue:: ERROR\n");
		printf("Should be 'E' or 'H' but has '%c'.\n",CHAMP);
		printf("Abort.\n\n");
		abort();
	}

	
}
