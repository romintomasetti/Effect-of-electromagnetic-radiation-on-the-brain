#include "ElectromagneticSource.h"

void ElectromagneticSource::setLengths(const double L_X, const double L_Y, const double L_Z){
	this->lengthX = L_X;
	this->lengthY = L_Y;
	this->lengthZ = L_Z;
}

void ElectromagneticSource::setCenter (const double C_X, const double C_Y, const double C_Z){
	this->centerX = C_X;
	this->centerY = C_Y;
	this->centerZ = C_Z;
}

void ElectromagneticSource::computeNodesInsideSource(const double L_dom_X,
													 const double L_dom_Y,
													 const double L_dom_Z,
													 const double deltaX,
													 const double deltaY,
													 const double deltaZ){
	// Number of the node of corner 1
	this->nbrNodeCorner1_X = (size_t) (this->centerX - this->lengthX / 2) / deltaX;
	this->nbrNodeCorner1_Y = (size_t) (this->centerY - this->lengthY / 2) / deltaY;
	this->nbrNodeCorner1_Z = (size_t) (this->centerZ - this->lengthZ / 2) / deltaZ;
	
	this->nodesInsideAlong_X = (size_t) this->lengthX / deltaX;
	this->nodesInsideAlong_Y = (size_t) this->lengthY / deltaY;
	this->nodesInsideAlong_Z = (size_t) this->lengthZ / deltaZ;
}

bool ElectromagneticSource::isInsideSource(const size_t x, const size_t y, const size_t z){
	if( ( x >= this->nbrNodeCorner1_X &&
		  x <= (this->nbrNodeCorner1_X + this->nodesInsideAlong_X))
		&&
		( y >= this->nbrNodeCorner1_Y &&
		  y <= (this->nbrNodeCorner1_Y + this->nodesInsideAlong_Y))
		&&
		( z >= this->nbrNodeCorner1_Z &&
		  z <= (this->nbrNodeCorner1_Z + this->nodesInsideAlong_Z))){
		return true;
	}
	// By default, return false:
	return false;
}