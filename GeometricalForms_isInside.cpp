#include "GeometricalForms_isInside.hpp"

bool is_inside_sphere(
	std::vector<double> const &coord,
	const double radius,
	std::vector<double> const &center)
{
	if(   (coord[0] - center[0])*(coord[0] - center[0]) 
		+ (coord[1] - center[1])*(coord[1] - center[1])
		+ (coord[2] - center[2])*(coord[2] - center[2])
		< radius*radius){
		return true;
	}else{
		return false;
	}
    return false;
}

bool is_inside_cube(
	std::vector<double> const &coord,
	const double side,
	std::vector<double> const &center
){
	if( coord[0] >= (center[0]-side/2) && coord[0] <= (center[0]+side/2)){
		if(coord[1] >= (center[1]-side/2) && coord[1] <= (center[1]+side/2)){
			if(coord[2] >= (center[2]-side/2) && coord[2] <= (center[2]+side/2)){
				return true;
			}
		}
	}
	return false;
}