#ifndef GEOMETRICALFORMS_ISINSIDE_HPP
#define GEOMETRICALFORMS_ISINSIDE_HPP


#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>

/**
 * @brief Check is a point is inside a sphere.
 * 
 * Args are: coordinates(x,y,z) of the point
 *           radius of the sphere
 *           center of the sphere (x,y,z)
 */
bool is_inside_sphere(
	std::vector<double> const &coord,
	const double radius,
	std::vector<double> const &center);

bool is_inside_cube(
	std::vector<double> const &coord,
	const double side,
	std::vector<double> const &center
);

#endif // !GEOMETRICALFORMS_ISINSIDE_HPP