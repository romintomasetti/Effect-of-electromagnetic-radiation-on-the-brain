#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
#include "readInputGeometryFile.h"

using namespace std;

bool is_inside_sphere(
	size_t I,
	size_t J,
	size_t K,
	size_t S,
	vector<double> radius,
	vector<double> distance_from_origin,
	vector<double> spatial_step)
{
	double coord_X = I * spatial_step[0];
	double coord_Y = J * spatial_step[1];
	double coord_Z = K * spatial_step[2];

	if( (coord_X - distance_from_origin[S+0])*(coord_X - distance_from_origin[S+0]) 
		+ (coord_Y - distance_from_origin[S+1])*(coord_Y - distance_from_origin[S+1])
		+ (coord_Z - distance_from_origin[S+2])*(coord_Z - distance_from_origin[S+2])
		< radius[S]*radius[S]){
		return true;
	}else{
		return false;
	}
}

int main(int argc, char *argv[]){
	/**
	 * Creates an input file with 3 "imbricated" spheres of given radii.
	 */

	/// Size of the domain:
	vector<double> size_dom = {10,10,10};
	vector<double> spatial_step = {0.1,0.1,0.1};
	size_t nodes_X = size_dom[0] / spatial_step[0] + 1;
	size_t nodes_Y = size_dom[1] / spatial_step[1] + 1;
	size_t nodes_Z = size_dom[2] / spatial_step[2] + 1;
	
	/// Radius of the sphere (decreasing order is preferred !)
	vector<double> radius = {4,3,2,1};

	/// Distance of the spheres' center w.r.t. the origin: {x1,y1,z1,x2,y2,z2,x3,y3,z3}
	vector<double> distance_from_origin = {5,5,5,5,5,5,5,5,5,5,5,5};

	if(distance_from_origin.size() != radius.size()*3){
		printf("You have %zu spheres but you provided %zu distances (needs %zu). Aborting.\n",
			radius.size(),distance_from_origin.size(),radius.size()*3);
		return EXIT_FAILURE;
	}

	/// Create the file's name
	string output_name = "imbricated_spheres_";
	output_name.append("_nbr_sph_");
	output_name.append(to_string(radius.size()));
	output_name.append("_R1_");
	output_name.append(to_string((size_t) radius[0]));
	output_name.append("_R2_");
	output_name.append(to_string((size_t) radius[1]));
	output_name.append("_R3_");
	output_name.append(to_string((size_t) radius[2]));
	output_name.append(".geometry");

	/// Open the file:
	FILE* output = NULL;
	if( NULL != (output = fopen(output_name.c_str(),"w+")) ){
		printf("File %s is now opened. Creating...\n",output_name.c_str());
	}else{
		printf("Cannot open the file %s. Aborting.\n",output_name.c_str());
		return EXIT_FAILURE;
	}

	/// Write some information on the first line:
	fprintf(output,"%s\n",output_name.c_str());
	string info_supp = "dx=";
	info_supp.append(to_string(spatial_step[0]));
	info_supp.append(";dy=");
	info_supp.append(to_string(spatial_step[1]));
	info_supp.append(";dz=");
	info_supp.append(to_string(spatial_step[2]));
	info_supp.append(";Lx=");
	info_supp.append(to_string(size_dom[0]));
	info_supp.append(";Ly=");
	info_supp.append(to_string(size_dom[1]));
	info_supp.append(";Lz=");
	info_supp.append(to_string(size_dom[2]));
	fprintf(output,"%s\n",info_supp.c_str());
	

	/// Create the geometry.
	/*
	 * Outside any sphere, we put 0.
	 * Inside the first sphere, we put 1.
	 * Inside the second sphere, we put 2.
	 * Etc...
	 */
	for(size_t I = 0 ; I < radius.size() ; I ++){
		printf(">>> Sphere %zu :: radius %lf :: distance from origin (%lf,%lf,%lf).\n",
			I+1,radius[I],
			distance_from_origin[I+0],
			distance_from_origin[I+1],
			distance_from_origin[I+2]);
	}

	for(size_t K = 0 ; K < nodes_Z ; K ++){
		fprintf(output,"-------\nSLICE %zu\n-------\n",K);
		for(size_t J = 0 ; J < nodes_Y ; J ++){
			for(size_t I = 0 ; I < nodes_X ; I ++){
		
				int value = -1;
				/// Determine if we are in a sphere:
				for(size_t S= 0 ; S < radius.size() ; S ++){
		
					if(is_inside_sphere(I,J,K,S,radius,distance_from_origin,spatial_step)){
						value = (int) S+1;
					}

				}

				/// The point is not inside any sphere, put 0;
				if(value == -1){
					value = 0;
				}

				fprintf(output,"%d",value);

			}
			fprintf(output,"\n");
		}
		fprintf(output,"\n");
	}

	/// Close the file:
	fclose(output);
	printf(">>> File %s is now closed... You can use your new geometry input file !\n",output_name.c_str());

	/// Exit the program:
	return EXIT_SUCCESS;
}
