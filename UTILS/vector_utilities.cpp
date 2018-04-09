#include "vector_utilities.hpp"

#include "omp.h"

UTILS_API void fill_double_vector_with_zeros(double *vec,size_t size){
	#pragma omp parallel firstprivate(vec)
	{
		#pragma omp for
		for(size_t I = 0 ; I < size ; I ++)
			vec[I] = 0.0;
	}
}
UTILS_API void fill_unchar_vector_with_zeros(unsigned char *vec,size_t size){
	#pragma omp parallel firstprivate(vec)
	{
		#pragma omp for
		for(size_t I = 0 ; I < size ; I ++)
			vec[I] = 0;
	}
}