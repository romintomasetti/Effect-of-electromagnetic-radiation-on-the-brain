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

UTILS_API std::string vector_string_to_one_string(
	std::vector<std::string> vec, 
	std::string sep_beg,
	std::string sep_end)
{
	std::string merged_str;
	for(size_t I = 0 ; I < vec.size() ; I ++){
		merged_str += sep_beg;
		merged_str += vec[I];
		merged_str += sep_end;
	}
	return merged_str;
}