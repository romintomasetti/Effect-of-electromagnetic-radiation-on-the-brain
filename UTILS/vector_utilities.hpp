#ifndef VECTOR_UTILITIES_HPP
#define VECTOR_UTILITIES_HPP

#if defined(WIN32)
#ifdef utils_EXPORTS
#define UTILS_API __declspec(dllexport)
#else
#define UTILS_API __declspec(dllimport)
#endif
#else
#define UTILS_API
#endif

#include <stdlib.h>
#include <string>
#include <vector>

UTILS_API void fill_double_vector_with_zeros(double *vec,size_t size);

UTILS_API void fill_unchar_vector_with_zeros(unsigned char *vec,size_t size);

UTILS_API std::string vector_string_to_one_string(
	std::vector<std::string> vec, 
	std::string sep_beg,
	std::string sep_end
);

#endif
