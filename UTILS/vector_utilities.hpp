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

#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split

UTILS_API void fill_double_vector_with_zeros(double *vec,size_t size);

UTILS_API void fill_unchar_vector_with_zeros(unsigned char *vec,size_t size);

UTILS_API std::string vector_string_to_one_string(
	std::vector<std::string> vec, 
	std::string sep_beg,
	std::string sep_end
);

UTILS_API std::vector<std::string> one_String_to_vector_string(
    std::string strToParse,
    std::string separator,
    bool display
);

UTILS_API void UTILS_parseVectorString_uint(
    std::vector<unsigned int> &vec_1,
    std::vector<std::string> &vec_2,
    std::vector<std::string> const& vec,
    std::string delimiter,
    bool display
);

#endif
