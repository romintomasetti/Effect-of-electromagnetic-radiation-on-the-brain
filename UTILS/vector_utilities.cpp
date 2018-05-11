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

UTILS_API std::vector<std::string> one_String_to_vector_string(
    std::string strToParse,
    std::string delimiter,
    bool display
)
{
    std::vector<std::string> splitedStr;
    boost::split(splitedStr, strToParse, boost::is_any_of(delimiter), boost::token_compress_on);
    if(display){
        for(size_t I = 0 ; I < splitedStr.size() ; I ++)
            printf("[..%s..]",splitedStr[I].c_str());
        printf("\n");
    }
    return splitedStr;
}

UTILS_API void UTILS_parseVectorString_uint(
    std::vector<unsigned int> &vec_1,
    std::vector<std::string> &vec_2,
    std::vector<std::string> const& vec,
    std::string delimiter,
    bool display
)
{
    for(size_t I = 0 ; I < vec.size() ; I ++){
        // Find separator:
        size_t pos_del = vec[I].find(delimiter);
        std::string temp1 = vec[I].substr(0,pos_del);
        std::string temp2 = vec[I].substr(pos_del+1,vec[I].length()-pos_del);
        printf("[%s] ==> [%s] + [%s]\n",vec[I].c_str(),temp1.c_str(),temp2.c_str());
        vec_1.push_back(std::stoul(temp1));
        vec_2.push_back(temp2);
    }
}