#ifndef MISCELLENEOUS_HPP
#define MISCELLENEOUS_HPP

#if defined(WIN32)
#ifdef utils_EXPORTS
#define UTILS_API __declspec(dllexport)
#else
#define UTILS_API __declspec(dllimport)
#endif
#else
#define UTILS_API
#endif

#include <boost/lexical_cast.hpp>
#include <string>

UTILS_API template<typename T> bool isValid(std::string& num) {
   bool flag = true;
   try {
      T tmp = boost::lexical_cast<T>(num);
      tmp += 1;
   }
   catch (boost::bad_lexical_cast &e) {
      flag = false;
   }
   return flag;
}

UTILS_API template <typename T>
size_t GetFileSize(std::string filename)
{
    FILE *f = NULL;
    f = fopen(filename.c_str() , "r");
    if( f == NULL ){
        printf("File %s cannot be opened.\n",filename.c_str());
        abort();
    }
    fseek(f, 0, SEEK_END);
    size_t len = (size_t)(ftell(f)/sizeof(T));
    fclose(f);
    return len;
}

UTILS_API template <typename T >
std::vector<T> readBinaryFile(std::string filename)
{

	std::ifstream in(filename, std::ios::binary | std::ios::in );

    size_t file_size = GetFileSize<T>(filename);
    
    std::vector<T> ret(file_size+1);

    in.read(reinterpret_cast<char*>(&ret[0]), ret.size()*sizeof(T));

    in.close();

    ret.pop_back();
    
    return ret;
}

#endif

