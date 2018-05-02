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

#endif

