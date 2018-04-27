#ifndef STRINGSEARCH_HPP
#define STRINGSEARCH_HPP

#if defined(WIN32)
#ifdef utils_EXPORTS
#define UTILS_API __declspec(dllexport)
#else
#define UTILS_API __declspec(dllimport)
#endif
#else
#define UTILS_API
#endif

#include <string>

UTILS_API size_t ComputeLevenshteinDistance(std::string str1, std::string str2);

#endif
