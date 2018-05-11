#ifndef JSON_H
#define JSON_H

#if defined(WIN32)
#ifdef json_EXPORTS
#define JSON_API __declspec(dllexport)
#else
#define JSON_API __declspec(dllimport)
#endif
#else
#define JSON_API
#endif

#include "rapidjson/document.h"
#include <string>
#include <vector>

/// reads a file
JSON_API void read_json(std::string fname, rapidjson::Document &d);

/// get parameters from a Document
JSON_API bool read_bool(rapidjson::Document const &d, char const *name, bool def);

JSON_API int read_int(rapidjson::Document const &d, char const *name, int def);

JSON_API double read_double(rapidjson::Document const &d, char const *name, double def);

JSON_API std::string read_string(
    rapidjson::Document const &d, 
    char const *name, 
    std::string const &def);

JSON_API std::vector<double> read_vector_double(
    rapidjson::Document const &d, 
    char const *name, 
    std::vector<double> const &def);

JSON_API std::vector<int>   read_vector_int(
    rapidjson::Document const &d, 
    char const *name, 
    std::vector<int> const &def);

JSON_API std::vector<std::string> read_vector_string(
    rapidjson::Document const &d,
    char const *name,
    std::vector<std::string> const &def
);

JSON_API std::vector<std::size_t> read_vector_size_t(
    rapidjson::Document const &d,
    char const *name,
    std::vector<std::size_t> const &def
);

#endif //JSON_H
