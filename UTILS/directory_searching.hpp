#ifndef DIRECTORY_SEARCHING_HPP
#define DIRECTORY_SEARCHING_HPP

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
#include <iostream>
#include <boost/filesystem.hpp>

UTILS_API int search_for_directory(
    std::string const &dir,
    boost::filesystem::path &foundPath
);

UTILS_API std::string what_is(
    std::string const &name
);

UTILS_API void read_directory_for_files(
  const std::string &name,
  std::vector<std::string> &vec
);

#endif
