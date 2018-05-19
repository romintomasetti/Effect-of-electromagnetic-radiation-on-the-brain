#ifndef VTL_H
#define VTL_H

#if defined(WIN32)
#ifdef vtl_EXPORTS
#define VTL_API __declspec(dllexport)
#else
#define VTL_API __declspec(dllimport)
#endif
#else
#define VTL_API
#endif

#ifdef _MSC_VER
#if !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS 1
#endif
#pragma warning(disable : 4251) // DLL/templates non exportes
#endif

#include <string>
#include <vector>
#include <fstream>
#include <map>
#define _USE_MATH_DEFINES // otherwise, M_PI undefined in VS
#include <math.h>

//#include "vtlSPoints.h"
/**
 * @brief VTL: "VTK Lite" - tools for exporting results to Paraview 
 */

class SPoints;

/// Format used by the export functions
enum class PFormat
{
    LEGACY_TXT = 0,
    LEGACY_BIN = 1,
    XML_BIN = 2,
    XML_BINZ = 3
};

/// Flag used to enable file compression (requires zlib)
enum class Zip
{
    UNZIPPED = 0,
    ZIPPED = 1
};

/// Text or binary output formats
enum class Mode
{
    TEXT = 0,
    BINARY = 1
};


void write_vector_LEGACY(std::ofstream &f, std::vector<double> const &pos, int nbp, int nj, bool binary);
std::string zlibstatus(int status);
size_t write_vectorXML(std::ofstream &f, std::vector<double> const &pos, bool usez);
size_t write_vectorXML(std::ofstream &f, std::vector<int> const &pos, bool usez);

#endif // VTL_H
