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
#include <map>
#define _USE_MATH_DEFINES // otherwise, M_PI undefined in VS
#include <math.h>

#include "vtlSPoints.h"
/**
 * @brief VTL: "VTK Lite" - tools for exporting results to Paraview 
 */
namespace vtl
{
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


VTL_API void export_polydata(std::string const &filename,
                             int step,
                             std::vector<double> const &pos,
                             std::map<std::string, std::vector<double> *> const &scalars,
                             std::map<std::string, std::vector<double> *> const &vectors,
                             PFormat format);

VTL_API void export_spoints_LEGACY(std::string const &filename,
                                   int step, 
                                   SPoints const &grid,
                                   Mode mode);

VTL_API void export_spoints_XML(std::string const &filename,
                                int step,
                                SPoints const &grid, 
                                SPoints const &mygrid,
                                Zip zip);

VTL_API void export_spoints_XMLP(std::string const &filename,
                                 int step,
                                 SPoints const &grid,
                                 SPoints const &mygrid,
                                 std::vector<SPoints> const &sgrids,
                                 Zip zip);

}

#endif // VTL_H
