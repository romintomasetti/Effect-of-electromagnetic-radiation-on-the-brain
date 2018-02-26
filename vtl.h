#ifndef VTL_H
#define VTL_H


#define VTL_API


#ifdef _MSC_VER
#if !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS 1
#endif
#pragma warning(disable : 4251) // DLL/templates non exportes
#endif                          //_MSC_VER

#include <string>
#include <vector>
#include <map>
#define _USE_MATH_DEFINES // otherwise, M_PI undefined in VS
#include <math.h>

#include "vtlSPoints.h"

#include "GridCreator.h"

namespace vtl
{
class SPoints;

enum PFormat
{
    LEGACY_TXT = 0,
    LEGACY_BIN = 1,
    XML_BIN = 2,
    XML_BINZ = 3
};

enum Zip
{
    UNZIPPED = 0,
    ZIPPED = 1
};

enum Mode
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
                                GridCreator &grid_creator,
                                Zip zip);

VTL_API void export_spoints_XMLP(std::string const &filename,
                                 int step,
                                 SPoints const &grid,
                                 SPoints const &mygrid,
                                 std::vector<SPoints> const &sgrids,
                                 Zip zip);

}

#endif // VTL_H
