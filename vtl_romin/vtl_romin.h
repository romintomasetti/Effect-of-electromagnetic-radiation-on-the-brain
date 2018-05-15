#ifndef VTL_ROMIN_H
#define VTL_ROMIN_H


#define VTL_API_ROMIN


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

#include "vtlSPoints_romin.h"

#include "GridCreator_NEW.h"

class GridCreator_NEW;

namespace vtl_romin
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


VTL_API_ROMIN void export_polydata(std::string const &filename,
                             int step,
                             std::vector<double> const &pos,
                             std::map<std::string, std::vector<double> *> const &scalars,
                             std::map<std::string, std::vector<double> *> const &vectors,
                             PFormat format);

VTL_API_ROMIN void export_spoints_LEGACY(std::string const &filename,
                                   int step, 
                                   SPoints const &grid,
                                   Mode mode);

VTL_API_ROMIN void export_spoints_XML(std::string const &filename,
                                int step,
                                SPoints const &grid, 
                                SPoints const &mygrid,
                                GridCreator &grid_creator,
                                Zip zip);

/**
 * OUR CUSTOM FUNCTIONS
 */

VTL_API_ROMIN void export_spoints_XMLP_custom_GridCreator(
    std::string type /* THERMAL or ELECTRO */,
    std::string outputFileName,
    size_t currentStep,
    vtl_romin::SPoints &grid,
    vtl_romin::SPoints &my_grid,
    std::vector<vtl_romin::SPoints> &subGrids,
    GridCreator &grid_Creator,
    Zip zip
);

VTL_API_ROMIN void export_spoints_XMLP_custom_GridCreator_NEW(
    std::string type /* THERMAL or ELECTRO */,
    std::string outputFileName,
    size_t currentStep,
    vtl_romin::SPoints &grid,
    vtl_romin::SPoints &my_grid,
    std::vector<vtl_romin::SPoints> &subGrids,
    //GridCreator_NEW &grid_Creator_NEW,
    Zip zip
);


VTL_API_ROMIN void export_spoints_XML_custom_GridCreator_NEW(
    std::string type /* THERMAL or ELECTRO */,
    std::string outputFileName,
    size_t currentStep,
    vtl_romin::SPoints &grid,
    vtl_romin::SPoints &my_grid,
    GridCreator_NEW &grid_Creator_NEW,
    Zip zip,
    algoElectroToVtlRomin pmlVersSauvegarde);

VTL_API_ROMIN void export_spoints_XML_GridCreatorNew(
    std::string const &filename,
    size_t step,
    SPoints const &grid, 
    SPoints const &mygrid,
    GridCreator_NEW &grid_creatorObj,
    Zip zip,
    algoElectroToVtlRomin pmlVersSauvegarde);
    

/**
 * END OF OUR CUSTOM FUNCTIONS
 */

VTL_API_ROMIN void export_spoints_XMLP(std::string const &filename,
                                 int step,
                                 SPoints const &grid,
                                 SPoints const &mygrid,
                                 std::vector<SPoints> const &sgrids,
                                 Zip zip);

}

#endif // VTL_H
