#ifndef VTL_SPOINTS_H
#define VTL_SPOINTS_H

#include "vtl.h"

VTL_API void export_spoints_LEGACY(std::string const &filename,
                                   int step, 
                                   SPoints const &grid,
                                   Mode mode, bool verb=true);

VTL_API void export_spoints_XML(std::string const &filename,
                                int step,
                                SPoints const &grid, 
                                SPoints const &mygrid,
                                Zip zip, bool verb=true);

VTL_API void export_spoints_XMLP(std::string const &filename,
                                 int step,
                                 SPoints const &grid,
                                 SPoints const &mygrid,
                                 std::vector<SPoints> const &sgrids,
                                 Zip zip, bool verb=true);

#endif // VTL_SPOINTS_H


