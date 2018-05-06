#include "vtl_spoints.h"
#include "vtlSPoints.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <stdint.h>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <cassert>
#include "swapbytes.h"

// --------------------------------------------------------------------------------------------

// export results to paraview (VTK polydata - legacy file fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]
//   binary:   'true' for binary format, 'false' for ASCII

VTL_API void export_spoints_LEGACY(std::string const &filename,
                                        int step, SPoints const &grid,
                                        Mode mode, bool verb)
{
    // build file name + stepno + vtk extension
    std::stringstream s;
    s << filename << std::setw(8) << std::setfill('0') << step << ".vtk";

    // open file
    if(verb)
        std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    f << std::scientific;
    // header
    f << "# vtk DataFile Version 3.0\n";
    f << "file written by vtl\n";
    f << ((mode == Mode::BINARY) ? "BINARY\n" : "ASCII\n");
    f << "DATASET STRUCTURED_POINTS\n";

    // dataset = STRUCTURED_POINTS
    f << "DIMENSIONS " << grid.np()[0] << ' ' << grid.np()[1] << ' ' << grid.np()[2] << '\n';
    f << "SPACING " << grid.dx[0] << ' ' << grid.dx[1] << ' ' << grid.dx[2] << '\n';
    f << "ORIGIN " << grid.o[0] << ' ' << grid.o[1] << ' ' << grid.o[2] << '\n';

    // fields - POINT_DATA
    int nbp = grid.nbp();
    if(nbp!=0)
    {
        f << "POINT_DATA " << nbp << '\n';
        f << "FIELD FieldData " << grid.scalars.size() + grid.vectors.size() << '\n';

        // scalar fields
        for (auto const &p : grid.scalars)
        {
            assert(p.second->size() == nbp);
            f << p.first << " 1 " << nbp << " float\n";
            write_vector_LEGACY(f, *p.second, nbp, 1, (mode == Mode::BINARY));
        }

        // vector fields
        for (auto const &p : grid.vectors)
        {
            assert(p.second->size() == 3 * nbp);
            f << p.first << " 3 " << nbp << " float\n";
            write_vector_LEGACY(f, *p.second, nbp, 3, (mode == Mode::BINARY));
        }
    }
    // fields - CELL_DATA
    int nbc = grid.nbc();
    if(nbc!=0)
    {
        f << "CELL_DATA " << nbc << '\n';
        f << "FIELD FieldData " << grid.cscalars.size() + grid.cvectors.size() << '\n';

        // scalar fields
        for (auto const &p : grid.cscalars)
        {
            assert(p.second->size() == nbc);
            f << p.first << " 1 " << nbc << " float\n";
            write_vector_LEGACY(f, *p.second, nbc, 1, (mode == Mode::BINARY));
        }

        // vector fields
        for (auto const &p : grid.cvectors)
        {
            assert(p.second->size() == 3 * nbc);
            f << p.first << " 3 " << nbc << " float\n";
            write_vector_LEGACY(f, *p.second, nbc, 3, (mode == Mode::BINARY));
        }
    }
    f.close();
}

// export results to paraview (VTK polydata - XML fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]

// see http://www.vtk.org/Wiki/VTK_XML_Formats

VTL_API void export_spoints_XML(std::string const &filename,
                                     int step,
                                     SPoints const &grid, SPoints const &mygrid,
                                     Zip zip, bool verb)
{
#if !defined(USE_ZLIB)
    if (zip == Zip::ZIPPED)
    {
        if(verb)
            std::cout << "INFO: zlib not present - vtk file will not be compressed!\n";
        zip = Zip::UNZIPPED;
    }
#endif

    // build file name (+rankno) + stepno + vtk extension
    std::stringstream s;
    s << filename;
    if (mygrid.id >= 0)
        s << "_r" << mygrid.id;
    s << '_' << std::setw(8) << std::setfill('0') << step << ".vti";
    std::stringstream s2;
    s2 << filename;
    if (mygrid.id >= 0)
        s2 << "r_" << mygrid.id;
    s2 << '_' << std::setw(8) << std::setfill('0') << step << ".vti.tmp";

    // open file
    if(verb)
        std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    std::ofstream f2(s2.str().c_str(), std::ios::binary | std::ios::out); // temp binary file
    f << std::scientific;

    size_t offset = 0;
    // header
    f << "<?xml version=\"1.0\"?>\n";

    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
    f << (isCpuLittleEndian() ? "LittleEndian" : "BigEndian") << "\" ";
    f << "header_type=\"UInt32\" "; // UInt64 could be better (?)
    if (zip == Zip::ZIPPED)
        f << "compressor=\"vtkZLibDataCompressor\" ";
    f << ">\n";

    f << "  <ImageData ";
    f << "WholeExtent=\""
      << grid.np1[0] << ' ' << grid.np2[0] << ' '
      << grid.np1[1] << ' ' << grid.np2[1] << ' '
      << grid.np1[2] << ' ' << grid.np2[2] << "\" ";
    f << "Origin=\"" << grid.o[0] << ' ' << grid.o[1] << ' ' << grid.o[2] << "\" ";
    f << "Spacing=\"" << grid.dx[0] << ' ' << grid.dx[1] << ' ' << grid.dx[2] << "\">\n";

    f << "    <Piece ";
    f << "Extent=\""
      << mygrid.np1[0] << ' ' << mygrid.np2[0] << ' '
      << mygrid.np1[1] << ' ' << mygrid.np2[1] << ' '
      << mygrid.np1[2] << ' ' << mygrid.np2[2] << "\">\n";

    // ------------------------------------------------------------------------------------
    // POINT DATA

    f << "      <PointData>\n";
    // scalar fields
    for (auto const &p : mygrid.scalars)
    {
        //assert(it->second->size() == nbp); // TODO
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *p.second, (zip == Zip::ZIPPED));
    }

    // vector fields
    for (auto const &p : mygrid.vectors)
    {
        //assert(it->second->size() == 3 * nbp); // TODO
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " NumberOfComponents=\"3\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *p.second, (zip == Zip::ZIPPED));
    }
    f << "      </PointData>\n";

    // ------------------------------------------------------------------------------------
    // CELL DATA

    f << "      <CellData>\n";
    // scalar fields
    for (auto const &p : mygrid.cscalars)
    {
        //assert(it->second->size() == nbc); // TODO
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *p.second, (zip == Zip::ZIPPED));
    }

    // vector fields
    for (auto const &p : mygrid.cvectors)
    {
        //assert(it->second->size() == 3 * nbc); // TODO
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " NumberOfComponents=\"3\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *p.second, (zip == Zip::ZIPPED));
    }
    f << "      </CellData>\n";

    f2.close();

    // ------------------------------------------------------------------------------------
    f << "    </Piece>\n";
    f << "  </ImageData>\n";
    // ------------------------------------------------------------------------------------
    f << "  <AppendedData encoding=\"raw\">\n";
    f << "    _";

    // copy temp binary file as "appended" data
    std::ifstream f3(s2.str().c_str(), std::ios::binary | std::ios::in);
    f << f3.rdbuf();
    f3.close();
    // remove temp file
    std::remove(s2.str().c_str());

    f << "  </AppendedData>\n";
    f << "</VTKFile>\n";

    f.close();
}

VTL_API void export_spoints_XMLP(std::string const &filename,
                                      int step,
                                      SPoints const &grid,
                                      SPoints const &mygrid,
                                      std::vector<SPoints> const &sgrids,
                                      Zip zip, bool verb)
{
#if !defined(USE_ZLIB)
    if (zip == ZIPPED)
    {
        if(verb)
            std::cout << "INFO: zlib not present - vtk file will not be compressed!\n";
        zip = UNZIPPED;
    }
#endif

    // build file name (+rankno) + stepno + vtk extension
    std::stringstream s;
    s << filename;
    s << '_' << std::setw(8) << std::setfill('0') << step << ".pvti";

    // open file
    if(verb)
        std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    f << std::scientific;

    // header
    f << "<?xml version=\"1.0\"?>\n";

    f << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"";
    f << (isCpuLittleEndian() ? "LittleEndian" : "BigEndian") << "\" ";
    f << "header_type=\"UInt32\" "; // UInt64 should be better
    if (zip == Zip::ZIPPED)
        f << "compressor=\"vtkZLibDataCompressor\" ";
    f << ">\n";

    f << "  <PImageData ";
    f << "WholeExtent=\""
      << grid.np1[0] << ' ' << grid.np2[0] << ' '
      << grid.np1[0] << ' ' << grid.np2[1] << ' '
      << grid.np1[0] << ' ' << grid.np2[2] << "\" ";
    f << "GhostLevel=\"0\" ";
    f << "Origin=\"" << grid.o[0] << ' ' << grid.o[1] << ' ' << grid.o[2] << "\" ";
    f << "Spacing=\"" << grid.dx[0] << ' ' << grid.dx[1] << ' ' << grid.dx[2] << "\">\n";

    // ------------------------------------------------------------------------------------
    // POINT DATA

    f << "      <PPointData>\n";
    // scalar fields
    for(auto const &p : mygrid.scalars)
    {
        f << "        <PDataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" />\n";
    }
    // vector fields
    for (auto const &p : mygrid.vectors)
    {
        f << "        <PDataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " NumberOfComponents=\"3\" />\n";
    }
    f << "      </PPointData>\n";

    // ------------------------------------------------------------------------------------
    // CELL DATA

    f << "      <PCellData>\n";
    // scalar fields
    for(auto const &p : mygrid.cscalars)
    {
        f << "        <PDataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" />\n";
    }
    // vector fields
    for (auto const &p : mygrid.cvectors)
    {
        f << "        <PDataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " NumberOfComponents=\"3\" />\n";
    }
    f << "      </PCellData>\n";

    // ------------------------------------------------------------------------------------

    for (auto const &g : sgrids)
    {
        f << "    <Piece ";
        f << " Extent=\"";
        f << g.np1[0] << ' ' << g.np2[0] << ' ';
        f << g.np1[1] << ' ' << g.np2[1] << ' ';
        f << g.np1[2] << ' ' << g.np2[2] << "\" ";

        f << "Source=\"";
        std::stringstream s;
        s << filename;
        s << "_r" << g.id;
        s << '_' << std::setw(8) << std::setfill('0') << step << ".vti";
        f << s.str() << "\" />\n";
    }
    // ------------------------------------------------------------------------------------
    f << "  </PImageData>\n";

    f << "</VTKFile>\n";

    f.close();
}
