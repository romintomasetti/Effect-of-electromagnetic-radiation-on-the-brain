#include "vtl.h"
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
#include "vtlSPoints.h"
#include "vtl_spoints.h"



/**
 * @brief export results to paraview (VTK polydata - legacy file fomat)
 * 
 * @param  filename: file name without vtk extension
 * @param  pos:     positions (vector of size 3*number of particles)
 * @param  step:    time step number
 * @param  scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
 * @param   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]
 * @param  binary:   'true' for binary format, 'false' for ASCII
 */
void export_polydata_LEGACY(std::string const &filename,
                            int step,
                            std::vector<double> const &pos,
                            std::map<std::string, std::vector<double> *> const &scalars,
                            std::map<std::string, std::vector<double> *> const &vectors,
                            bool binary, bool verb)
{
    int nbp = (int)pos.size() / 3;
    assert(pos.size() == nbp * 3); // should be multiple of 3

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
    f << "file written by sph.exe\n";
    f << (binary ? "BINARY\n" : "ASCII\n");
    f << "DATASET POLYDATA\n";

    // points
    f << "POINTS " << nbp << " float\n";
    write_vector_LEGACY(f, pos, nbp, 3, binary);

    // vertices
    f << "VERTICES " << nbp << " " << 2 * nbp << "\n";
    if (!binary)
    {
        for (int i = 0; i < nbp; ++i)
            f << "1 " << i << '\n';
        f << '\n'; // empty line (required)
    }
    else
    {
        bool lendian = isCpuLittleEndian();
        int32_t type = lendian ? swap_int32(1) : 1;
        for (int i = 0; i < nbp; ++i)
        {
            int32_t ii = lendian ? swap_int32(i) : i;
            f.write((char *)&type, sizeof(int));
            f.write((char *)&ii, sizeof(int));
        }
    }

    // fields
    f << "POINT_DATA " << nbp << '\n';
    f << "FIELD FieldData " << scalars.size() + vectors.size() << '\n';

    // scalar fields
    std::map<std::string, std::vector<double> *>::const_iterator it = scalars.begin();
    for (; it != scalars.end(); ++it)
    {
        assert(it->second->size() == nbp);
        f << it->first << " 1 " << nbp << " float\n";
        write_vector_LEGACY(f, *it->second, nbp, 1, binary);
    }

    // vector fields
    it = vectors.begin();
    for (; it != vectors.end(); ++it)
    {
        assert(it->second->size() == 3 * nbp);
        f << it->first << " 3 " << nbp << " float\n";
        write_vector_LEGACY(f, *it->second, nbp, 3, binary);
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

void export_polydata_XML(std::string const &filename,
                         int step,
                         std::vector<double> const &pos,
                         std::map<std::string, std::vector<double> *> const &scalars,
                         std::map<std::string, std::vector<double> *> const &vectors,
                         bool binary,
                         bool usez,
                         bool verb)
{
#if !defined(USE_ZLIB)
    if (binary && usez)
    {
        if(verb)
            std::cout << "INFO: zlib not present - vtk file will not be compressed!\n";
        usez = false;
    }
#endif

    int nbp = (int)pos.size() / 3;
    assert(pos.size() == nbp * 3); // should be multiple of 3

    // build file name + stepno + vtk extension
    std::stringstream s;
    s << filename << std::setw(8) << std::setfill('0') << step << ".vtp";
    std::stringstream s2;
    s2 << filename << std::setw(8) << std::setfill('0') << step << ".vtp.tmp";

    // open file
    if(verb)
        std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    std::ofstream f2(s2.str().c_str(), std::ios::binary | std::ios::out); // temp binary file
    f << std::scientific;

    size_t offset = 0;
    // header
    f << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"";
    f << (isCpuLittleEndian() ? "LittleEndian" : "BigEndian") << "\" ";
    f << "header_type=\"UInt32\" "; // UInt64 should be better
    if (usez)
        f << "compressor=\"vtkZLibDataCompressor\" ";
    f << ">\n";
    f << "  <PolyData>\n";
    f << "    <Piece NumberOfPoints=\"" << nbp << "\" ";
    f << "NumberOfVerts=\"" << nbp << "\" ";
    f << "NumberOfLines=\"0\" ";
    f << "NumberOfStrips=\"0\" ";
    f << "NumberOfPolys=\"0\">\n";

    // ------------------------------------------------------------------------------------
    f << "      <PointData>\n";
    // scalar fields
    //for (auto it = scalars.begin(); it != scalars.end(); ++it)
    for (auto const &p : scalars)
    {
        assert(p.second->size() == nbp);
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *p.second, usez);
    }
    // vector fields
    //for (auto it = vectors.begin(); it != vectors.end(); ++it)
    for (auto const &p : vectors)
    {
        assert(p.second->size() == 3 * nbp);
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " NumberOfComponents=\"3\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *p.second, usez);
    }
    f << "      </PointData>\n";

    // ------------------------------------------------------------------------------------
    f << "      <CellData>\n";
    f << "      </CellData>\n";

    // ------------------------------------------------------------------------------------
    f << "      <Points>\n";
    f << "        <DataArray type=\"Float32\" ";
    f << " Name=\"Points\" ";
    f << " NumberOfComponents=\"3\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    offset += write_vectorXML(f2, pos, usez);
    f << "      </Points>\n";
    // ------------------------------------------------------------------------------------
    f << "      <Verts>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"" << nbp - 1 << "\" ";
    f << " offset=\"" << offset << "\" />\n";

    std::vector<int> connectivity(nbp); // <= hard to avoid if zlib is used
    for (int i = 0; i < nbp; ++i)
        connectivity[i] = i;
    offset += write_vectorXML(f2, connectivity, usez);

    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"1\" ";
    f << " RangeMax=\"" << nbp << "\" ";
    f << " offset=\"" << offset << "\" />\n";

    // reuse "connectivity" for offsets
    for (int i = 0; i < nbp; ++i)
        connectivity[i] = i + 1;
    offset += write_vectorXML(f2, connectivity, usez);

    f << "      </Verts>\n";

    std::vector<double> empty;
    // ------------------------------------------------------------------------------------
    f << "      <Lines>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    offset += write_vectorXML(f2, empty, usez);

    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    offset += write_vectorXML(f2, empty, usez);
    f << "      </Lines>\n";

    // ------------------------------------------------------------------------------------
    f << "      <Strips>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    offset += write_vectorXML(f2, empty, usez);
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    offset += write_vectorXML(f2, empty, usez);
    f << "      </Strips>\n";

    // ------------------------------------------------------------------------------------
    f << "      <Polys>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    offset += write_vectorXML(f2, empty, usez);
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    offset += write_vectorXML(f2, empty, usez);
    f << "      </Polys>\n";

    f2.close();

    // ------------------------------------------------------------------------------------
    f << "    </Piece>\n";
    f << "  </PolyData>\n";
    // ------------------------------------------------------------------------------------
    f << "  <AppendedData encoding=\"raw\">\n";
    f << "    _";

    // copy temp binary file
    std::ifstream f3(s2.str().c_str(), std::ios::binary | std::ios::in);
    f << f3.rdbuf();
    f3.close();
    // remove temp file
    std::remove(s2.str().c_str());

    f << "  </AppendedData>\n";
    f << "</VTKFile>\n";

    f.close();
}

// interface

VTL_API void export_polydata(std::string const &filename,
                                  int step,
                                  std::vector<double> const &pos,
                                  std::map<std::string, std::vector<double> *> const &scalars,
                                  std::map<std::string, std::vector<double> *> const &vectors,
                                  PFormat format, bool verb)
{
    switch (format)
    {
    case PFormat::LEGACY_TXT:
        export_polydata_LEGACY(filename, step, pos, scalars, vectors, false, verb);
        break;
    case PFormat::XML_BIN:
        export_polydata_XML(filename, step, pos, scalars, vectors, true, false, verb);
        break;
    case PFormat::XML_BINZ:
        export_polydata_XML(filename, step, pos, scalars, vectors, true, true, verb);
        break;
    case PFormat::LEGACY_BIN:
    default:
        export_polydata_LEGACY(filename, step, pos, scalars, vectors, true, verb);
        break;
    }
}

