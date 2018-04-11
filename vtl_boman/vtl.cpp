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

using namespace vtl;

#ifdef USE_ZLIB
#include <zlib.h>
#else
#define Z_OK 0
#define uLong size_t
#define uLongf size_t
#endif

static const int __one__ = 1;
static const bool isCpuLittleEndian = 1 == *(char *)(&__one__); // CPU endianness

/**
 * @brief this routine writes a vector of double as a vector of float to a legacy VTK file f.
 * 
 * the vector is converted to float (32bits) and to "big endian" format
 *  (required by the legacy VTK format)
 */

void write_vectorLEGACY(std::ofstream &f, std::vector<double> const &pos, int nbp, int nj, bool binary)
{
    assert(pos.size() == nbp * nj);
    if (!binary)
    {
        for (int i = 0; i < nbp; ++i)
        {
            for (int j = 0; j < nj; ++j)
                f << pos[nj * i + j] << " ";
            f << '\n';
        }
    }
    else
    {
        if (isCpuLittleEndian)
            for (int i = 0; i < nbp; ++i)
            {
                // float+little endian => double should be converted to float, then swapped
                for (int j = 0; j < nj; ++j)
                {
                    float fx = (float)pos[nj * i + j];
                    uint32_t x = swap_uint32(*(uint32_t *)&fx); // convert if CPU is little endian
                    f.write((char *)&x, sizeof(uint32_t));
                }
            }
        else
        {
            // double+bigendian => vector can be written as in memory
            //f.write(reinterpret_cast<char const*>(&(pos[0])), pos.size()*sizeof(double));
            // float+bigendian => vector should be converted to float
            for (int i = 0; i < nbp; ++i)
                for (int j = 0; j < nj; ++j)
                {
                    float fx = (float)pos[nj * i + j];
                    f.write((char *)&fx, sizeof(uint32_t));
                }
        }
    }
}

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
                            bool binary)
{
    int nbp = (int)pos.size() / 3;
    assert(pos.size() == nbp * 3); // should be multiple of 3

    // build file name + stepno + vtk extension
    std::stringstream s;
    s << filename << std::setw(8) << std::setfill('0') << step << ".vtk";

    // open file
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
    write_vectorLEGACY(f, pos, nbp, 3, binary);

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
        int32_t type = isCpuLittleEndian ? swap_int32(1) : 1;
        for (int i = 0; i < nbp; ++i)
        {
            int32_t ii = isCpuLittleEndian ? swap_int32(i) : i;
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
        write_vectorLEGACY(f, *it->second, nbp, 1, binary);
    }

    // vector fields
    it = vectors.begin();
    for (; it != vectors.end(); ++it)
    {
        assert(it->second->size() == 3 * nbp);
        f << it->first << " 3 " << nbp << " float\n";
        write_vectorLEGACY(f, *it->second, nbp, 3, binary);
    }
    f.close();
}

// converts zlib status to a human-readable string

std::string zlibstatus(int status)
{
#ifdef USE_ZLIB
    switch (status)
    {
    case Z_OK:
        return "Z_OK";
    case Z_BUF_ERROR:
        return "Z_BUF_ERROR";
    case Z_MEM_ERROR:
        return "Z_MEM_ERROR";
    case Z_STREAM_ERROR:
        return "Z_STREAM_ERROR";
    default:
        std::stringstream str;
        str << "Unknown (" << status << ")";
        return str.str();
    }
#else
    return "zlib missing";
#endif
}

// sends a vector of "doubles" in binary XML/format into filestream f
//  f: destination filestream
//  pos: the vector to be sent
//  usez: true if zlib should be used

size_t write_vectorXML(std::ofstream &f, std::vector<double> const &pos, bool usez)
{
    size_t written = 0;

    // convert doubles to floats
    std::vector<float> buffer(pos.size());
    for (int i = 0; i < pos.size(); ++i)
        buffer[i] = (float)pos[i];

    if (!usez)
    {
        // data block size
        uint32_t sz = (uint32_t)pos.size() * sizeof(float);
        f.write((char *)&sz, sizeof(uint32_t));
        written += sizeof(uint32_t);
        // data
        f.write((char *)&buffer[0], sz);
        written += sz;
    }
    else
    {
        uLong sourcelen = (uLong)pos.size() * sizeof(float);
        uLongf destlen = uLongf(sourcelen * 1.001) + 12; // see doc
        char *destbuffer = new char[destlen];
#ifdef USE_ZLIB
        int status = compress2((Bytef *)destbuffer, &destlen,
                               (Bytef *)&(buffer[0]), sourcelen, Z_DEFAULT_COMPRESSION);
#else
        int status = Z_OK + 1;
#endif
        if (status != Z_OK)
        {
            std::cout << "ERROR: zlib Error status=" << zlibstatus(status) << "\n";
        }
        else
        {
            //std::cout << "block of size " << sourcelen << " compressed to " << destlen << '\n';
            // blocks description
            uint32_t nblocks = 1;
            f.write((char *)&nblocks, sizeof(uint32_t));
            written += sizeof(uint32_t);
            uint32_t srclen = (uint32_t)sourcelen;
            f.write((char *)&srclen, sizeof(uint32_t));
            written += sizeof(uint32_t);
            uint32_t lastblocklen = 0;
            f.write((char *)&lastblocklen, sizeof(uint32_t));
            written += sizeof(uint32_t);
            uint32_t szblocki = (uint32_t)destlen;
            f.write((char *)&szblocki, sizeof(uint32_t));
            written += sizeof(uint32_t);
            // data
            f.write(destbuffer, destlen);
            written += destlen;
        }

        delete[] destbuffer;
    }

    return written;
}

// sends a vector of "integers" in binary XML/format into filestream f
//  f: destination filestream
//  pos: the vector to be sent
//  usez: true if zlib should be used

// TODO: merge both versions using a template?

size_t write_vectorXML(std::ofstream &f, std::vector<int> const &pos, bool usez)
{
    size_t written = 0;

    if (!usez)
    {
        // data block size
        uint32_t sz = (uint32_t)pos.size() * sizeof(int);
        f.write((char *)&sz, sizeof(uint32_t));
        written += sizeof(uint32_t);
        // data
        f.write((char *)&pos[0], sz);
        written += sz;
    }
    else
    {
        uLong sourcelen = (uLong)pos.size() * sizeof(int);
        uLongf destlen = uLongf(sourcelen * 1.001) + 12; // see doc
        char *destbuffer = new char[destlen];
#ifdef USE_ZLIB
        int status = compress2((Bytef *)destbuffer, &destlen,
                               (Bytef *)&(pos[0]), sourcelen, Z_DEFAULT_COMPRESSION);
#else
        int status = Z_OK + 1;
#endif
        if (status != Z_OK)
        {
            std::cout << "ERROR: zlib Error status=" << zlibstatus(status) << "\n";
        }
        else
        {
            //std::cout << "block of size " << sourcelen << " compressed to " << destlen << '\n';
            // blocks description
            uint32_t nblocks = 1;
            f.write((char *)&nblocks, sizeof(uint32_t));
            written += sizeof(uint32_t);
            uint32_t srclen = (uint32_t)sourcelen;
            f.write((char *)&srclen, sizeof(uint32_t));
            written += sizeof(uint32_t);
            uint32_t lastblocklen = 0;
            f.write((char *)&lastblocklen, sizeof(uint32_t));
            written += sizeof(uint32_t);
            uint32_t szblocki = (uint32_t)destlen;
            f.write((char *)&szblocki, sizeof(uint32_t));
            written += sizeof(uint32_t);
            // data
            f.write(destbuffer, destlen);
            written += destlen;
        }

        delete[] destbuffer;
    }

    return written;
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
                         bool usez)
{
#if !defined(USE_ZLIB)
    if (binary && usez)
    {
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
    std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    std::ofstream f2(s2.str().c_str(), std::ios::binary | std::ios::out); // temp binary file
    f << std::scientific;

    size_t offset = 0;
    // header
    f << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"";
    f << (isCpuLittleEndian ? "LittleEndian" : "BigEndian") << "\" ";
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

VTL_API void vtl::export_polydata(std::string const &filename,
                                  int step,
                                  std::vector<double> const &pos,
                                  std::map<std::string, std::vector<double> *> const &scalars,
                                  std::map<std::string, std::vector<double> *> const &vectors,
                                  PFormat format)
{
    switch (format)
    {
    case PFormat::LEGACY_TXT:
        export_polydata_LEGACY(filename, step, pos, scalars, vectors, false);
        break;
    case PFormat::XML_BIN:
        export_polydata_XML(filename, step, pos, scalars, vectors, true, false);
        break;
    case PFormat::XML_BINZ:
        export_polydata_XML(filename, step, pos, scalars, vectors, true, true);
        break;
    case PFormat::LEGACY_BIN:
    default:
        export_polydata_LEGACY(filename, step, pos, scalars, vectors, true);
        break;
    }
}

// --------------------------------------------------------------------------------------------

// export results to paraview (VTK polydata - legacy file fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]
//   binary:   'true' for binary format, 'false' for ASCII

VTL_API void vtl::export_spoints_LEGACY(std::string const &filename,
                                        int step, SPoints const &grid,
                                        Mode mode)
{
    // build file name + stepno + vtk extension
    std::stringstream s;
    s << filename << std::setw(8) << std::setfill('0') << step << ".vtk";

    // open file
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
    f << "POINT_DATA " << nbp << '\n';
    f << "FIELD FieldData " << grid.scalars.size() + grid.vectors.size() << '\n';

    // scalar fields
    //for (auto it = grid.scalars.begin(); it != grid.scalars.end(); ++it)
    for (auto const &p : grid.scalars)
    {
        assert(p.second->size() == nbp);
        f << p.first << " 1 " << nbp << " float\n";
        write_vectorLEGACY(f, *p.second, nbp, 1, (mode == Mode::BINARY));
    }

    // vector fields
    //for (auto it = grid.vectors.begin(); it != grid.vectors.end(); ++it)
    for (auto const &p : grid.vectors)
    {
        assert(p.second->size() == 3 * nbp);
        f << p.first << " 3 " << nbp << " float\n";
        write_vectorLEGACY(f, *p.second, nbp, 3, (mode == Mode::BINARY));
    }

    // fields - CELL_DATA
    // [TODO]

    f.close();
}

// export results to paraview (VTK polydata - XML fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]

// see http://www.vtk.org/Wiki/VTK_XML_Formats

VTL_API void vtl::export_spoints_XML(std::string const &filename,
                                     int step,
                                     SPoints const &grid, SPoints const &mygrid,
                                     Zip zip)
{
#if !defined(USE_ZLIB)
    if (zip == ZIPPED)
    {
        std::cout << "INFO: zlib not present - vtk file will not be compressed!\n";
        zip = UNZIPPED;
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
    std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    std::ofstream f2(s2.str().c_str(), std::ios::binary | std::ios::out); // temp binary file
    f << std::scientific;

    size_t offset = 0;
    // header
    f << "<?xml version=\"1.0\"?>\n";

    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
    f << (isCpuLittleEndian ? "LittleEndian" : "BigEndian") << "\" ";
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
    f << "      <PointData>\n";

    // scalar fields
    //for (auto it = mygrid.scalars.begin(); it != mygrid.scalars.end(); ++it)
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
    //for (auto it = mygrid.vectors.begin(); it != mygrid.vectors.end(); ++it)
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
    f << "      <CellData>\n";
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

VTL_API void vtl::export_spoints_XMLP(std::string const &filename,
                                      int step,
                                      SPoints const &grid,
                                      SPoints const &mygrid,
                                      std::vector<SPoints> const &sgrids,
                                      Zip zip)
{
#if !defined(USE_ZLIB)
    if (zip == ZIPPED)
    {
        std::cout << "INFO: zlib not present - vtk file will not be compressed!\n";
        zip = UNZIPPED;
    }
#endif

    // build file name (+rankno) + stepno + vtk extension
    std::stringstream s;
    s << filename;
    s << '_' << std::setw(8) << std::setfill('0') << step << ".pvti";

    // open file
    std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    f << std::scientific;

    // header
    f << "<?xml version=\"1.0\"?>\n";

    f << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"";
    f << (isCpuLittleEndian ? "LittleEndian" : "BigEndian") << "\" ";
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
    f << "      <PPointData>\n";
    // scalar fields
    //for (auto it = mygrid.scalars.begin(); it != mygrid.scalars.end(); ++it)
    for(auto const &p : mygrid.scalars)
    {
        f << "        <PDataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" />\n";
    }
    // vector fields
    //for (auto it = mygrid.vectors.begin(); it != mygrid.vectors.end(); ++it)
    for (auto const &p : mygrid.vectors)
    {
        f << "        <PDataArray type=\"Float32\" ";
        f << " Name=\"" << p.first << "\" ";
        f << " NumberOfComponents=\"3\" />\n";
    }
    f << "      </PPointData>\n";

    // ------------------------------------------------------------------------------------
    f << "      <PCellData>\n";
    f << "      </PCellData>\n";

    // ------------------------------------------------------------------------------------

    //for (auto it = sgrids.begin(); it != sgrids.end(); ++it)
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
