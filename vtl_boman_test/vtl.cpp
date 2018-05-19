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

#ifdef USE_ZLIB
#include <zlib.h>
#else
#define Z_OK 0
#define uLong size_t
#define uLongf size_t
#endif

/**
 * @brief this routine writes a vector of double as a vector of float to a legacy VTK file f.
 * 
 * the vector is converted to float (32bits) and to "big endian" format
 *  (required by the legacy VTK format)
 */

void write_vector_LEGACY(std::ofstream &f, std::vector<double> const &pos, int nbp, int nj, bool binary)
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
        if (isCpuLittleEndian())
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
