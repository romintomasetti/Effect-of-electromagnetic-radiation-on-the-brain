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

#include <sys/types.h>
#include <sys/stat.h>

#include "omp.h"

using namespace vtl;

#ifdef USE_ZLIB
#include <zlib.h>
#else
#define Z_OK 0
#define uLong size_t
#define uLongf size_t
#endif

const int __one__ = 1;
const bool isCpuLittleEndian = 1 == *(char *)(&__one__); // CPU endianness

/**
 * @brief Check that nested omp region is enabled to speed up the buffer filling during writing.
 * 
 * Custom function to ensure that nested parallel regions are enabled.
 * 
 */
void check_omp_nested_enabled(void){
    
	if(const char *omp_dynamic_env = std::getenv("OMP_NESTED")){
		// Already declared. Check it is false.
		if(std::strcmp(omp_dynamic_env,"true") == 0){
			printf("OMP_NESTED=%s.\n",std::getenv("OMP_NESTED"));
		}else{
			std::string set_env = "OMP_NESTED=true";
			putenv(&set_env[0]);
			printf("OMP_NESTED=%s.\n",std::getenv("OMP_NESTED"));
		}
	}else{
		// OMP_DYNAMIC was not declared. Declare it.
		std::string set_env = "OMP_NESTED=true";
		putenv(&set_env[0]);
		printf("OMP_NESTED=%s.\n",std::getenv("OMP_NESTED"));
	}
}


// this routine writes a vector of double as a vector of float to a legacy VTK file f.
//  the vector is converted to float (32bits) and to "big endian" format
//  (required by the legacy VTK format)

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

// export results to paraview (VTK polydata - legacy file fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]
//   binary:   'true' for binary format, 'false' for ASCII

void export_polydata_LEGACY(std::string const &filename,
                            int step,
                            std::vector<double> const &pos,
                            std::map<std::string, std::vector<double> *> const &scalars,
                            std::map<std::string, std::vector<double> *> const &vectors,
                            bool binary)
{
    //std::cout << "system is " << (isCpuLittleEndian? "little" : "big") << " endian\n";
    //bool binary=false;

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

size_t write_vectorXML_custom_GridCreatorNew(
    std::ofstream &f, 
    GridCreator_NEW &grid, 
    std::string fieldName, 
    char vecORsca, 
    bool usez)
{
    check_omp_nested_enabled();
    // Size written in file:
    size_t written = 0;

    // Buffer of floats:
    std::vector<float> buffer;

    // Size of the written field:
    size_t size_field = 0;
    
    if(vecORsca == 's'){
        // Example: temperature field or coordinates x y z.
        if(fieldName == "Temperature"){

            std::vector<size_t> size = grid.size_Thermal;

            printf("MPI %d :: Has field Temperature and size(%ld,,%ld,%ld)\n",
                grid.MPI_communicator.getRank(),size[0],size[1],size[2]);

            size_field = size[0]*size[1]*size[2];
            buffer.resize(size_field);

            size_t counter = 0 ;
            size_t index   = 0;

            size_t I,J,K;

            #pragma omp private(index,I,J,K) shared(buffer,grid) for schedule(static) 
            for(K = 0 ; K < size[2] ; K ++){
                for(J = 0 ; J < size[1] ; J ++){
                    for(I = 0 ; I < size[0] ; I ++){

                        index = I + size[0] * (J + K * size[1]);

                        buffer[index] = (float)grid.temperature[index];
                        
                    }
                }
            }

        }else{
            printf("vtl::write_vectorXML_custom::ERROR in scalar field name. Has %s\n",fieldName.c_str());
            std::abort();
        }
    }else if(vecORsca == 'v'){
        // Example: electric and magnetic fields:
        if(fieldName == "ElectricField"){

            std::cout << "\n\t>>> Writing VTI electric field...\n";

            // Get size M, N, P:
            std::vector<size_t> size = grid.sizes_EH;
            
            // Get the sizes:
            // Ex of size (M − 1) × N × P
            // Ey of size M × (N − 1) × P
            // Ez of size M × N × (P − 1)

            size_field = ( size[0] * size[1] * size[2] ) * 3;
            
            buffer.resize(size_field);

            size_t index   = 0;
            size_t buff_index = 0;

            size_t I,J,K;


            // EX field:
            #pragma omp private(index,I,J,K,buff_index) shared(grid,buffer) for schedule(static) 
            for(K = 1 ; K < grid.size_Ex[2]-1 ; K ++){
                for(J = 1 ; J < grid.size_Ex[1]-1 ; J ++ ){
                    for(I = 1 ; I < grid.size_Ex[0]-1 ; I ++){

                        index = I + grid.size_Ex[0] * ( J + grid.size_Ex[1] * K );

                        buff_index = I-1 + grid.sizes_EH[0] * ( J-1 + grid.sizes_EH[1] * (K-1) );

                        buffer[3*buff_index] = grid.E_x[index];
                        
                    }
                }
            }

            // EY field:
            #pragma omp private(index,I,J,K,buff_index) shared(grid,buffer) for schedule(static) 
            for(size_t K = 1 ; K < grid.size_Ey[2]-1 ; K ++){
                for(size_t J = 1 ; J < grid.size_Ey[1]-1 ; J ++ ){
                    for(size_t I = 1 ; I < grid.size_Ey[0]-1 ; I ++){

                        index = I + grid.size_Ey[0] * ( J + grid.size_Ey[1] * K );

                        buff_index = I-1 + grid.sizes_EH[0] * ( J-1 + grid.sizes_EH[1] * (K-1) );

                        buffer[3*buff_index+1] = grid.E_y[index];
                        
                    }
                }
            }

            // EZ field:
            #pragma omp private(index,I,J,K,buff_index) shared(grid,buffer) for schedule(static) 
            for(size_t K = 1 ; K < grid.size_Ez[2]-1 ; K ++){
                for(size_t J = 1 ; J < grid.size_Ez[1]-1 ; J ++ ){
                    for(size_t I = 1 ; I < grid.size_Ez[0]-1 ; I ++){
                        index = I + grid.size_Ez[0] * ( J + grid.size_Ez[1] * K );

                        buff_index = I-1 + grid.sizes_EH[0] * ( J-1 + grid.sizes_EH[1] * (K-1) );

                        buffer[3*buff_index+2] = grid.E_z[index];
                
                    }
                }
            }

            // END OF ELECTRIC FIELD

        }else if(fieldName == "MagneticField"){

            const size_t DECALAGE_H = 0;

            std::vector<size_t> size = grid.sizes_EH;

            std::cout << "\n\t>>> Writing VTI electric field...\n";

            // Get the sizes:
            // Hx of size M × (N − 1) × (P − 1)
            // Hy of size (M − 1) × N × (P − 1)
            // Hz of size (M − 1) × (N − 1) × P

            size_field = (size[0])*(size[1])*(size[2])*3;

            buffer.resize(size_field);

            size_t index   = 0;
            size_t buff_index = 0;

            size_t I,J,K;

            // HX field:
            #pragma omp private(index,I,J,K,buff_index) shared(grid,buffer) for schedule(static) 
            for(K = 0 ; K < grid.size_Hx[2]-DECALAGE_H ; K ++){
                for(J = 0 ; J < grid.size_Hx[1]-DECALAGE_H ; J ++ ){
                    for(I = 0 ; I < grid.size_Hx[0]-DECALAGE_H ; I ++){

                        index = I + grid.size_Hx[0] * ( J + grid.size_Hx[1] * K );

                        buff_index = I-DECALAGE_H 
                                + grid.sizes_EH[0] * ( J-DECALAGE_H 
                                    + grid.sizes_EH[1] * (K-DECALAGE_H) );

                        buffer[3*buff_index] = grid.H_x[index];
                        
                    }
                }
            }

            // HY field:
            #pragma omp private(index,I,J,K,buff_index) shared(grid,buffer) for schedule(static) 
            for(K = 0 ; K < grid.size_Hy[2]-DECALAGE_H ; K ++){
                for(J = 0 ; J < grid.size_Hy[1]-DECALAGE_H ; J ++ ){
                    for(I = 0 ; I < grid.size_Hy[0]-DECALAGE_H ; I ++){

                        index = I + grid.size_Hy[0] * ( J + grid.size_Hy[1] * K );

                        buff_index = I-DECALAGE_H 
                                + grid.sizes_EH[0] * ( J-DECALAGE_H
                                     + grid.sizes_EH[1] * (K-DECALAGE_H) );

                        buffer[3*buff_index+1] = grid.H_y[index];
                        
                    }
                }
            }

            // HZ field:
            #pragma omp private(index,I,J,K,buff_index) shared(grid,buffer) for schedule(static) 
            for(size_t K = 0 ; K < grid.size_Hz[2]-DECALAGE_H ; K ++){
                for(size_t J = 0 ; J < grid.size_Hz[1]-DECALAGE_H ; J ++ ){
                    for(size_t I = 0 ; I < grid.size_Hz[0]-DECALAGE_H ; I ++){

                        index = I + grid.size_Hz[0] * ( J + grid.size_Hz[1] * K );

                        buff_index = I-DECALAGE_H
                                + grid.sizes_EH[0] * ( J-DECALAGE_H
                                     + grid.sizes_EH[1] * (K-DECALAGE_H) );

                        buffer[3*buff_index+2] = grid.H_z[index];
                        
                    }
                }
            }

            /* END OF MAGNETIC FIELD */

        }else{
            printf("vtl::write_vectorXML_custom::ERROR in vector field name. Has %s\n",fieldName.c_str());
            std::abort();
        }
    }else{
        printf("vtl::write_vectorXML_custom::ERROR on vecORsca\n");
        std::abort();
    }

    if (!usez)
    {
        // data block size
        //uint32_t sz = (uint32_t)pos.size() * sizeof(float);
        uint32_t sz = (uint32_t)size_field * sizeof(float);
        f.write((char *)&sz, sizeof(uint32_t));
        written += sizeof(uint32_t);
        // data
        f.write((char *)&buffer[0], sz);
        written += sz;
    }
    else
    {
        //uLong sourcelen = (uLong)pos.size() * sizeof(float);
        uLong sourcelen = (uLong) size_field * sizeof(float);
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

size_t write_vectorXML_custom(
    std::ofstream &f, 
    GridCreator &grid, 
    std::string fieldName, 
    char vecORsca, 
    bool usez)
{
    size_t written = 0;

    std::vector<float> buffer;
    size_t size_field = 0;
    if(vecORsca == 's'){
        // Example: temperature field or coordinates x y z.
        if(fieldName == "Temperature"){
            std::vector<size_t> size = grid.nodesTemp.get_size_data();
            printf("MPI %d :: Has field Temperature and size(%ld,,%ld,%ld)\n",
                grid.MPI_communicator.getRank(),size[0],size[1],size[2]);
            size_field = size[0]*size[1]*size[2];
            buffer.resize(size_field);
            size_t counter = 0 ;
            for(size_t K = 0 ; K < size[2] ; K ++){
                for(size_t J = 0 ; J < size[1] ; J ++){
                    for(size_t I = 0 ; I < size[0] ; I ++){
                        buffer[counter++] = (float)grid.nodesTemp(I,J,K).field;
                        //printf("nodesTemp(%ld,%ld,%ld)=%f.\n",I,J,K,(float)grid.nodesTemp(I,J,K).field);
                    }
                }
            }

        }else{
            printf("vtl::write_vectorXML_custom::ERROR in scalar field name. Has %s\n",fieldName.c_str());
            std::abort();
        }
    }else if(vecORsca == 'v'){
        // Example: electric and magnetic fields:
        if(fieldName == "ElectricField"){
            std::vector<size_t> size = grid.nodesElec.get_size_data();
            printf("MPI %d :: Has field Electric and size(%ld,%ld,%ld)*3\n",
                grid.MPI_communicator.getRank(),(size[0]-2),size[1]-2,size[2]-2);
            // Do not account for nodes of the neighboors:
            size_field = (size[0]-2)*(size[1]-2)*(size[2]-2)*3;
            buffer.resize(size_field);
            size_t counter = 0;
            for(size_t K = 1 ; K < size[2]-1 ; K ++){
                for(size_t J = 1 ; J < size[1]-1 ; J ++){
                    for(size_t I = 1 ; I < size[0]-1 ; I ++){
                        buffer[counter++] = (float)grid.nodesElec(I,J,K).field[0];
                        buffer[counter++] = (float)grid.nodesElec(I,J,K).field[1];
                        buffer[counter++] = (float)grid.nodesElec(I,J,K).field[2];
                        unsigned long local[3] = {I,J,K};
                        unsigned long global[3] = {0,0,0};
                        grid.LocalToGlobal(local,global);
                        /*printf(">> ELEC(%ld,%ld,%ld) is (%f,%f,%f).\n",global[0],
                            global[1],global[2],
                            (float)grid.nodesElec(I,J,K).field[0],
                            (float)grid.nodesElec(I,J,K).field[1],
                            (float)grid.nodesElec(I,J,K).field[2]);*/
                    }
                }
            }

        }else if(fieldName == "MagneticField"){
            std::vector<size_t> size = grid.nodesMagn.get_size_data();
            printf("MPI %d :: Has field Magentic and size(%ld,,%ld,%ld)*3\n",
                grid.MPI_communicator.getRank(),size[0]-1,size[1]-1,size[2]-1);
            size_field = (size[0]-1)*(size[1]-1)*(size[2]-1)*3;
            buffer.resize(size_field);
            size_t counter = 0;
            for(size_t K = 0 ; K < size[2]-1 ; K ++){
                for(size_t J = 0 ; J < size[1]-1 ; J ++){
                    for(size_t I = 0 ; I < size[0]-1 ; I ++){
                        buffer[counter++] = (float)grid.nodesMagn(I,J,K).field[0];
                        buffer[counter++] = (float)grid.nodesMagn(I,J,K).field[1];
                        buffer[counter++] = (float)grid.nodesMagn(I,J,K).field[2];
                        /*printf(">> MAGN(%ld,%ld,%ld) is (%f,%f,%f).\n",I,J,K,
                            (float)grid.nodesMagn(I,J,K).field[0],
                            (float)grid.nodesMagn(I,J,K).field[1],
                            (float)grid.nodesMagn(I,J,K).field[2]);*/
                    }
                }
            }

        }else{
            printf("vtl::write_vectorXML_custom::ERROR in vector field name. Has %s\n",fieldName.c_str());
            std::abort();
        }
    }else{
        printf("vtl::write_vectorXML_custom::ERROR on vecORsca\n");
        std::abort();
    }

    if (!usez)
    {
        // data block size
        //uint32_t sz = (uint32_t)pos.size() * sizeof(float);
        uint32_t sz = (uint32_t)size_field * sizeof(float);
        f.write((char *)&sz, sizeof(uint32_t));
        written += sizeof(uint32_t);
        // data
        f.write((char *)&buffer[0], sz);
        written += sz;
    }
    else
    {
        //uLong sourcelen = (uLong)pos.size() * sizeof(float);
        uLong sourcelen = (uLong) size_field * sizeof(float);
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
//std::cout << "system is " << (isCpuLittleEndian? "little" : "big") << " endian\n";
//bool binary=false;
//bool usez=true;
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
    for (auto it = scalars.begin(); it != scalars.end(); ++it)
    {
        assert(it->second->size() == nbp);
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << it->first << "\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *it->second, usez);
    }
    // vector fields
    for (auto it = vectors.begin(); it != vectors.end(); ++it)
    {
        assert(it->second->size() == 3 * nbp);
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << it->first << "\" ";
        f << " NumberOfComponents=\"3\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, *it->second, usez);
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
    case LEGACY_TXT:
        export_polydata_LEGACY(filename, step, pos, scalars, vectors, false);
        break;
    case XML_BIN:
        export_polydata_XML(filename, step, pos, scalars, vectors, true, false);
        break;
    case XML_BINZ:
        export_polydata_XML(filename, step, pos, scalars, vectors, true, true);
        break;
    case LEGACY_BIN:
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
    f << ((mode==BINARY) ? "BINARY\n" : "ASCII\n");
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
    for (auto it = grid.scalars.begin(); it != grid.scalars.end(); ++it)
    {
        assert(it->second->size() == nbp);
        f << it->first << " 1 " << nbp << " float\n";
        write_vectorLEGACY(f, *it->second, nbp, 1, (mode==BINARY));
    }

    // vector fields
    for (auto it = grid.vectors.begin(); it != grid.vectors.end(); ++it)
    {
        assert(it->second->size() == 3 * nbp);
        f << it->first << " 3 " << nbp << " float\n";
        write_vectorLEGACY(f, *it->second, nbp, 3, (mode==BINARY));
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
                                SPoints const &grid, 
                                SPoints const &mygrid,
							 	GridCreator &grid_creatorObj,
                                Zip zip)
{
#if !defined(USE_ZLIB)
    if (zip==ZIPPED)
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
    else{
        printf("vtl::export_spoints_XML::ERROR\n");
        printf("\tmygrid.id=%d is smaller than 0. Aborting.\n",mygrid.id);
        abort();
    }
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
    if (zip==ZIPPED)
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
    f << "      <CellData>\n";

    //////////////////////////////////////////////////////////
    //////         OVERWRITTING BOMAN FUNCTION          //////
    //////////////////////////////////////////////////////////
    ////// USING GRIDCREATOR OBJECT INSTREAD OF SPOINTS //////
    // scalar fields
    for (auto it = mygrid.scalars.begin(); it != mygrid.scalars.end(); ++it)
    {
        //assert(it->second->size() == nbp); // TODO
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << it->first << "\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML_custom(f2, grid_creatorObj, it->first ,'s', (zip==ZIPPED));
    }

    // vector fields
    for (auto it = mygrid.vectors.begin(); it != mygrid.vectors.end(); ++it)
    {
        //assert(it->second->size() == 3 * nbp); // TODO
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << it->first << "\" ";
        f << " NumberOfComponents=\"3\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML_custom(f2, grid_creatorObj, it->first, 'v', (zip==ZIPPED));
    }
    f << "      </CellData>\n";

    // ------------------------------------------------------------------------------------
    f << "      <PointData>\n";
    f << "      </PointData>\n";

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
    if (zip==ZIPPED)
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
    if (zip==ZIPPED)
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
    f << "      <PCellData>\n";
    // scalar fields
    if(filename.find("THERMAL") != std::string::npos){
        std::cout << ">>> PVTI :: THERMAL :: IN" << std::endl;
        if(mygrid.scalars.size() == 0){
            fprintf(stderr,"mygrid.vectors.size() == %d\n",mygrid.vectors.size());
            fprintf(stderr,"mygrid.scalars.size() == 0 :: ABORTING\n");
            fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
            abort();
        }
        for (auto it = mygrid.scalars.begin(); it != mygrid.scalars.end(); ++it)
        {
            std::cout << "COUCOU " << it->first << std::endl;
            f << "        <PDataArray type=\"Float32\" ";
            f << " Name=\"" << it->first << "\" />\n";
        }
    }else if(filename.find("ELECTRO") != std::string::npos){
        // vector fields
        for (auto it = mygrid.vectors.begin(); it != mygrid.vectors.end(); ++it)
        {
            std::cout << "COUCOU " << it->first << std::endl;
            f << "        <PDataArray type=\"Float32\" ";
            f << " Name=\"" << it->first << "\" ";
            f << " NumberOfComponents=\"3\" />\n";
        }
    }else{
        fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
        abort();
    }
    f << "      </PCellData>\n";

    // ------------------------------------------------------------------------------------
    f << "      <PPointData>\n";
    f << "      </PPointData>\n";

    // ------------------------------------------------------------------------------------

    for (auto it = sgrids.begin(); it != sgrids.end(); ++it)
    {
        f << "    <Piece ";
        f << " Extent=\"";
        f << it->np1[0] << ' ' << it->np2[0] << ' ';
        f << it->np1[1] << ' ' << it->np2[1] << ' ';
        f << it->np1[2] << ' ' << it->np2[2] << "\" ";

        f << "Source=\"";
        std::stringstream s;
        s << filename;
        
        s << "_r" << it->id;
        s << '_' << std::setw(8) << std::setfill('0') << step << ".vti";
        f << s.str() << "\" />\n";
    }
    // ------------------------------------------------------------------------------------
    f << "  </PImageData>\n";

    f << "</VTKFile>\n";

    f.close();
}



///////////////////////////////////////////////////////////////////////////////
/// CUSTOM FUNCTIONS

/**
 * @brief This function is called by the root MPI process, with GridCreator type.
 */
VTL_API void vtl::export_spoints_XMLP_custom_GridCreator(
    std::string type /* THERMAL or ELECTRO */,
    std::string outputFileName,
    size_t currentStep,
    vtl::SPoints &grid,
    vtl::SPoints &my_grid,
    std::vector<vtl::SPoints> &subGrids,
    GridCreator &grid_Creator,
    Zip zip
){
    if(strcmp(type.c_str(),"ELECTRO") == 0){
        /**
         * @brief The root MPI writes the file for EM grid.
         */
        outputFileName.append("_ELECTRO");
        export_spoints_XMLP(
            outputFileName,
            currentStep,
            grid,
            my_grid,
            subGrids, 
            zip);
    }else if(strcmp(type.c_str(),"THERMAL") == 0){
        /**
         * @brief The root MPI writes the file for the thermal grid.
         */
        outputFileName.append("_THERMAL");
        export_spoints_XMLP(
            outputFileName,
            currentStep,
            grid,
            my_grid,
            subGrids, 
            zip);
    }else{
        fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
        abort();
    }
}

/**
 * @brief This function is called by the root MPI process, with GridCreator_NEW type.
 */
VTL_API void vtl::export_spoints_XMLP_custom_GridCreator_NEW(
    std::string type /* THERMAL or ELECTRO */,
    std::string outputFileName,
    size_t currentStep,
    vtl::SPoints &grid,
    vtl::SPoints &my_grid,
    std::vector<vtl::SPoints> &subGrids,
    GridCreator_NEW &grid_Creator_NEW,
    Zip zip
){
    if(strcmp(type.c_str(),"ELECTRO") == 0){
        /**
         * @brief The root MPI writes the file for EM grid.
         */
        outputFileName.append("_ELECTRO");
        export_spoints_XMLP(
            outputFileName,
            currentStep,
            grid,
            my_grid,
            subGrids, 
            zip);
    }else if(strcmp(type.c_str(),"THERMAL") == 0){
        /**
         * @brief The root MPI writes the file for the thermal grid.
         */
        outputFileName.append("_THERMAL");
        std::cout << "Call::export_spoints_XMLP::GridCreator_NEW::IN" << std::endl;
        export_spoints_XMLP(
            outputFileName,
            currentStep,
            grid,
            my_grid,
            subGrids, 
            zip);
    }else{
        fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
        abort();
    }
}

/**
 * @brief This function is called by all MPI process, with GridCreator type.
 */
VTL_API void vtl::export_spoints_XML_custom_GridCreator(
    std::string type /* THERMAL or ELECTRO */,
    std::string outputFileName,
    size_t currentStep,
    vtl::SPoints &grid,
    vtl::SPoints &my_grid,
    GridCreator &grid_creator,
    Zip zip
){
    if(strcmp(type.c_str(),"ELECTRO") == 0){
        /**
         * @brief The MPI process writes its electromagnetic grid.
         */
        fprintf(stderr,"In %s: not yet coded. We don't need that any more !\n",__FUNCTION__);
        fprintf(stderr,"In file %s:%d\n",__FILE__,__LINE__);
        #ifdef MPI_COMM_WORLD
        MPI_Abort(MPI_COMM_WORLD,-1);
        #else
        abort();
        #endif
    }else if(strcmp(type.c_str(),"THERMAL") == 0){
        /**
         * @brief The MPI process writes its thermal grid.
         */
        fprintf(stderr,"In %s: not yet coded. We don't need that any more !\n",__FUNCTION__);
        fprintf(stderr,"In file %s:%d\n",__FILE__,__LINE__);
        #ifdef MPI_COMM_WORLD
        MPI_Abort(MPI_COMM_WORLD,-1);
        #else
        abort();
        #endif
    }else{
        fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
        abort();
    }
}

/**
 * @brief This function is called by all MPI process, with GridCreator type.
 * 
 */
VTL_API void vtl::export_spoints_XML_custom_GridCreator_NEW(
    std::string type /* THERMAL or ELECTRO */,
    std::string outputFileName,
    size_t currentStep,
    vtl::SPoints &grid,
    vtl::SPoints &my_grid,
    GridCreator_NEW &grid_Creator_NEW,
    Zip zip
){
    if(strcmp(type.c_str(),"ELECTRO") == 0){
        /**
         * @brief The MPI process writes its electromagnetic grid.
         */
        outputFileName.append("_ELECTRO");
        export_spoints_XML_GridCreatorNew(
            outputFileName,
            currentStep,
            grid, 
            my_grid,
            grid_Creator_NEW,
            zip);
    }else if(strcmp(type.c_str(),"THERMAL") == 0){
        /**
         * @brief The MPI process writes its thermal grid.
         */
        outputFileName.append("_THERMAL");
        export_spoints_XML_GridCreatorNew(
            outputFileName,
            currentStep,
            grid, 
            my_grid,
            grid_Creator_NEW,
            zip);
    }else{
        fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
        abort();
    }
}

VTL_API void vtl::export_spoints_XML_GridCreatorNew(
    std::string const &filename,
    size_t step,
    SPoints const &grid, 
    SPoints const &mygrid,
    GridCreator_NEW &grid_creatorObj,
    Zip zip)
{
#if !defined(USE_ZLIB)
    if (zip==ZIPPED)
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
    else{
        printf("vtl::export_spoints_XML::ERROR\n");
        printf("\tmygrid.id=%d is smaller than 0. Aborting.\n",mygrid.id);
        abort();
    }
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
    if (zip==ZIPPED)
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
    f << "      <CellData>\n";

    //////////////////////////////////////////////////////////
    //////         OVERWRITTING BOMAN FUNCTION          //////
    //////////////////////////////////////////////////////////
    ////// USING GRIDCREATOR OBJECT INSTREAD OF SPOINTS //////
    // scalar fields
    if(filename.find("THERMAL") != std::string::npos){
        for (auto it = mygrid.scalars.begin(); it != mygrid.scalars.end(); ++it)
        {
            //assert(it->second->size() == nbp); // TODO
            f << "        <DataArray type=\"Float32\" ";
            f << " Name=\"" << it->first << "\" ";
            f << " format=\"appended\" ";
            f << " RangeMin=\"0\" ";
            f << " RangeMax=\"1\" ";
            f << " offset=\"" << offset << "\" />\n";
            offset += write_vectorXML_custom_GridCreatorNew(
                f2, grid_creatorObj, it->first ,'s', (zip==ZIPPED));
        }
    }else if(filename.find("ELECTRO") != std::string::npos){

        // vector fields
        for (auto it = mygrid.vectors.begin(); it != mygrid.vectors.end(); ++it)
        {
            //assert(it->second->size() == 3 * nbp); // TODO
            f << "        <DataArray type=\"Float32\" ";
            f << " Name=\"" << it->first << "\" ";
            f << " NumberOfComponents=\"3\" ";
            f << " format=\"appended\" ";
            f << " RangeMin=\"0\" ";
            f << " RangeMax=\"1\" ";
            f << " offset=\"" << offset << "\" />\n";
            offset += write_vectorXML_custom_GridCreatorNew(f2, 
                        grid_creatorObj, 
                        it->first, 'v', 
                        (zip==ZIPPED));
        }
    }else{
        fprintf(stderr,"Cannot find thermal or electro in the filename (has %s)\n",
            filename.c_str());
        fprintf(stderr,"File %s:%d\n",__FILE__,__LINE__);
        abort();
    }

    f << "      </CellData>\n";

    // ------------------------------------------------------------------------------------
    f << "      <PointData>\n";
    f << "      </PointData>\n";

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

