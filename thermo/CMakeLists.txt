FILE(GLOB SRCS *.h *.cpp)
ADD_EXECUTABLE(TEST_MUMPS ${SRCS})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/vtl_boman)
TARGET_LINK_LIBRARIES(TEST_MUMPS vtl)



# MPI --
# --------------------------------------------------------------
# Finds MPI (including MSMPI)
# --------------------------------------------------------------
if(WIN32)  
    FIND_PATH(MPI_INCLUDE_PATH NAMES  mpi.h HINTS "$ENV{MSMPI_INC}")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
        set(MS_MPI_ARCH_DIR x64)
        find_library(MPI_LIBRARIES msmpi PATHS "$ENV{MSMPI_LIB64}")
    else()
        set(MS_MPI_ARCH_DIR x86)
        find_library(MPI_LIBRARIES msmpi PATHS "$ENV{MSMPI_LIB32}")
    endif()
else()
    find_package(MPI REQUIRED)
endif()
# --------------------------------------------------------------

include_directories(${MPI_INCLUDE_PATH})
#target_link_libraries(TEST_MUMPS ${MPI_LIBRARIES})
IF( MINGW )
    add_definitions(-DMSMPI_NO_SAL)
ENDIF()


# OpenMP --
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    #MESSAGE ("OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# MUMPS -- 

SET(MUMPSDIR $ENV{HOME}/MUMPS_5.1.2)
#SET(MUMPSDIR /home/boman/src/mumps-5.1.2)
#SET(MUMPSDIR C:/mingw-w64/src/mumps-5.1.2)


find_path(MUMPS_INCLUDE_DIRS NAMES "dmumps_c.h" PATHS ${MUMPSDIR}/include NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "MUMPS_INCLUDE_DIRS=" ${MUMPS_INCLUDE_DIRS})
if(NOT MUMPS_INCLUDE_DIRS)
    message(FATAL_ERROR "MUMPS include dir not found!")
endif()
include_directories(${MUMPS_INCLUDE_DIRS})

set(MUMPS_LIBRARIES "")
find_library(DMUMPSLIB  NAMES dmumps PATHS ${MUMPSDIR}/lib NO_SYSTEM_ENVIRONMENT_PATH)
find_library(COMMONLIB  NAMES mumps_common PATHS ${MUMPSDIR}/lib NO_SYSTEM_ENVIRONMENT_PATH)
find_library(PORDLIB  NAMES pord PATHS ${MUMPSDIR}/lib NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "DMUMPSLIB=" ${DMUMPSLIB})
message(STATUS "COMMONLIB=" ${COMMONLIB})
message(STATUS "PORDLIB=" ${PORDLIB})

target_link_libraries(TEST_MUMPS ${DMUMPSLIB} ${COMMONLIB} ${PORDLIB})
# ubuntu
#target_link_libraries(TEST_MUMPS -lgfortran)
#target_link_libraries(TEST_MUMPS -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp -lpthread -lm -ldl)
#target_link_libraries(TEST_MUMPS -pthread -I/usr/lib/openmpi/lib -Wl,-rpath -Wl,/usr/lib/openmpi/lib -Wl,--enable-new-dtags -L/usr/lib/openmpi/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi)

# nic4
target_link_libraries(TEST_MUMPS -lgfortran)
target_link_libraries(TEST_MUMPS  $ENV{MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_gf_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_gnu_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a $ENV{MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lgomp -lpthread -lm -ldl)
target_link_libraries(TEST_MUMPS -pthread -I/cm/shared/apps/openmpi/1.6.4/gcc-4.9.2/lib -L/cm/shared/apps/openmpi/1.6.4/gcc-4.9.2/lib -lmpi_f90 -lmpi_f77 -lmpi -losmcomp -lrdmacm -libverbs -lrt -lnsl -lutil -lm -lpsm_infinipath -ldl -lm -lnuma -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ldl)



# mingw
#target_link_libraries(TEST_MUMPS  "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64/mkl_scalapack_lp64_dll.lib" )
#target_link_libraries(TEST_MUMPS  "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64/mkl_blacs_lp64_dll.lib" )
#target_link_libraries(TEST_MUMPS  "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64/mkl_intel_lp64_dll.lib" )
#target_link_libraries(TEST_MUMPS  "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64/mkl_intel_thread_dll.lib" )
#target_link_libraries(TEST_MUMPS  "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64/mkl_core_dll.lib" )
#target_link_libraries(TEST_MUMPS  "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/compiler/lib/intel64/libiomp5md.lib")
#target_link_libraries(TEST_MUMPS  "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64/msmpifec.lib")
#target_link_libraries(TEST_MUMPS  "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64/msmpi.lib")
#target_link_libraries(TEST_MUMPS  "-lgfortran")




# common
target_link_libraries(TEST_MUMPS ${MPI_LIBRARIES})
