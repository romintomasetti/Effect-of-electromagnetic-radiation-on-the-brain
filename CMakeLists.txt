

cmake_minimum_required (VERSION 2.8)


project (Effect-of-electromagnetic-radiation-on-the-brain CXX)

include(CheckCXXCompilerFlag)

set(EXEC_NAME main)

check_cxx_compiler_flag("-std=c++11" STDCXX11)
  if(STDCXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  endif(STDCXX11)

set (CPPFILES main.cpp Array_3D.cpp Materials.cpp GridCreator.cpp MPI_Initializer.cpp
			InputParser.cpp ElectromagneticSource.cpp)
add_executable(${EXEC_NAME} ${CPPFILES})

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
target_link_libraries(${EXEC_NAME} ${MPI_LIBRARIES})
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
