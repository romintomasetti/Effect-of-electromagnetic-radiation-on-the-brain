#Use a recent version of CMake:
set(MinimalCmakeVersion 3.5.0)
if(${CMAKE_VERSION} VERSION_LESS ${MinimalCmakeVersion}) 
    message(STATUS "Please consider to switch to CMake ${MinimalCmakeVersion}")
	message(STATUS "> You can easily install a second version of cmake in parallel.")
	message(STATUS "> Go to your home directory.")
	message(STATUS "> First, get the sources with: wget https://cmake.org/files/v3.5/cmake-3.5.0.tar.gz")
	message(STATUS "> Then, do: tar xzf cmake-${MinimalCmakeVersion}.tar.gz")
	message(STATUS "> Then, do: cd cmake-${MinimalCmakeVersion}")
	message(STATUS "> Then, do: ./bootstrap --system-curl")
	message(STATUS "> Then, do: gmake")
	message(STATUS "> Note that 'gmake' might take a while, so you can speed it up with gmake - j NBR_PROC")
	message(STATUS "> Use an alias for this new installation (i.e.: $HOME/cmake-${MinimalCmakeVersion}/bin/cmake)")
	message(FATAL_ERROR "Exit.")
endif()

cmake_minimum_required (VERSION ${MinimalCmakeVersion})

#Project name:
project(Effect-of-electromagnetic-radiation-on-the-brain CXX)

#Build documentation:
option(BUILD_DOC "Build documentation" ON)

#Set the build type (Release/Debug):
IF(NOT CMAKE_BUILD_TYPE)
    SET( CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#Set executable name:
set(EXEC_NAME main)

#Add CXX compiler flags:
include(CheckCXXCompilerFlag)
# OPTIMIZATION (only if release mode)
string(COMPARE EQUAL "${CMAKE_BUILD_TYPE}" Release _cmp)
if (_cmp)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
	message(STATUS "Build type is 'Release'. Using -Ofast flag for C++ compiler")
endif()
#Use C++11:
check_cxx_compiler_flag("-std=c++11" STDCXX11)
if(STDCXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(STDCXX11)

#Suppresses a warning on the NIC4 cluster:
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_definitions(-Wno-literal-suffix)
ENDIF()

# ADDING WARNINGS
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long")
endif()

#Add files to be compiled:
FILE(GLOB SRCS *.h *.hpp *.cpp *.c)
ADD_EXECUTABLE(${EXEC_NAME} ${SRCS})

# OpenMP --
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
	MESSAGE(FATAL_ERROR "OpenMP was not found !")
endif()

#Boost:
set(BoostVersion 1.66.0)
set(Boost_VERSION 106600)
find_package(Boost ${BoostVersion} COMPONENTS system filesystem QUIET) 
if(Boost_FOUND)
	message(STATUS "Boost found in ${Boost_INCLUDE_DIRS} / ${BoostSourceDir}")
else()
	################################################################################
	#Install required boost libs:
	include(cmake_modules/add_boost.cmake)
	################################################################################
	# Boost
	set(Boost_USE_STATIC_LIBS OFF) 
	set(Boost_USE_MULTITHREADED ON)  
	set(Boost_USE_STATIC_RUNTIME OFF) 
	message("In main cmakelists boostsourcedir is ${BoostSourceDir}")
	set( BOOST_ROOT "${BoostSourceDir}" )
	set( Boost_NO_SYSTEM_PATHS on CACHE BOOL "Do not search system for Boost" )
	find_package(Boost ${BoostVersion} REQUIRED COMPONENTS system filesystem QUIET) 
	if(NOT Boost_FOUND)
		message(FATAL_ERROR "Boost not found")
	else()
		message(STATUS "After install : Boost ${BoostVersion} found in ${Boost_INCLUDE_DIRS}.")
	endif()
endif()
include_directories(${Boost_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${EXEC_NAME} ${Boost_LIBRARIES})

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

#Zlib:
FIND_PACKAGE(ZLIB)
IF(ZLIB_FOUND)
    add_definitions(-DUSE_ZLIB)
    INCLUDE_DIRECTORIES( ${ZLIB_INCLUDE_DIR} )
    TARGET_LINK_LIBRARIES(${EXEC_NAME} ${ZLIB_LIBRARY} )
ENDIF()


# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${PROJECT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    #message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    # "make dox" rebuilds the doc
    add_custom_target( dox #ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
# ------------------------------------------------------------------------

#Enable testing:
enable_testing()

###########################
## ADDING SUBDIRECTORIES ##
###########################

############ JSON #########
ADD_SUBDIRECTORY( JSON )

find_path(
    RAPIDJSON_INCLUDE_DIRS 
    NAMES "rapidjson/rapidjson.h"
    PATHS ${PROJECT_SOURCE_DIR})
message(STATUS "RAPIDJSON_INCLUDE_DIRS=" ${RAPIDJSON_INCLUDE_DIRS})
if(NOT RAPIDJSON_INCLUDE_DIRS)
    message(FATAL_ERROR "RapidJSON include dir not found!")
endif()
include_directories(${RAPIDJSON_INCLUDE_DIRS})

######## DISC_INTEGR ######
ADD_SUBDIRECTORY( DISC_INTEGR )
INCLUDE_DIRECTORIES( DISC_INTEGR)

########### UTILS #########
ADD_SUBDIRECTORY   ( UTILS )
INCLUDE_DIRECTORIES( UTILS )

##### THERMO ##############
ADD_SUBDIRECTORY(thermo)
INCLUDE_DIRECTORIES(thermo)

##### VTL_BOMAN #################
ADD_SUBDIRECTORY(vtl_boman)
INCLUDE_DIRECTORIES(vtl_boman)

##### VTL_ROMIN #################
ADD_SUBDIRECTORY(vtl_romin)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/vtl_romin)

###########################
## LINKING LIBRARIES ######
###########################
TARGET_LINK_LIBRARIES(${EXEC_NAME} json disc_integr utils vtl vtl_romin)
TARGET_LINK_LIBRARIES(${EXEC_NAME} ${Boost_LIBRARIES})