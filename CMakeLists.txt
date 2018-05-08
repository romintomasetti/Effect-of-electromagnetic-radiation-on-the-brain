cmake_minimum_required (FATAL_ERROR VERSION 2.8)

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

string(COMPARE EQUAL "${CMAKE_BUILD_TYPE}" Debug _cmp)
if(_cmp)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
	message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
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
message(${SRCS})

# OpenMP --
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
	MESSAGE(FATAL_ERROR "OpenMP was not found !")
endif()

#Boost:
set(Boost_NO_BOOST_CMAKE ON) #Workaround for old boost versions like the one on NIC4.
find_package(Boost REQUIRED COMPONENTS system filesystem) 
if(Boost_FOUND)
	message(STATUS "Boost found in ${Boost_INCLUDE_DIRS} / ${BoostSourceDir}")
else()
	message(FATAL_ERROR "Cannot find any installation of boost.")
endif()
include_directories(${Boost_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${EXEC_NAME} ${Boost_LIBRARIES})

# vtl_romin.h depends on GridCreator_NEW.h
include_directories(${PROJECT_SOURCE_DIR})

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

##### Brain tests #########
macro (CHECK_IF_DEBUG)
  if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
    get_filename_component(THIS_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    STRING(REGEX REPLACE " " "_" THIS_DIR ${THIS_DIR}) #Replace spaces with underscores
    if (DEFINED DEBUG_${THIS_DIR})
      message(STATUS "Note: Targets in directory ${THIS_DIR} will be built Debug") #A reminder
      set (CMAKE_BUILD_TYPE Debug)
    endif()
  endif()
endmacro()
ADD_SUBDIRECTORY(BrainScans)
INCLUDE_DIRECTORIES(BrainScans)

##### CMake_tests #########
ADD_SUBDIRECTORY(CMake_tests)

##### VTL_BOMAN #################
ADD_SUBDIRECTORY(vtl_boman)
INCLUDE_DIRECTORIES(vtl_boman)

##### CREATE_GEOMETRY ###########
ADD_SUBDIRECTORY(CREATE_GEOMETRY)
INCLUDE_DIRECTORIES(CREATE_GREOMTRY)

##### VTL_ROMIN #################
ADD_SUBDIRECTORY(vtl_romin)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/vtl_romin)

#############################################################
############## THERMO #######################################

# MUMPS -- 
SET(MUMPSDIR /home/users/r/t/rtoma/src/MUMPS_5.1.2)
set(METISDIR /home/users/r/t/rtoma/local/)
set(METIS_LIBRARIES /home/users/r/t/rtoma/local/lib/libparmetis.a /home/users/r/t/rtoma/local/lib/libmetis.a)



find_path(MUMPS_INCLUDE_DIRS NAMES "dmumps_c.h" PATHS ${MUMPSDIR}/include NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "MUMPS_INCLUDE_DIRS=" ${MUMPS_INCLUDE_DIRS})
if(NOT MUMPS_INCLUDE_DIRS)
    message(FATAL_ERROR "MUMPS include dir not found!")
endif()
include_directories(${MUMPS_INCLUDE_DIRS})

find_path(METIS_INCLUDE_DIR NAMES "metis.h" PATHS ${METISDIR}/include NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "METIS_INCLUDE_DIR=" ${METIS_INCLUDE_DIR})
if(NOT METIS_INCLUDE_DIR)
    message(FATAL_ERROR "METIS include dir not found!")
endif()
message(STATUS "METIS_INCLUDE_DIR = ${METIS_INCLUDE_DIR}")
include_directories(${METIS_INCLUDE_DIR})

message(STATUS "METIS_LIBRARIES == ${METIS_LIBRARIES}")


set(MUMPS_LIBRARIES "")
find_library(DMUMPSLIB  NAMES dmumps PATHS ${MUMPSDIR}/lib NO_SYSTEM_ENVIRONMENT_PATH)
find_library(COMMONLIB  NAMES mumps_common PATHS ${MUMPSDIR}/lib NO_SYSTEM_ENVIRONMENT_PATH)
find_library(PORDLIB  NAMES pord PATHS ${MUMPSDIR}/lib NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "DMUMPSLIB=" ${DMUMPSLIB})
message(STATUS "COMMONLIB=" ${COMMONLIB})
message(STATUS "PORDLIB=" ${PORDLIB})

target_link_libraries(${EXEC_NAME} ${DMUMPSLIB} ${COMMONLIB} ${PORDLIB} ${METIS_LIBRARIES})

# nic4
target_link_libraries(${EXEC_NAME} -lgfortran)
target_link_libraries(${EXEC_NAME}  $ENV{MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_gf_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_gnu_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a $ENV{MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lgomp -lpthread -lm -ldl)
target_link_libraries(${EXEC_NAME} -pthread -I/cm/shared/apps/openmpi/1.6.4/gcc-4.9.2/lib -L/cm/shared/apps/openmpi/1.6.4/gcc-4.9.2/lib -lmpi_f90 -lmpi_f77 -lmpi -losmcomp -lrdmacm -libverbs -lrt -lnsl -lutil -lm -lpsm_infinipath -ldl -lm -lnuma -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ldl)

###########################
## LINKING LIBRARIES ######
###########################
TARGET_LINK_LIBRARIES(${EXEC_NAME} json disc_integr utils vtl vtl_romin create_geometry)
TARGET_LINK_LIBRARIES(${EXEC_NAME} ${Boost_LIBRARIES})

