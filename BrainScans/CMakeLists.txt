FILE(GLOB SRCS *.cpp *.hpp)
ADD_EXECUTABLE(brain ${SRCS})

ADD_SUBDIRECTORY( vtl_new_boman )
INCLUDE_DIRECTORIES( vtl_new_boman )
TARGET_LINK_LIBRARIES( brain vtl_new_boman )

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
TARGET_LINK_LIBRARIES(brain ${Boost_LIBRARIES})

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
target_link_libraries(brain ${MPI_LIBRARIES})

