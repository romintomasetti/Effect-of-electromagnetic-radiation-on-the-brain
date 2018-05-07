FILE(GLOB SRCS *.hpp *.cpp *.c *.h)
ADD_LIBRARY(create_geometry SHARED ${SRCS})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
TARGET_LINK_LIBRARIES(create_geometry)

ADD_EXECUTABLE(createInoutFile_3spheres ${SRCS})
