
FILE(GLOB SRCS *.cpp)
ADD_EXECUTABLE(filesCompare ${SRCS})

enable_testing()


add_test(NAME "Test_IS_1D_FACE_EX_Electric_along_Y"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_EX_Electric_along_Y.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS/OUTPUT_FILES_ELECTRO_r0_00000020.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EX_Electric_along_Y/OUTPUT_FILES_ELECTRO_r0_00000020.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)


#
