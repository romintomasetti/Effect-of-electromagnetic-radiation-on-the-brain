
FILE(GLOB SRCS *.cpp)
ADD_EXECUTABLE(filesCompare ${SRCS})

enable_testing()


################# 1D ABC Face e_x - Electric field along Y #########
add_test(NAME "Test_IS_1D_FACE_EX_Electric_along_Y_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_EX_Electric_along_Y.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_EX_Y/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EX_Electric_along_Y_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face e_x - Electric field along Z #########
add_test(NAME "Test_IS_1D_FACE_EX_Electric_along_Z_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_EX_Electric_along_Z.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_EX_Z/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EX_Electric_along_Z_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face -e_x - Electric field along Y #########
add_test(NAME "Test_IS_1D_FACE_Minus_EX_Electric_along_Y_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_Minus_EX_Electric_along_Y.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_Minus_EX_Y/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EX_Electric_along_Y_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face -e_x - Electric field along Z #########
add_test(NAME "Test_IS_1D_FACE_Minus_EX_Electric_along_Z_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_Minus_EX_Electric_along_Z.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_Minus_EX_Z/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EX_Electric_along_Z_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face e_y - Electric field along X #########
add_test(NAME "Test_IS_1D_FACE_EY_Electric_along_X_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_EY_Electric_along_X.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_EY_X/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EY_Electric_along_X_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face e_y - Electric field along Z #########
add_test(NAME "Test_IS_1D_FACE_EY_Electric_along_Z_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_EY_Electric_along_Z.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_EY_Z/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EY_Electric_along_Z_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face -e_y - Electric field along X #########
add_test(NAME "Test_IS_1D_FACE_Minus_EY_Electric_along_X_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_Minus_EY_Electric_along_X.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_Minus_EY_X/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EY_Electric_along_X_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face -e_y - Electric field along Z #########
add_test(NAME "Test_IS_1D_FACE_Minus_EY_Electric_along_Z_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_Minus_EY_Electric_along_Z.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_Minus_EY_Z/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EY_Electric_along_Z_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face e_z - Electric field along X #########
add_test(NAME "Test_IS_1D_FACE_EZ_Electric_along_X_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_EZ_Electric_along_X.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_EZ_X/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EZ_Electric_along_X_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face e_z - Electric field along Y #########
add_test(NAME "Test_IS_1D_FACE_EZ_Electric_along_Y_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_EZ_Electric_along_Y.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_EZ_Y/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EZ_Electric_along_Y_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face -e_z - Electric field along X #########
add_test(NAME "Test_IS_1D_FACE_Minus_EZ_Electric_along_X_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_Minus_EZ_Electric_along_X.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_Minus_EZ_X/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EZ_Electric_along_X_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D ABC Face -e_z - Electric field along X #########
add_test(NAME "Test_IS_1D_FACE_Minus_EZ_Electric_along_Y_ABC"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TESTS_1D_ABC/FACE_Minus_EZ_Electric_along_Y.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_1D_ABC_Minus_EZ_Y/OUTPUT_FILES_ELECTRO_r0_00003500.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EZ_Electric_along_Y_ABC/OUTPUT_FILES_ELECTRO_r0_00003500.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)


################# 1D PML Face e_x - Electric field along Y - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_EX_Electric_along_Y_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EX_Electric_along_Y_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EX_Y_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EX_Electric_along_Y_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face e_x - Electric field along Z - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_EX_Electric_along_Z_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EX_Electric_along_Z_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EX_Z_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EX_Electric_along_Z_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_x - Electric field along Y - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EX_Electric_along_Y_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EX_Electric_along_Y_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EX_Y_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EX_Electric_along_Y_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_x - Electric field along Z - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EX_Electric_along_Z_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EX_Electric_along_Z_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EX_Z_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EX_Electric_along_Z_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)



################# 1D PML Face e_y - Electric field along X - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_EY_Electric_along_X_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EY_Electric_along_X_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EY_X_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EY_Electric_along_X_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face e_y - Electric field along Z - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_EY_Electric_along_Z_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EY_Electric_along_Z_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EY_Z_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EY_Electric_along_Z_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)


################# 1D PML Face -e_y - Electric field along X - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EY_Electric_along_X_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EY_Electric_along_X_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EY_X_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EY_Electric_along_X_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_y - Electric field along Z - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EY_Electric_along_Z_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EY_Electric_along_Z_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EY_Z_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EY_Electric_along_Z_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face e_z - Electric field along X - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_EZ_Electric_along_X_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EZ_Electric_along_X_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EZ_X_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EZ_Electric_along_X_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face e_z - Electric field along Y - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_EZ_Electric_along_Y_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EZ_Electric_along_Y_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EZ_Y_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EZ_Electric_along_Y_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_z - Electric field along X - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EZ_Electric_along_X_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EZ_Electric_along_X_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EZ_X_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EZ_Electric_along_X_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_z - Electric field along Y - Order 0 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EZ_Electric_along_Y_PML_ORDER_0"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EZ_Electric_along_Y_PML_ORDER_0.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EZ_Y_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EZ_Electric_along_Y_PML_ORDER_0/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)





################# 1D PML Face e_x - Electric field along Y - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_EX_Electric_along_Y_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EX_Electric_along_Y_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EX_Y_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EX_Electric_along_Y_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face e_x - Electric field along Z - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_EX_Electric_along_Z_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EX_Electric_along_Z_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EX_Z_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EX_Electric_along_Z_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_x - Electric field along Y - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EX_Electric_along_Y_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EX_Electric_along_Y_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EX_Y_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EX_Electric_along_Y_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_x - Electric field along Z - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EX_Electric_along_Z_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EX_Electric_along_Z_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EX_Z_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EX_Electric_along_Z_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)



################# 1D PML Face e_y - Electric field along X - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_EY_Electric_along_X_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EY_Electric_along_X_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EY_X_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EY_Electric_along_X_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face e_y - Electric field along Z - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_EY_Electric_along_Z_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EY_Electric_along_Z_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EY_Z_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EY_Electric_along_Z_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)


################# 1D PML Face -e_y - Electric field along X - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EY_Electric_along_X_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EY_Electric_along_X_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EY_X_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EY_Electric_along_X_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_y - Electric field along Z - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EY_Electric_along_Z_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EY_Electric_along_Z_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EY_Z_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EY_Electric_along_Z_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face e_z - Electric field along X - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_EZ_Electric_along_X_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EZ_Electric_along_X_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EZ_X_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EZ_Electric_along_X_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face e_z - Electric field along Y - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_EZ_Electric_along_Y_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_EZ_Electric_along_Y_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_EZ_Y_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_EZ_Electric_along_Y_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_z - Electric field along X - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EZ_Electric_along_X_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EZ_Electric_along_X_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EZ_X_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EZ_Electric_along_X_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)

################# 1D PML Face -e_z - Electric field along Y - Order 2 #########
add_test(NAME "Test_IS_1D_FACE_Minus_EZ_Electric_along_Y_PML_ORDER_2"
         COMMAND ${CMAKE_COMMAND}
         -Dtest_cmd1=mpirun
	 -Dtest_cmd2=-n
	 -Dtest_cmd3=1
	 -Dtest_cmd4=${PROJECT_SOURCE_DIR}/build/main
	 -Dtest_cmd5=-inputfile
	 -Dtest_cmd6=${PROJECT_SOURCE_DIR}/CMake_tests/TESTS/TEST_1D_PML/FACE_Minus_EZ_Electric_along_Y_PML_ORDER_2.input
	 -Dtest_cmd7=-v
	 -Dtest_cmd8=5
         -Doutput_blessed=${PROJECT_SOURCE_DIR}/build/CMake_tests/RESULTS_Minus_EZ_Y_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
         -Doutput_test=${PROJECT_SOURCE_DIR}/CMake_tests/VALIDATING_DATA/IS_1D_FACE_Minus_EZ_Electric_along_Y_PML_ORDER_2/OUTPUT_FILES_ELECTRO_r0_00002000.vti
	 -Dcomparator=${PROJECT_SOURCE_DIR}/build/CMake_tests/filesCompare
         -P ${PROJECT_SOURCE_DIR}/cmake/run_test.cmake
)
