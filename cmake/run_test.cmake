# adapted from https://cmake.org/pipermail/cmake/2009-July/030595.html


SET(ENV{CTEST_OUTPUT_ON_FAILURE} 1)


# some argument checking:
# test_cmd is the command to run with all its arguments
if( NOT test_cmd1 )
   message( FATAL_ERROR "Variable test_cmd not defined" )
endif()
# output_blessed contains the name of the "blessed" output file
if( NOT output_blessed )
   message( FATAL_ERROR "Variable output_blessed not defined" )
endif()
# output_test contains the name of the output file the test_cmd will produce
if( NOT output_test )
   message( FATAL_ERROR "Variable output_test not defined" )
endif()

# -------------------------------------------------------------------

# !probleme si output_blessed/output_test comportent des guillemets

message(STATUS "RUNNING ${test_cmd1} ${test_cmd2} ${test_cmd3} ${test_cmd4} ${test_cmd5} ${test_cmd6} ${test_cmd7} ${test_cmd8}")
execute_process(
     COMMAND ${test_cmd1} ${test_cmd2} ${test_cmd3} ${test_cmd4} ${test_cmd5} ${test_cmd6} ${test_cmd7} ${test_cmd8}
     RESULT_VARIABLE test_successful
)
message("test_successful=${test_successful}")
if( NOT test_successful EQUAL 0 )
   message( SEND_ERROR "${test_cmd} did not run properly!" )
endif()

message("COMPARING output files:")
message("\toutput_blessed=${output_blessed}")
message("\toutput_test=${output_test}")

execute_process(
     COMMAND ${comparator} ${output_blessed} ${output_test}
     RESULT_VARIABLE test_successful
)
message("test_successful=${test_successful}")

if( NOT test_successful EQUAL 0 )
   message( SEND_ERROR "${output_test} does not match ${output_blessed}!" )
endif()
