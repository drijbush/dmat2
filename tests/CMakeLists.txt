file(GLOB sources *.f90)
add_executable(test_prog ${sources})
target_link_libraries(test_prog dmat2)

add_test(NAME build_tests
  COMMAND ${CMAKE_MAKE_PROGRAM} test_prog
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
foreach(TEST_SET add subtract multiply diag ks_multiply ks_extract ks_unary ks_add ks_subtract ks_multiply ks_diag)
  add_test(${TEST_SET}
    ${BIN_OUTPUT_DIR}/test_prog ${TEST_SET})
  set_tests_properties(${TEST_SET}
    PROPERTIES
    DEPENDS build_tests
    FAIL_REGULAR_EXPRESSION "FAILED")
endforeach()
