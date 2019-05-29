file(GLOB sources *.f90)
add_executable(test_example ${sources})
target_link_libraries(test_example dmat2)

add_test(NAME build_example
  COMMAND ${CMAKE_MAKE_PROGRAM} test_example
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
