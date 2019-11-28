# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(dmat Fortran C) 

set(CMAKE_Fortran_FORMAT FREE)

set(CMAKE_MODULE_PATH "${dmat_SOURCE_DIR}/cmake")

find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(ScaLAPACK REQUIRED MODULE)
include_directories(${MPI_Fortran_INCLUDE_PATH})
message(STATUS ${MPI_Fortran_INCLUDE_PATH})
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
#  set(CMAKE_Fortran_FLAGS_DEBUG "-O -g -std=f2008 -Wall -Wextra -Wimplicit-interface -Wuse-without-only -fcheck=all -finit-real=snan" CACHE STRING "Flags used in debug compiles" FORCE)
  set(CMAKE_Fortran_FLAGS_DEBUG "-O -g -std=f2008 -Wall -Wextra -Wuse-without-only -fcheck=all -finit-real=snan" CACHE STRING "Flags used in debug compiles" FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops" CACHE STRING "Flags used in release compiles" FORCE)
endif()
if (CMAKE_Fortran_COMPILER_ID STREQUAL "INTEL")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O -g -stand f08 -check all" CACHE STRING "Flags used in debug compiles" FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "Flags used in release compiles" FORCE)
endif()

set(SOURCE_DIR_NAME "src")
set(TESTS_DIR_NAME "tests")
set(TOOLS_DIR_NAME "tools")
set(EXAMPLES_DIR_NAME "examples")

set(SOURCE_DIR_SOURCE_PATH "${dmat_SOURCE_DIR}/${SOURCE_DIR_NAME}")
set(TESTS_DIR_SOURCE_PATH "${dmat_SOURCE_DIR}/${TESTS_DIR_NAME}")
set(TOOLS_DIR_SOURCE_PATH "${dmat_SOURCE_DIR}/${TOOLS_DIR_NAME}")
set(EXAMPLES_DIR_SOURCE_PATH "${dmat_SOURCE_DIR}/${EXAMPLES_DIR_NAME}")

set(SOURCE_DIR_BINARY_PATH "${dmat_BINARY_DIR}/${SOURCE_DIR_NAME}")
set(TESTS_DIR_BINARY_PATH "${dmat_BINARY_DIR}/${TESTS_DIR_NAME}")
set(TOOLS_DIR_BINARY_PATH "${dmat_BINARY_DIR}/${TOOLS_DIR_NAME}")
set(EXAMPLES_DIR_BINARY_PATH "${dmat_BINARY_DIR}/${EXAMPLES_DIR_NAME}")

set(BIN_OUTPUT_DIR "${dmat_BINARY_DIR}/bin")
set(LIB_OUTPUT_DIR "${dmat_BINARY_DIR}/lib")

file(MAKE_DIRECTORY ${BIN_OUTPUT_DIR} ${LIB_OUTPUT_DIR})
set(EXECUTABLE_OUTPUT_PATH ${BIN_OUTPUT_DIR})
set(LIBRARY_OUTPUT_PATH ${LIB_OUTPUT_DIR})

include_directories(${SOURCE_DIR_BINARY_PATH})
add_subdirectory(src)
enable_testing()
add_subdirectory(tests)
add_subdirectory(tools)
add_subdirectory(examples)
