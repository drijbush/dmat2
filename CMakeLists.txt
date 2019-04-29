# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(dmat Fortran) 

set(CMAKE_Fortran_FORMAT FREE)

find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})
find_library(blacs NAMES blacs-openmpi blacs)
find_library(blacsf77 NAMES blacsF77init-openmpi)
find_library(scalapack scalapack-openmpi)
if (NOT EXISTS blacs OR NOT EXISTS blacsf77 OR NOT EXISTS scalapack)
  message(WARNING "Cannot find following parallel libraries:" ${blacs}, ${blacsf77}, ${scalapack})
  return()
endif()


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_DEBUG_FLAGS "-O -g -std=f2008 -Wall -Wextra -Wimplicit-interface -Wuse-without-only -fcheck=all -finit-real=snan" CACHE STRING "Flags used in debug compiles" FORCE)
  set(CMAKE_Fortran_RELEASE_FLAGS "-O3 -funroll-loops" CACHE STRING "Flags used in debug compiles" FORCE)
endif()

add_subdirectory(src)
add_subdirectory(tools)
