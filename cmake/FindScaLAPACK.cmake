# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#.rst:
# FindScaLAPACK
# ----------
#
# Find ScaLAPACK library
#
# This module finds an installed library that implements the
# ScaLAPACK parallel linear-algebra interface (see http://www.netlib.org/lapack/).
#
# The approach follows that taken for the autoconf macro file,
# acx_lapack.m4 (distributed at
# http://ac-archive.sourceforge.net/ac-archive/acx_lapack.html).
#
# This module sets the following variables:
#
# ::
#
#

#Try to use default ScaLAPACK
# find_package(scalapack CONFIG QUIET)
# if (scalapack_FOUND)
#   return()
# endif()

# Check the language being used
if( NOT (CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED OR CMAKE_Fortran_COMPILER_LOADED) )
  if(SCALAPACK_FIND_REQUIRED)
    message(FATAL_ERROR "FindScaLAPACK requires Fortran, C, or C++ to be enabled.")
  else()
    message(STATUS "Looking for ScaLAPACK... - NOT found (Unsupported languages)")
    return()
  endif()
endif()

# Check MPI exists
if (NOT MPI_FOUND)
  if (ScaLAPACK_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for ScaLAPACK... - NOT found MPI")
  else()
    message(STATUS "Looking for ScaLAPACK... - NOT found MPI")
    return()
  endif()
endif()

execute_process(
  COMMAND ${MPIEXEC_EXECUTABLE} --version
  COMMAND sed -r "s/.*(open\\s*mpi|mpich).*/\\1/I"
  OUTPUT_VARIABLE MPI_PROVIDER
  )

string(TOLOWER "${MPI_PROVIDER}" MPI_PROVIDER)

if (MPI_PROVIDER MATCHES "open *mpi|openrte")
  set(MPI_PROVIDER "openmpi")
elseif (MPI_PROVIDER MATCHES "mpich")
  set(MPI_PROVIDER "mpich")
else()
  message(FATAL_ERROR "Unrecognised MPI Provider ${MPI_PROVIDER}")
endif()


find_library(scalapack NAMES scalapack-${MPI_PROVIDER})

if (EXISTS ${scalapack})
  set(SCALAPACK_FOUND TRUE)
  set(SCALAPACK_LIBRARIES ${scalapack})
else()
  set(SCALAPACK_FOUND FALSE)
endif()

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.1.1 .so.1)
endif()


find_library(blacs NAMES blacs-${MPI_PROVIDER})
if (EXISTS ${blacs})
  set(BLACS_FOUND TRUE)
  set(BLACS_LIBRARIES ${blacs})
else()
  set(BLACS_FOUND FALSE)
endif()

set(BLACS_INIT_FOUND TRUE) # Assume true so all must be true

if (CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED)
  find_library(blacsC NAMES blacsCinit-${MPI_PROVIDER})
  if (EXISTS ${blacsC})
    set(BLACS_C_INIT_FOUND TRUE)
    set(BLACS_C_LIBRARIES BLACS_LIBRARIES ${blacsC})
  else()
    set(BLACS_C_INIT_FOUND FALSE)
    set(BLACS_INIT_FOUND FALSE)
  endif()
endif()

if (CMAKE_Fortran_COMPILER_LOADED)
  find_library(blacsFort NAMES blacsF77init-${MPI_PROVIDER})
  if (EXISTS ${blacsFort})
    set(BLACS_Fortran_INIT_FOUND TRUE)
    set(BLACS_Fortran_LIBRARIES ${BLACS_LIBRARIES} ${blacsFort})
  else()
    set(BLACS_Fortran_INIT_FOUND FALSE)
    set(BLACS_INIT_FOUND FALSE)
  endif()
endif()

if (SCALAPACK_FOUND AND BLACS_FOUND AND BLACS_INIT_FOUND)
  set(ScaLAPACK_FOUND TRUE)
  set(ScaLAPACK_C_LIBRARIES ${SCALAPACK_LIBRARIES} ${BLACS_C_LIBRARIES} CACHE STRING "Location of ScaLAPACK C libraries")
  set(ScaLAPACK_Fortran_LIBRARIES ${SCALAPACK_LIBRARIES} ${BLACS_Fortran_LIBRARIES} CACHE STRING "Location of ScaLAPACK C libraries")
  mark_as_advanced(ScaLAPACK_C_LIBRARIES ScaLAPACK_Fortran_LIBRARIES)
  message(STATUS "Found compatible ScaLAPACK Libraries")
else()
  if (SCALAPACK_FIND_REQUIRED)
    message(FATAL_ERROR "ScaLAPACK libraries not found: 
    Found ScaLAPACK : ${SCALAPACK_FOUND}
    Found blacs : ${BLACS_FOUND}
    Found blacs_init : ${BLACS_INIT_FOUND}")
  else()
    message(SEND_ERROR "ScaLAPACK libraries not found: 
    Found ScaLAPACK : ${SCALAPACK_FOUND}
    Found blacs : ${BLACS_FOUND}
    Found blacs_init : ${BLACS_INIT_FOUND}")
    set(ScaLAPACK_FOUND FALSE)
  endif()
endif()
