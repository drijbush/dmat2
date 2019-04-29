
add_executable(dmat2 ks_array_module.f90 ks_matrix_module.f90 main.f90 matrix.f90
  matrix_mapping.f90 numbers.f90 proc_mapping.f90
  scalapack_interfaces.f90 blacs_interfaces.f90)

target_link_libraries(dmat2 ${MPI_Fortran_LIBRARIES})
target_link_libraries(dmat2 blacs)
target_link_libraries(dmat2 blacsf77)
target_link_libraries(dmat2 scalapck)
target_link_libraries(dmat2 ${LAPACK_LIBRARIES})
