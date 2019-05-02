add_library(dmat2 ks_array_module.f90 ks_matrix_module.f90 matrix.f90
  matrix_mapping.f90 numbers.f90 proc_mapping.f90
  scalapack_interfaces.f90 blacs_interfaces.f90)
#target_link_libraries(dmat2 ${MPI_Fortran_LIBRARIES} ${ScaLAPACK_Fortran_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(dmat2 ${MPI_Fortran_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})
# message(${LAPACK_LIBRARIES})
# target_link_libraries(dmat2 ${LAPACK_LIBRARIES})#${ScaLAPACK_Fortran_LIBRARIES})
#"scalapack-openmpi" blacs-openmpi blacsF77init-openmpi openblas)
