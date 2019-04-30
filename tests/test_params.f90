module test_params

  use numbers_module, only : wp
  implicit none
  Logical, Parameter :: verbose = .False.

  Integer :: me, nproc
  Integer :: n, m, k
  Integer :: n_block
  Integer :: ns, nk
  Integer :: error

  Integer, parameter :: n_def = 1000
  Integer, parameter :: m_def = 1000
  Integer, parameter :: k_def = 1000
  Integer, parameter :: n_block_def = 100
  Integer, parameter :: ns_def = 100
  integer, parameter :: nk_def = 100
  
  Character( Len = * ), Parameter :: error_format = '( "--> ", a, t40, g26.20, t70, a )'
  Character( Len = * ), Parameter ::   run_format = '( a, t20, i0 )'
  Character( Len = * ), Parameter :: title_format = '( t5, a )'
  Character( Len = * ), Parameter :: passed       = 'Passed      '
  Character( Len = * ), Parameter :: FAILED       = '      FAILED'
  
  Real( wp ), Parameter :: tol = 1.0e-11_wp

end module test_params

