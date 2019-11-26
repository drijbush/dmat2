module diag_tests
  use test_params
  implicit none
contains

  ! Diag tests
  Subroutine test_diag_real()

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision
    Use mpi, Only : mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Qm, QTm, Cm
    Real( wp ), Dimension( : ) , Allocatable :: E

    Real( wp ), Dimension( :, : ), Allocatable :: A
    Real( wp ), Dimension( :, : ), Allocatable :: tmp

    Integer :: n
    Integer :: i

    n = m

    Allocate( A( 1:n, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      A = A + Transpose( A )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Allocate( real_distributed_matrix :: Qm )
    Call Am%create( n, n, base )
    Call Am%set_by_global( 1, n, 1, n, A )

    Call Am%diag( Qm, E )

    QTm = .Dagger. Qm
    Cm = QTm * Am * Qm 
    Allocate( tmp( 1:n, 1:n ) )
    Call Cm%get_by_global( 1, n, 1, n, tmp )
    Do i = 1, n
      tmp( i, i ) = tmp( i, i ) - E( i )
    End Do

    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real diag ', Maxval( Abs( tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( tmp ) ) < tol )
    End If

  End Subroutine test_diag_real

  Subroutine test_diag_complex()

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex
    Use mpi, Only : mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Qm, QTm, Cm
    Real( wp ), Dimension( : ) , Allocatable :: E

    Complex( wp ), Dimension( :, : ), Allocatable :: A
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Integer :: n
    Integer :: i

    n = m

    Allocate( A( 1:n, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:n, 1:n ), rand2( 1:n, 1:n ) )
      Call Random_number( rand1 ); Call Random_number( rand2 )
      A = Cmplx( rand1, rand2, wp )
      A = A + Conjg( Transpose( A ) )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Qm )
    Call Am%create( n, n, base )
    Call Am%set_by_global( 1, n, 1, n, A )

    Call Am%diag( Qm, E )

    QTm = .Dagger. Qm
    Cm = QTm * Am * Qm 
    Allocate( tmp( 1:n, 1:n ) )
    Call Cm%get_by_global( 1, n, 1, n, tmp )
    Do i = 1, n
      tmp( i, i ) = tmp( i, i ) - E( i )
    End Do

    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex diag ', Maxval( Abs( tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( tmp ) ) < tol )
    End If

  End Subroutine test_diag_complex
end module diag_tests
