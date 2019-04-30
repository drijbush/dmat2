module subtract_tests
use test_params
implicit none
contains
  ! Subtract tests
  Subroutine test_real_subtract_NN

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm

    Allocate( A( 1:m, 1:n ) )
    Allocate( B( 1:m, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      Call Random_number( B )
      C = A - B
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Allocate( real_distributed_matrix :: Bm )
    Call Am%create( m, n, base )
    Call Am%set_by_global( 1, m, 1, n, A )
    Call Bm%create( m, n, Am )
    Call Bm%set_by_global( 1, m, 1, n, B )
    Cm = Am - Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real subtract NN ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_subtract_NN

  Subroutine test_real_subtract_TN

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm, AmT

    Allocate( A( 1:n, 1:m ) )
    Allocate( B( 1:m, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      Call Random_number( B )
      C = Transpose( A ) - B
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Allocate( real_distributed_matrix :: Bm )
    Call Am%create( n, m, base )
    Call Am%set_by_global( 1, n, 1, m, A )
    Call Bm%create( m, n, Am )
    Call Bm%set_by_global( 1, m, 1, n, B )
    AmT = .Dagger. Am
    Cm = AmT - Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real subtract TN ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_subtract_TN

  Subroutine test_real_subtract_NT

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm, BmT

    Allocate( A( 1:m, 1:n ) )
    Allocate( B( 1:n, 1:m ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      Call Random_number( B )
      C = A - Transpose( B )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Allocate( real_distributed_matrix :: Bm )
    Call Am%create( m, n, base )
    Call Am%set_by_global( 1, m, 1, n, A )
    Call Bm%create( n, m, Am )
    Call Bm%set_by_global( 1, n, 1, m, B )
    BmT = .Dagger. Bm
    Cm = Am - BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real subtract NT ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_subtract_NT

  Subroutine test_real_subtract_TT

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm, AmT, BmT

    Allocate( A( 1:n, 1:m ) )
    Allocate( B( 1:n, 1:m ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      Call Random_number( B )
      C = Transpose( A ) - Transpose( B )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Allocate( real_distributed_matrix :: Bm )
    Call Am%create( n, m, base )
    Call Am%set_by_global( 1, n, 1, m, A )
    Call Bm%create( n, m, Am )
    Call Bm%set_by_global( 1, n, 1, m, B )
    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT - BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real subtract TT ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_subtract_TT

  Subroutine test_complex_subtract_NN

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( complex_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm

    Allocate( A( 1:m, 1:n ) )
    Allocate( B( 1:m, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      A = Cmplx( rand1, rand2, wp )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      B = Cmplx( rand1, rand2, wp )
      C = A - B
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Bm )
    Call Am%create( m, n, base )
    Call Am%set_by_global( 1, m, 1, n, A )
    Call Bm%create( m, n, Am )
    Call Bm%set_by_global( 1, m, 1, n, B )
    Cm = Am - Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex subtract NN ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_subtract_NN

  Subroutine test_complex_subtract_TN

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( complex_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm, AmT

    Allocate( A( 1:n, 1:m ) )
    Allocate( B( 1:m, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      A = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      B = Cmplx( rand1, rand2, wp )
      C = Transpose( Conjg( A ) ) - B
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Bm )
    Call Am%create( n, m, base )
    Call Am%set_by_global( 1, n, 1, m, A )
    Call Bm%create( m, n, Am )
    Call Bm%set_by_global( 1, m, 1, n, B )
    AmT = .Dagger. Am
    Cm = AmT - Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex subtract TN ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_subtract_TN

  Subroutine test_complex_subtract_NT

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( complex_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm, BmT

    Allocate( A( 1:m, 1:n ) )
    Allocate( B( 1:n, 1:m ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      A = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      B = Cmplx( rand1, rand2, wp )
      C = A - Transpose( Conjg( B ) )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Bm )
    Call Am%create( m, n, base )
    Call Am%set_by_global( 1, m, 1, n, A )
    Call Bm%create( n, m, Am )
    Call Bm%set_by_global( 1, n, 1, m, B )
    BmT = .Dagger. Bm
    Cm = Am - BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex subtract NT ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_subtract_NT

  Subroutine test_complex_subtract_TT

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( complex_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm, AmT, BmT

    Allocate( A( 1:n, 1:m ) )
    Allocate( B( 1:n, 1:m ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      A = Cmplx( rand1, rand2, wp )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      B = Cmplx( rand1, rand2, wp )
      C = Transpose( Conjg( A ) ) - Transpose( Conjg( B ) )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Bm )
    Call Am%create( n, m, base )
    Call Am%set_by_global( 1, n, 1, m, A )
    Call Bm%create( n, m, Am )
    Call Bm%set_by_global( 1, n, 1, m, B )
    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT - BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex subtract TT ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_subtract_TT
end module subtract_tests

