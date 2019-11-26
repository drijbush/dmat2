module multiply_tests
  use test_params
  use mpi
  implicit none

contains
  ! Multiply tests
  Subroutine test_real_pre_scale_real

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision
    Use mpi, Only : mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Cm

    Allocate( A( 1:m, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      C = 3.0_wp * A
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Call Am%create( m, n, base )
    Call Am%set_by_global( 1, m, 1, n, A )
    Cm = 3.0_wp * Am
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real pre_scale real ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_pre_scale_real

  Subroutine test_real_post_scale_real

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision
    Use mpi, Only : mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Cm

    Allocate( A( 1:m, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      C = 3.0_wp * A
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Call Am%create( m, n, base )
    Call Am%set_by_global( 1, m, 1, n, A )
    Cm =  Am * 3.0_wp
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real post_scale real ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_post_scale_real

  Subroutine test_complex_pre_scale_real

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex
    Use mpi, Only : mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Cm

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Allocate( A( 1:m, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
      Call Random_number( rand1 ); Call Random_number( rand2 )
      A = Cmplx( rand1, rand2, wp )
      C = 3.0_wp * A
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Call Am%create( m, n, base )
    Call Am%set_by_global( 1, m, 1, n, A )
    Cm = 3.0_wp * Am
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex pre_scale real ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_pre_scale_real

  Subroutine test_complex_post_scale_real

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex
    Use mpi, Only : mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Cm

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Allocate( A( 1:m, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
      Call Random_number( rand1 ); Call Random_number( rand2 )
      A = Cmplx( rand1, rand2, wp )
      C = 3.0_wp * A
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp, Mold = C )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Call Am%create( m, n, base )
    Call Am%set_by_global( 1, m, 1, n, A )
    Cm = Am * 3.0_wp
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex post_scale real ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_post_scale_real

  Subroutine test_real_matmul_NN

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision
    Use mpi, Only : mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm

    Allocate( A( 1:m, 1:k ) )
    Allocate( B( 1:k, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      Call Random_number( B )
      C = Matmul( A, B )
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
    Call Am%create( m, k, base )
    Call Am%set_by_global( 1, m, 1, k, A )
    Call Bm%create( k, n, Am )
    Call Bm%set_by_global( 1, k, 1, n, B )
    Cm = Am * Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real matmul NN ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_NN

  Subroutine test_real_matmul_TN

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision
    Use mpi, Only : mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, AmT, Bm, Cm

    Allocate( A( 1:k, 1:m ) )
    Allocate( B( 1:k, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      Call Random_number( B )
      C = Matmul( Transpose( A ), B )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp( 1:m, 1:n ) )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Allocate( real_distributed_matrix :: Bm )
    Call Am%create( k, m, base )
    Call Am%set_by_global( 1, k, 1, m, A )
    Call Bm%create( k, n, Am )
    Call Bm%set_by_global( 1, k, 1, n, B )
    AmT = .Dagger. Am 
    Cm  = AmT * Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real matmul TN ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_TN

  Subroutine test_real_matmul_NT

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision
    Use mpi, Only : mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, BmT, Cm

    Allocate( A( 1:m, 1:k ) )
    Allocate( B( 1:n, 1:k ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      Call Random_number( B )
      C = Matmul( A, Transpose( B ) )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp( 1:m, 1:n ) )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Allocate( real_distributed_matrix :: Bm )
    Call Am%create( m, k, base )
    Call Am%set_by_global( 1, m, 1, k, A )
    Call Bm%create( n, k, Am )
    Call Bm%set_by_global( 1, n, 1, k, B )
    BmT = .Dagger. Bm
    Cm = Am * BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real matmul NT ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_NT

  Subroutine test_real_matmul_TT

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision
    Use mpi, Only : mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, AmT, Bm, BmT, Cm

    Allocate( A( 1:k, 1:m ) )
    Allocate( B( 1:n, 1:k ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Call Random_number( A )
      Call Random_number( B )
      C = Matmul( Transpose( A ), Transpose( B ) )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp( 1:m, 1:n ) )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Allocate( real_distributed_matrix :: Bm )
    Call Am%create( k, m, base )
    Call Am%set_by_global( 1, k, 1, m, A )
    Call Bm%create( n, k, Am )
    Call Bm%set_by_global( 1, n, 1, k, B )
    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT * BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in real matmul TT ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_TT

  Subroutine test_complex_matmul_NN

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex
    Use mpi, Only : mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( complex_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, Cm

    Allocate( A( 1:m, 1:k ) )
    Allocate( B( 1:k, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:m, 1:k ), rand2( 1:m, 1:k ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      A = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      Allocate( rand1( 1:k, 1:n ), rand2( 1:k, 1:n ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      B = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      C = Matmul( A, B )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp( 1:m, 1:n ) )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Bm )
    Call Am%create( m, k, base )
    Call Am%set_by_global( 1, m, 1, k, A )
    Call Bm%create( k, n, Am )
    Call Bm%set_by_global( 1, k, 1, n, B )
    Cm = Am * Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex matmul NN ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_NN

  Subroutine test_complex_matmul_TN

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex
    Use mpi, Only : mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( complex_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, AmT, Bm, Cm

    Allocate( A( 1:k, 1:m ) )
    Allocate( B( 1:k, 1:n ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:k, 1:m ), rand2( 1:k, 1:m ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      A = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      Allocate( rand1( 1:k, 1:n ), rand2( 1:k, 1:n ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      B = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      C = Matmul( Conjg( Transpose( A ) ), B )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp( 1:m, 1:n ) )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Bm )
    Call Am%create( k, m, base )
    Call Am%set_by_global( 1, k, 1, m, A )
    Call Bm%create( k, n, Am )
    Call Bm%set_by_global( 1, k, 1, n, B )
    AmT = .Dagger. Am
    Cm  = AmT * Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex matmul TN ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_TN

  Subroutine test_complex_matmul_NT

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex
    Use mpi, Only : mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( complex_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Bm, BmT, Cm

    Allocate( A( 1:m, 1:k ) )
    Allocate( B( 1:n, 1:k ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:m, 1:k ), rand2( 1:m, 1:k ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      A = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      Allocate( rand1( 1:n, 1:k ), rand2( 1:n, 1:k ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      B = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      C = Matmul( A, Conjg( Transpose( B ) ) )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp( 1:m, 1:n ) )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Bm )
    Call Am%create( m, k, base )
    Call Am%set_by_global( 1, m, 1, k, A )
    Call Bm%create( n, k, Am )
    Call Bm%set_by_global( 1, n, 1, k, B )
    BmT = .Dagger. Bm
    Cm = Am * BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex matmul NT ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_NT

  Subroutine test_complex_matmul_TT

!!$    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex
    Use mpi, Only : mpi_comm_world, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, &
      distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
      distributed_matrix_set_default_blocking

    Implicit None

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( complex_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, AmT, Bm, BmT, Cm

    Allocate( A( 1:k, 1:m ) )
    Allocate( B( 1:n, 1:k ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
      Allocate( rand1( 1:k, 1:m ), rand2( 1:k, 1:m ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      A = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      Allocate( rand1( 1:n, 1:k ), rand2( 1:n, 1:k ) )
      Call Random_number( rand1 ); Call Random_number( rand2 ) 
      B = Cmplx( rand1, rand2, wp )
      Deallocate( rand1, rand2 )
      C = Matmul( Conjg( Transpose( A ) ), Conjg( Transpose( B ) ) )
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( B, Size( B ), mpi_double_complex, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex, 0, mpi_comm_world, error )

    Allocate( tmp( 1:m, 1:n ) )

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Allocate( complex_distributed_matrix :: Bm )
    Call Am%create( k, m, base )
    Call Am%set_by_global( 1, k, 1, m, A )
    Call Bm%create( n, k, Am )
    Call Bm%set_by_global( 1, n, 1, k, B )
    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT * BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in complex matmul TT ', Maxval( Abs( C - tmp ) ), &
        Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_TT
end module multiply_tests
