Program test_distributed_matrix

  Use numbers_module, Only : wp
  Use mpi, Only : mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, mpi_bcast, &
       mpi_comm_world, mpi_integer
  
  Implicit None

  Integer :: me, nproc
  Integer :: n, m, k
  Integer :: n_block
  Integer :: error

  Character( Len = 26 ), Parameter :: format = '( a, t35, g26.20, t65, a )'

  Real( wp ), Parameter :: tol = 1.0e-12_wp
  
  Call mpi_init( error )

  Call mpi_comm_size( mpi_comm_world, nproc, error ) 
  Call mpi_comm_rank( mpi_comm_world,    me, error )

  If( me == 0 ) Then
     Write( *, * ) 'm, n, k ?'
     Read ( *, * )  m, n, k
     Write( *, * ) 'n_block ?'
     Read ( *, * )  n_block
  End If
  Call mpi_bcast( m , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( n , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( k , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( n_block, 1, mpi_integer, 0, mpi_comm_world, error )

  Call test_real_matmul_NN
  Call test_real_matmul_TN
  Call test_real_matmul_NT
  Call test_real_matmul_TT
  Call test_complex_matmul_NN
  Call test_complex_matmul_TN
  Call test_complex_matmul_NT
  Call test_complex_matmul_TT

  Call test_ks_matrix_matmul_real_NN
  Call test_ks_matrix_matmul_real_TN
  Call test_ks_matrix_matmul_real_NT
  Call test_ks_matrix_matmul_real_TT

  Call test_ks_matrix_matmul_complex_NN
  Call test_ks_matrix_matmul_complex_TN
  
  Call mpi_finalize( error )

Contains

  !DISTRIBUTED_MATRIX tests
  
  Subroutine test_real_matmul_NN
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

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
       Write( *, format ) 'Error in real matmul NN ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_NN

  Subroutine test_real_matmul_TN
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

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
       Write( *, format ) 'Error in real matmul TN ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_TN

  Subroutine test_real_matmul_NT
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

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
       Write( *, format ) 'Error in real matmul NT ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_NT

  Subroutine test_real_matmul_TT
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

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
       Write( *, format ) 'Error in real matmul TT ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_TT

  Subroutine test_complex_matmul_NN
    
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
       Write( *, format ) 'Error in complex matmul NN ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_NN

  Subroutine test_complex_matmul_TN
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

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
       Write( *, format ) 'Error in complex matmul TN ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_TN

  Subroutine test_complex_matmul_NT
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

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
       Write( *, format ) 'Error in complex matmul NT ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_NT

  Subroutine test_complex_matmul_TT
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

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
       Write( *, format ) 'Error in complex matmul TT ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_TT

  !KS_MATRIX tests
  
  Subroutine test_ks_matrix_matmul_real_NN

    Use numbers_module  , Only : wp
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_init, ks_matrix_comm_to_base, &
         ks_matrix_finalise
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( ks_matrix ) :: base
    Type( ks_matrix ) :: Am, Bm, Cm

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
    
    Allocate( tmp( 1:m, 1:n ) )

    Call ks_matrix_init( n_block )
    Call ks_matrix_comm_to_base( mpi_comm_world, base )
    Call Am%create( .False., m, k, base )
    Call Am%set_by_global( 1, m, 1, k, A )
    Call Bm%create( .False., k, n, Am )
    Call Bm%set_by_global( 1, k, 1, n, B )
    Cm = Am * Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, format ) 'Error in real ks_matmul NN ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise

  End Subroutine test_ks_matrix_matmul_real_NN

  Subroutine test_ks_matrix_matmul_real_TN

    Use numbers_module  , Only : wp
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_init, ks_matrix_comm_to_base, &
         ks_matrix_finalise
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( ks_matrix ) :: base
    Type( ks_matrix ) :: Am, Bm, Cm, AmT

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

    Call ks_matrix_init( n_block )
    Call ks_matrix_comm_to_base( mpi_comm_world, base )
    Call Am%create( .False., k, m, base )
    Call Am%set_by_global( 1, k, 1, m, A )
    Call Bm%create( .False., k, n, Am )
    Call Bm%set_by_global( 1, k, 1, n, B )
    AmT = .Dagger. Am
    Cm = AmT * Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, format ) 'Error in real ks_matmul TN ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise
    
  End Subroutine test_ks_matrix_matmul_real_TN

  Subroutine test_ks_matrix_matmul_real_NT

    Use numbers_module  , Only : wp
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_init, ks_matrix_comm_to_base, &
         ks_matrix_finalise
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( ks_matrix ) :: base
    Type( ks_matrix ) :: Am, Bm, BmT, Cm

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

    Call ks_matrix_init( n_block )
    Call ks_matrix_comm_to_base( mpi_comm_world, base )
    Call Am%create( .False., m, k, base )
    Call Am%set_by_global( 1, m, 1, k, A )
    Call Bm%create( .False., n, k, Am )
    Call Bm%set_by_global( 1, n, 1, k, B )
    BmT = .Dagger. Bm
    Cm = Am * BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, format ) 'Error in real ks_matmul NT ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise

  End Subroutine test_ks_matrix_matmul_real_NT

  Subroutine test_ks_matrix_matmul_real_TT

    Use numbers_module  , Only : wp
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_init, ks_matrix_comm_to_base, &
         ks_matrix_finalise
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Real( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Type( ks_matrix ) :: base
    Type( ks_matrix ) :: Am, AmT, Bm, BmT, Cm

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

    Call ks_matrix_init( n_block )
    Call ks_matrix_comm_to_base( mpi_comm_world, base )
    Call Am%create( .False., k, m, base )
    Call Am%set_by_global( 1, k, 1, m, A )
    Call Bm%create( .False., n, k, Am )
    Call Bm%set_by_global( 1, n, 1, k, B )
    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT * BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, format ) 'Error in real ks_matmul TT ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise

  End Subroutine test_ks_matrix_matmul_real_TT

  Subroutine test_ks_matrix_matmul_complex_NN

    Use numbers_module  , Only : wp
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_init, ks_matrix_comm_to_base, &
         ks_matrix_finalise
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( ks_matrix ) :: base
    Type( ks_matrix ) :: Am, Bm, Cm

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

    Call ks_matrix_init( n_block )
    Call ks_matrix_comm_to_base( mpi_comm_world, base )
    Call Am%create( .True., m, k, base )
    Call Am%set_by_global( 1, m, 1, k, A )
    Call Bm%create( .True., k, n, Am )
    Call Bm%set_by_global( 1, k, 1, n, B )
    Cm = Am * Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, format ) 'Error in complex ks_matmul NN ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise
    
  End Subroutine test_ks_matrix_matmul_complex_NN

  Subroutine test_ks_matrix_matmul_complex_TN

    Use numbers_module  , Only : wp
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_init, ks_matrix_comm_to_base, &
         ks_matrix_finalise
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( ks_matrix ) :: base
    Type( ks_matrix ) :: Am, Bm, Cm, AmT

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

    Call ks_matrix_init( n_block )
    Call ks_matrix_comm_to_base( mpi_comm_world, base )
    Call Am%create( .True., k, m, base )
    Call Am%set_by_global( 1, k, 1, m, A )
    Call Bm%create( .True., k, n, Am )
    Call Bm%set_by_global( 1, k, 1, n, B )
    AmT = .Dagger. Am
    Cm = AmT * Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, format ) 'Error in complex ks_matmul TN ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise
    
  End Subroutine test_ks_matrix_matmul_complex_TN

End Program test_distributed_matrix
