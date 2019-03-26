Program test_distributed_matrix

  Use numbers_module, Only : wp
  Use mpi, Only : mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, mpi_bcast, &
       mpi_comm_world, mpi_integer
  
  Implicit None

  Integer :: me, nproc
  Integer :: n, m, k
  Integer :: n_block
  Integer :: ns, nk
  Integer :: error

  Character( Len = 34 ), Parameter :: format = '( "--> ", a, t35, g26.20, t65, a )'

  Real( wp ), Parameter :: tol = 1.0e-11_wp
  
  Call mpi_init( error )

  Call mpi_comm_size( mpi_comm_world, nproc, error ) 
  Call mpi_comm_rank( mpi_comm_world,    me, error )

  If( me == 0 ) Then
     Write( *, * ) 'm, n, k ?'
     Read ( *, * )  m, n, k
     Write( *, * ) 'n_block ?'
     Read ( *, * )  n_block
     Write( *, * ) 'ns, nk ?'
     Read ( *, * )  ns, nk
  End If
  Call mpi_bcast( m , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( n , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( k , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( n_block, 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( ns, 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( nk, 1, mpi_integer, 0, mpi_comm_world, error )

  ! Distributed matrix tests
  Call test_real_matmul_NN
  Call test_real_matmul_TN
  Call test_real_matmul_NT
  Call test_real_matmul_TT
  Call test_complex_matmul_NN
  Call test_complex_matmul_TN
  Call test_complex_matmul_NT
  Call test_complex_matmul_TT

  ! KS_matrix tests
  Call test_ks_matrix_matmul_real_NN
  Call test_ks_matrix_matmul_real_TN
  Call test_ks_matrix_matmul_real_NT
  Call test_ks_matrix_matmul_real_TT
  Call test_ks_matrix_matmul_complex_NN
  Call test_ks_matrix_matmul_complex_TN
  Call test_ks_matrix_matmul_complex_NT
  Call test_ks_matrix_matmul_complex_TT

  ! KS_array tests
  Call test_ks_array_matmul_NN
  Call test_ks_array_matmul_TN
  Call test_ks_array_matmul_NT
  Call test_ks_array_matmul_TT

  Call test_ks_split_matmul_NN
  Call test_ks_split_matmul_TN
  Call test_ks_split_matmul_NT
  Call test_ks_split_matmul_TT

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

  Subroutine test_ks_matrix_matmul_complex_NT

    Use numbers_module  , Only : wp
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_init, ks_matrix_comm_to_base, &
         ks_matrix_finalise
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( ks_matrix ) :: base
    Type( ks_matrix ) :: Am, Bm, Cm, BmT

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

    Call ks_matrix_init( n_block )
    Call ks_matrix_comm_to_base( mpi_comm_world, base )
    Call Am%create( .True., m, k, base )
    Call Am%set_by_global( 1, m, 1, k, A )
    Call Bm%create( .True., n, k, Am )
    Call Bm%set_by_global( 1, n, 1, k, B )
    BmT = .Dagger. Bm
    Cm = Am * BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, format ) 'Error in complex ks_matmul NT ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise
    
  End Subroutine test_ks_matrix_matmul_complex_NT

  Subroutine test_ks_matrix_matmul_complex_TT

    Use numbers_module  , Only : wp
    Use ks_matrix_module, Only : ks_matrix, ks_matrix_init, ks_matrix_comm_to_base, &
         ks_matrix_finalise
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

    Complex( wp ), Dimension( :, : ), Allocatable :: A, B, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Type( ks_matrix ) :: base
    Type( ks_matrix ) :: Am, Bm, Cm, AmT, BmT

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

    Call ks_matrix_init( n_block )
    Call ks_matrix_comm_to_base( mpi_comm_world, base )
    Call Am%create( .True., k, m, base )
    Call Am%set_by_global( 1, k, 1, m, A )
    Call Bm%create( .True., n, k, Am )
    Call Bm%set_by_global( 1, n, 1, k, B )
    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT * BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, format ) 'Error in complex ks_matmul TT ', Maxval( Abs( C - tmp ) ), &
            Merge( 'Passed', 'FAILED', Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise
    
  End Subroutine test_ks_matrix_matmul_complex_TT

  Subroutine test_ks_array_matmul_NN

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, B_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, B_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Bm, Cm
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:m, 1:k, 1:nk, 1:ns ) )
    Allocate( B_r( 1:k, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:k, 1:nk, 1:ns ) )
    Allocate( B_c( 1:k, 1:n, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    B_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    B_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    If( me == 0 ) Then
       Do kpoint = 1, nk
          k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do
       Do spin = 1, ns
          Do kpoint = 1, nk
             If( k_types( kpoint ) == K_POINT_REAL ) Then
                ! Real
                Call Random_number( A_r( :, :, kpoint, spin ) )
                Call Random_number( B_r( :, :, kpoint, spin ) )
                C_r( :, :, kpoint, spin ) = &
                     Matmul( A_r( :, :, kpoint, spin ), B_r( :, :, kpoint, spin ) )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:k ), rand2( 1:m, 1:k ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:k, 1:n ), rand2( 1:k, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Matmul( A_c( :, :, kpoint, spin ), B_c( :, :, kpoint, spin ) )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B_r, Size( B_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( B_c, Size( B_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( m, k, base )
    Call Bm%create( k, n, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am * Bm

    Call Cm%print_info( 'Cm - the result matrix', 9999 )
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, format ) 'Error in ks_array matmul NN ', max_diff, &
            Merge( 'Passed', 'FAILED', max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_array_matmul_NN

  Subroutine test_ks_array_matmul_TN

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, B_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, B_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Bm, Cm, AmT
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:k, 1:m, 1:nk, 1:ns ) )
    Allocate( B_r( 1:k, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:k, 1:m, 1:nk, 1:ns ) )
    Allocate( B_c( 1:k, 1:n, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    B_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    B_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    If( me == 0 ) Then
       Do kpoint = 1, nk
          k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do
       Do spin = 1, ns
          Do kpoint = 1, nk
             If( k_types( kpoint ) == K_POINT_REAL ) Then
                ! Real
                Call Random_number( A_r( :, :, kpoint, spin ) )
                Call Random_number( B_r( :, :, kpoint, spin ) )
                C_r( :, :, kpoint, spin ) = &
                     Matmul( Transpose( A_r( :, :, kpoint, spin ) ), B_r( :, :, kpoint, spin ) )
             Else
                ! Complex
                Allocate( rand1( 1:k, 1:m ), rand2( 1:k, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:k, 1:n ), rand2( 1:k, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Matmul( Conjg( Transpose( A_c( :, :, kpoint, spin ) ) ), B_c( :, :, kpoint, spin ) )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B_r, Size( B_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( B_c, Size( B_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( k, m, base )
    Call Bm%create( k, n, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    Cm = AmT * Bm

    Call Cm%print_info( 'Cm - the result matrix', 9999 )
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, format ) 'Error in ks_array matmul TN ', max_diff, &
            Merge( 'Passed', 'FAILED', max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_array_matmul_TN

  Subroutine test_ks_array_matmul_NT

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, B_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, B_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Bm, Cm, BmT
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:m, 1:k, 1:nk, 1:ns ) )
    Allocate( B_r( 1:n, 1:k, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:k, 1:nk, 1:ns ) )
    Allocate( B_c( 1:n, 1:k, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    B_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    B_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    If( me == 0 ) Then
       Do kpoint = 1, nk
          k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do
       Do spin = 1, ns
          Do kpoint = 1, nk
             If( k_types( kpoint ) == K_POINT_REAL ) Then
                ! Real
                Call Random_number( A_r( :, :, kpoint, spin ) )
                Call Random_number( B_r( :, :, kpoint, spin ) )
                C_r( :, :, kpoint, spin ) = &
                     Matmul( A_r( :, :, kpoint, spin ), Transpose( B_r( :, :, kpoint, spin ) ) )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:k ), rand2( 1:m, 1:k ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:n, 1:k ), rand2( 1:n, 1:k ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Matmul( A_c( :, :, kpoint, spin ), Conjg( Transpose( B_c( :, :, kpoint, spin ) ) ) )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B_r, Size( B_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( B_c, Size( B_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( m, k, base )
    Call Bm%create( n, k, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, n, 1, k, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, n, 1, k, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    BmT = .Dagger. Bm
    Cm = Am * BmT

    Call Cm%print_info( 'Cm - the result matrix', 9999 )
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, format ) 'Error in ks_array matmul NN ', max_diff, &
            Merge( 'Passed', 'FAILED', max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_array_matmul_NT

  Subroutine test_ks_array_matmul_TT

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, B_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, B_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Bm, Cm, AmT, BmT
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:k, 1:m, 1:nk, 1:ns ) )
    Allocate( B_r( 1:n, 1:k, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:k, 1:m, 1:nk, 1:ns ) )
    Allocate( B_c( 1:n, 1:k, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    B_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    B_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    If( me == 0 ) Then
       Do kpoint = 1, nk
          k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do
       Do spin = 1, ns
          Do kpoint = 1, nk
             If( k_types( kpoint ) == K_POINT_REAL ) Then
                ! Real
                Call Random_number( A_r( :, :, kpoint, spin ) )
                Call Random_number( B_r( :, :, kpoint, spin ) )
                C_r( :, :, kpoint, spin ) = &
                     Matmul( Transpose( A_r( :, :, kpoint, spin ) ), Transpose( B_r( :, :, kpoint, spin ) ) )
             Else
                ! Complex
                Allocate( rand1( 1:k, 1:m ), rand2( 1:k, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:n, 1:k ), rand2( 1:n, 1:k ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Matmul( Conjg( Transpose( A_c( :, :, kpoint, spin ) ) ), &
                             Conjg( Transpose( B_c( :, :, kpoint, spin ) ) ) )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B_r, Size( B_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( B_c, Size( B_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( k, m, base )
    Call Bm%create( n, k, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, n, 1, k, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, n, 1, k, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT * BmT

    Call Cm%print_info( 'Cm - the result matrix', 9999 )
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, format ) 'Error in ks_array matmul TT ', max_diff, &
            Merge( 'Passed', 'FAILED', max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_array_matmul_TT

  Subroutine test_ks_split_matmul_NN

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, B_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, B_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Bm, Cm
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:m, 1:k, 1:nk, 1:ns ) )
    Allocate( B_r( 1:k, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:k, 1:nk, 1:ns ) )
    Allocate( B_c( 1:k, 1:n, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    B_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    B_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    If( me == 0 ) Then
       k_types = K_POINT_REAL
       Do While( All( k_types == K_POINT_REAL ) .Or. All( k_types == K_POINT_COMPLEX ) )
          Do kpoint = 1, nk
             k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
             Call Random_number( rand )
             k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
          End Do
       End Do
       Do spin = 1, ns
          Do kpoint = 1, nk
             If( k_types( kpoint ) == K_POINT_REAL ) Then
                ! Real
                Call Random_number( A_r( :, :, kpoint, spin ) )
                Call Random_number( B_r( :, :, kpoint, spin ) )
                C_r( :, :, kpoint, spin ) = &
                     Matmul( A_r( :, :, kpoint, spin ), B_r( :, :, kpoint, spin ) )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:k ), rand2( 1:m, 1:k ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:k, 1:n ), rand2( 1:k, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Matmul( A_c( :, :, kpoint, spin ), B_c( :, :, kpoint, spin ) )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B_r, Size( B_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( B_c, Size( B_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, k, base )
    Call Am_base%split( 2.0_wp, Am )
    Call Am%print_info( 'Am - the split matrix', 9999 )

    Call Bm%create( k, n, Am )
    Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am * Bm
    Call Cm%print_info( 'Cm_split - the result matrix', 9999 )

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, format ) 'Error in ks_split matmul NN ', max_diff, &
            Merge( 'Passed', 'FAILED', max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_matmul_NN

  Subroutine test_ks_split_matmul_TN

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, B_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, B_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Bm, Cm, AmT
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:k, 1:m, 1:nk, 1:ns ) )
    Allocate( B_r( 1:k, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:k, 1:m, 1:nk, 1:ns ) )
    Allocate( B_c( 1:k, 1:n, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    B_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    B_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    If( me == 0 ) Then
       Do kpoint = 1, nk
          k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do
       Do spin = 1, ns
          Do kpoint = 1, nk
             If( k_types( kpoint ) == K_POINT_REAL ) Then
                ! Real
                Call Random_number( A_r( :, :, kpoint, spin ) )
                Call Random_number( B_r( :, :, kpoint, spin ) )
                C_r( :, :, kpoint, spin ) = &
                     Matmul( Transpose( A_r( :, :, kpoint, spin ) ), B_r( :, :, kpoint, spin ) )
             Else
                ! Complex
                Allocate( rand1( 1:k, 1:m ), rand2( 1:k, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:k, 1:n ), rand2( 1:k, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Matmul( Conjg( Transpose( A_c( :, :, kpoint, spin ) ) ), B_c( :, :, kpoint, spin ) )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B_r, Size( B_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( B_c, Size( B_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( k, m, base )

    Call Am_base%split( 2.0_wp, Am )
    Call Bm%create( k, n, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    Cm = AmT * Bm

    Call Cm%print_info( 'Cm - the result matrix', 9999 )
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, format ) 'Error in ks_array matmul TN ', max_diff, &
            Merge( 'Passed', 'FAILED', max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_matmul_TN

  Subroutine test_ks_split_matmul_NT

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, B_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, B_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Bm, Cm, BmT
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:m, 1:k, 1:nk, 1:ns ) )
    Allocate( B_r( 1:n, 1:k, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:k, 1:nk, 1:ns ) )
    Allocate( B_c( 1:n, 1:k, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    B_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    B_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    If( me == 0 ) Then
       Do kpoint = 1, nk
          k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do
       Do spin = 1, ns
          Do kpoint = 1, nk
             If( k_types( kpoint ) == K_POINT_REAL ) Then
                ! Real
                Call Random_number( A_r( :, :, kpoint, spin ) )
                Call Random_number( B_r( :, :, kpoint, spin ) )
                C_r( :, :, kpoint, spin ) = &
                     Matmul( A_r( :, :, kpoint, spin ), Transpose( B_r( :, :, kpoint, spin ) ) )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:k ), rand2( 1:m, 1:k ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:n, 1:k ), rand2( 1:n, 1:k ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Matmul( A_c( :, :, kpoint, spin ), Conjg( Transpose( B_c( :, :, kpoint, spin ) ) ) )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B_r, Size( B_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( B_c, Size( B_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, k, base )
    Call Am_base%split( 2.0_wp, Am )
    Call Bm%create( n, k, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, n, 1, k, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, n, 1, k, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    BmT = .Dagger. Bm
    Cm = Am * BmT

    Call Cm%print_info( 'Cm - the result matrix', 9999 )
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, format ) 'Error in ks_array matmul NN ', max_diff, &
            Merge( 'Passed', 'FAILED', max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_matmul_NT

  Subroutine test_ks_split_matmul_TT

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, B_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, B_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Bm, Cm, AmT, BmT
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:k, 1:m, 1:nk, 1:ns ) )
    Allocate( B_r( 1:n, 1:k, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:k, 1:m, 1:nk, 1:ns ) )
    Allocate( B_c( 1:n, 1:k, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    B_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    B_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    If( me == 0 ) Then
       Do kpoint = 1, nk
          k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do
       Do spin = 1, ns
          Do kpoint = 1, nk
             If( k_types( kpoint ) == K_POINT_REAL ) Then
                ! Real
                Call Random_number( A_r( :, :, kpoint, spin ) )
                Call Random_number( B_r( :, :, kpoint, spin ) )
                C_r( :, :, kpoint, spin ) = &
                     Matmul( Transpose( A_r( :, :, kpoint, spin ) ), Transpose( B_r( :, :, kpoint, spin ) ) )
             Else
                ! Complex
                Allocate( rand1( 1:k, 1:m ), rand2( 1:k, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:n, 1:k ), rand2( 1:n, 1:k ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Matmul( Conjg( Transpose( A_c( :, :, kpoint, spin ) ) ), &
                             Conjg( Transpose( B_c( :, :, kpoint, spin ) ) ) )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( B_r, Size( B_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( B_c, Size( B_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( k, m, base )
    Call Am_base%split( 2.0_wp, Am )
    Call Bm%create( n, k, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, n, 1, k, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( spin, k_points( :, kpoint ), 1, k, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( spin, k_points( :, kpoint ), 1, n, 1, k, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT * BmT

    Call Cm%print_info( 'Cm - the result matrix', 9999 )
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( spin, k_points( :, kpoint ), 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, format ) 'Error in ks_array matmul TT ', max_diff, &
            Merge( 'Passed', 'FAILED', max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_matmul_TT

End Program test_distributed_matrix
