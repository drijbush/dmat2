Program test_distributed_matrix

  Use mpi, Only : mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, mpi_bcast, &
       mpi_comm_world, mpi_integer
  
!!$  Use numbers_module           , Only : wp
!!$  Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
!!$       distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
!!$       distributed_matrix_set_default_blocking

  Implicit None

  Integer :: me, nproc
  Integer :: n, m, k
  Integer :: n_block
  Integer :: error

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
  Call test_complex_matmul_NN
  
  Call mpi_finalize( error )

Contains

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

    Integer :: ma, na
    Integer :: mb, nb
    Integer :: mc, nc

    ma = m
    na = k

    mb = k
    nb = n

    mc = m
    nc = n
    
    Allocate( A( 1:ma, 1:na ) )
    Allocate( B( 1:mb, 1:nb ) )
    Allocate( C( 1:mc, 1:nc ) )
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
    Call Am%create( ma, na, base )
    Call Am%set_by_global( 1, ma, 1, na, A )
    Call Bm%create( mb, nb, Am )
    Call Bm%set_by_global( 1, mb, 1, nb, B )
    Cm = Am * Bm
    Call Cm%get_by_global( 1, mc, 1, nc, tmp )
    If( me == 0 ) Then
       Write( *, * ) 'Error in real matmul NN ', Maxval( Abs( C - tmp ) )
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
       Write( *, * ) 'Error in real matmul TN ', Maxval( Abs( C - tmp ) )
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

    Integer :: ma, na
    Integer :: mb, nb
    Integer :: mc, nc

    ma = m
    na = k

    mb = n
    nb = k

    mc = m
    nc = n
    
    Allocate( A( 1:m, 1:k ) )
    Allocate( B( 1:n, 1:k ) )
    Allocate( C( 1:m, 1:n ) )
    If( me == 0 ) Then
       Call Random_number( A )
       Call Random_number( B )
       C = Matmul( A, Transpose( B ) )
       Write( *, * ) '!!!!!!', Shape( Matmul( A, Transpose( B ) ) )
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
    Write( *, * ) mc, nc, Cm%size( 1 ), Cm%size( 2 )
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, * ) 'Error in real matmul NT ', Maxval( Abs( C - tmp ) )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_NT

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
       Write( *, * ) 'Error in complex matmul NN ', Maxval( Abs( C - tmp ) )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_NN

End Program test_distributed_matrix
