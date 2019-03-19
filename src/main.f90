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
  Integer :: nb
  Integer :: error

  Call mpi_init( error )

  Call mpi_comm_size( mpi_comm_world, nproc, error ) 
  Call mpi_comm_rank( mpi_comm_world,    me, error )

  If( me == 0 ) Then
     Write( *, * ) 'm, n, k ?'
     Read ( *, * )  m, n, k
     Write( *, * ) 'nb ?'
     Read ( *, * )  nb
  End If
  Call mpi_bcast( m , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( n , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( k , 1, mpi_integer, 0, mpi_comm_world, error )
  Call mpi_bcast( nb, 1, mpi_integer, 0, mpi_comm_world, error )

  Call test_real_matmul_NN
  
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

    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( nb )
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
       Write( *, * ) Maxval( Abs( C - tmp ) )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_matmul_NN
  
End Program test_distributed_matrix
