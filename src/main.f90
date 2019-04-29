Program test_distributed_matrix

  Use numbers_module, Only : wp
  Use mpi, Only : mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, mpi_bcast, &
       mpi_comm_world, mpi_integer
  
  Implicit None

  Logical, Parameter :: verbose = .False.

  Integer :: me, nproc
  Integer :: n, m, k
  Integer :: n_block
  Integer :: ns, nk
  Integer :: error

  Character( Len = * ), Parameter :: error_format = '( "--> ", a, t40, g26.20, t70, a )'
  Character( Len = * ), Parameter ::   run_format = '( a, t20, i0 )'
  Character( Len = * ), Parameter :: title_format = '( t5, a )'
  Character( Len = * ), Parameter :: passed       = 'Passed      '
  Character( Len = * ), Parameter :: FAILED       = '      FAILED'
  
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

  If( me == 0 ) Then
     Write( *, * ) 'Running tests for: '
     Write( *, run_format ) ' m = ', m
     Write( *, run_format ) ' n = ', n
     Write( *, run_format ) ' k = ', k
     Write( *, run_format ) ' n_block = ', n_block
     Write( *, run_format ) ' ns = ', ns
     Write( *, run_format ) ' nk = ', nk
     Write( *, * )
  End If

  ! Distributed matrix tests
  If( me == 0 ) Then
     Write( *, * ) 'Distributed Matrix tests:'
  End If
  ! Adds
  If( me == 0 ) Then
     Write( *, title_format ) 'Additions'
  End If
  Call test_real_add_NN
  Call test_real_add_TN
  Call test_real_add_NT
  Call test_real_add_TT
  Call test_real_post_add_diagonal
  Call test_complex_add_NN
  Call test_complex_add_NT
  Call test_complex_add_TN
  Call test_complex_add_TT
  Call test_complex_post_add_diagonal
  ! Subtracts
  If( me == 0 ) Then
     Write( *, title_format ) 'Subtractions'
  End If
  Call test_real_subtract_NN
  Call test_real_subtract_TN
  Call test_real_subtract_NT
  Call test_real_subtract_TT
  Call test_complex_subtract_NN
  Call test_complex_subtract_NT
  Call test_complex_subtract_TN
  Call test_complex_subtract_TT
  ! Multiplies
  If( me == 0 ) Then
     Write( *, title_format ) 'Multiplies'
  End If
  Call test_real_pre_scale_real
  Call test_real_post_scale_real
  Call test_complex_pre_scale_real
  Call test_complex_post_scale_real
  Call test_real_matmul_NN
  Call test_real_matmul_TN
  Call test_real_matmul_NT
  Call test_real_matmul_TT
  Call test_complex_matmul_NN
  Call test_complex_matmul_TN
  Call test_complex_matmul_NT
  Call test_complex_matmul_TT
  ! Diagonalisation tests
  If( me == 0 ) Then
     Write( *, title_format ) 'Diagonalisations'
  End If
  Call test_diag_real
  Call test_diag_complex

  If( me == 0 ) Then
     Write( *, * )
  End If

  ! KS_matrix tests
  If( me == 0 ) Then
     Write( *, * ) 'KS_matrix tests'
  End If
  ! Multiplies
  If( me == 0 ) Then
     Write( *, title_format ) 'Multiplies'
  End If
  Call test_ks_matrix_matmul_real_NN
  Call test_ks_matrix_matmul_real_TN
  Call test_ks_matrix_matmul_real_NT
  Call test_ks_matrix_matmul_real_TT
  Call test_ks_matrix_matmul_complex_NN
  Call test_ks_matrix_matmul_complex_TN
  Call test_ks_matrix_matmul_complex_NT
  Call test_ks_matrix_matmul_complex_TT

  If( me == 0 ) Then
     Write( *, * )
  End If

  ! KS_array tests
  If( me == 0 ) Then
     Write( *, * ) 'KS_array tests'
  End If
  ! Extract
  If( me == 0 ) Then
     Write( *, title_format ) 'Split distribution extract'
  End If
  Call test_ks_array_extract
  Call test_ks_array_extract_transpose
  ! Unary +/-
  If( me == 0 ) Then
     Write( *, title_format ) 'Split distribution unary +/-'
  End If
  Call test_ks_array_plus
  Call test_ks_array_minus
  If( me == 0 ) Then
     Write( *, title_format ) 'Split distribution adds'
  End If
  ! Adds
  Call test_ks_split_add_NN
  Call test_ks_split_add_TN
  Call test_ks_split_add_NT
  Call test_ks_split_add_TT
  Call test_ks_split_pre_add_diagonal
  Call test_ks_split_post_add_diagonal
  ! Subtracts
  If( me == 0 ) Then
     Write( *, title_format ) 'Split distribution subtractions'
  End If
  Call test_ks_split_subtract_NN
  Call test_ks_split_subtract_TN
  Call test_ks_split_subtract_NT
  Call test_ks_split_subtract_TT
  Call test_ks_split_pre_subtract_diagonal
  Call test_ks_split_post_subtract_diagonal
  ! Multiplies
  If( me == 0 ) Then
     Write( *, title_format ) 'All distribution multiplies'
  End If
  Call test_ks_array_matmul_NN
  Call test_ks_array_matmul_TN
  Call test_ks_array_matmul_NT
  Call test_ks_array_matmul_TT
  If( me == 0 ) Then
     Write( *, title_format ) 'Split distribution multiplies'
  End If
  Call test_ks_split_real_pre_scale
  Call test_ks_split_real_post_scale
  Call test_ks_split_matmul_NN
  Call test_ks_split_matmul_TN
  Call test_ks_split_matmul_NT
  Call test_ks_split_matmul_TT

  Call test_ks_rejoin_matmul_NN  

  If( me == 0 ) Then
     Write( *, title_format ) 'Split distribution diags'
  End If
  Call test_ks_array_diag
  If( me == 0 ) Then
     Write( *, title_format ) 'Split distribution Choleskis'
  End If
  Call test_ks_array_choleski
  If( me == 0 ) Then
     Write( *, title_format ) 'Split distribution Triangular inverts'
  End If
  Call test_ks_array_tr_inv
  Call test_ks_array_tr_inv_with_iterator

  Call mpi_finalize( error )

Contains

  ! DISTRIBUTED_MATRIX tests

  ! Add tests
  Subroutine test_real_add_NN
    
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
       C = A + B
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
    Cm = Am + Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in real add NN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_add_NN

  Subroutine test_real_add_TN
    
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
       C = Transpose( A ) + B
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
    Cm = AmT + Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in real add TN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_add_TN

  Subroutine test_real_add_NT
    
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
       C = A + Transpose( B )
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
    Cm = Am + BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in real add NT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_add_NT

  Subroutine test_real_add_TT
    
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
       C = Transpose( A ) + Transpose( B )
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
    Cm = AmT + BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in real add TT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_add_TT

  Subroutine test_real_post_add_diagonal
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, real_distributed_matrix, &
         distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
         distributed_matrix_set_default_blocking

    Implicit None
    
    Real( wp ), Dimension( :, : ), Allocatable :: A, C, tmp

    Real( wp ), Dimension( : ), Allocatable :: d

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Cm

    Integer :: i
    
    Allocate( A( 1:m, 1:m ) )
    Allocate( d( 1:m      ) )
    Allocate( C( 1:m, 1:m ) )
    If( me == 0 ) Then
       Call Random_number( A )
       Call Random_number( d )
       C = A
       Do i = 1, m
          C( i, i ) = A( i, i ) + d( i )
       End Do
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( d, Size( d ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_precision, 0, mpi_comm_world, error )
    
    Allocate( tmp, Mold = C )
    
    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( real_distributed_matrix :: Am )
    Call Am%create( m, m, base )
    Call Am%set_by_global( 1, m, 1, m, A )
    Cm =  Am + d
    Call Cm%get_by_global( 1, m, 1, m, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in real post_add diagonal ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_real_post_add_diagonal

  Subroutine test_complex_post_add_diagonal
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_double_complex

    Use numbers_module           , Only : wp
    Use distributed_matrix_module, Only : distributed_matrix, complex_distributed_matrix, real_distributed_matrix, &
         distributed_matrix_init, distributed_matrix_comm_to_base, distributed_matrix_finalise, &
         distributed_matrix_set_default_blocking

    Implicit None
    
    Complex( wp ), Dimension( :, : ), Allocatable :: A, C, tmp

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ), Dimension( : ), Allocatable :: d

    Type( real_distributed_matrix ) :: base
    Class( distributed_matrix ), Allocatable :: Am, Cm

    Integer :: i
    
    Allocate( A( 1:m, 1:m ) )
    Allocate( d( 1:m      ) )
    Allocate( C( 1:m, 1:m ) )
    If( me == 0 ) Then
       Allocate( rand1( 1:m, 1:m ), rand2( 1:m, 1:m ) )
       Call Random_number( rand1 ); Call Random_number( rand2 ) 
       A = Cmplx( rand1, rand2, wp )
       Call Random_number( d )
       C = A
       Do i = 1, m
          C( i, i ) = A( i, i ) + d( i )
       End Do
    End If
    Call mpi_bcast( A, Size( A ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( d, Size( d ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C, Size( C ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp, Mold = C )
    
    Call distributed_matrix_init
    Call distributed_matrix_set_default_blocking( n_block )
    Call distributed_matrix_comm_to_base( mpi_comm_world, base )
    Allocate( complex_distributed_matrix :: Am )
    Call Am%create( m, m, base )
    Call Am%set_by_global( 1, m, 1, m, A )
    Cm =  Am + d
    Call Cm%get_by_global( 1, m, 1, m, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in real post_add diagonal ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_post_add_diagonal

  Subroutine test_complex_add_NN
    
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
       C = A + B
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
    Cm = Am + Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in complex add NN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_add_NN

  Subroutine test_complex_add_TN
    
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
       C = Transpose( Conjg( A ) ) + B
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
    Cm = AmT + Bm
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in complex add TN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_add_TN

  Subroutine test_complex_add_NT
    
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
       C = A+ Transpose( Conjg( B ) )
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
    Cm = Am + BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in complex add NT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_add_NT

  Subroutine test_complex_add_TT
    
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
       C = Transpose( Conjg( A ) ) + Transpose( Conjg( B ) )
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
    Cm = AmT + BmT
    Call Cm%get_by_global( 1, m, 1, n, tmp )
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in complex add TT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_add_TT

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

  ! Multiply tests
  Subroutine test_real_pre_scale_real
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

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
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

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
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

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
    
    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

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
       Write( *, error_format ) 'Error in real matmul NN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in real matmul TN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in real matmul NT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in real matmul TT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in complex matmul NN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in complex matmul TN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in complex matmul NT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in complex matmul TT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call distributed_matrix_finalise

  End Subroutine test_complex_matmul_TT

  ! Diag tests
  Subroutine test_diag_real()

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_precision

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

    Use mpi, Only : mpi_bcast, mpi_comm_world, mpi_double_complex

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
       Write( *, error_format ) 'Error in real ks_matmul NN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in real ks_matmul TN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in real ks_matmul NT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in real ks_matmul TT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in complex ks_matmul NN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in complex ks_matmul TN ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in complex ks_matmul NT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
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
       Write( *, error_format ) 'Error in complex ks_matmul TT ', Maxval( Abs( C - tmp ) ), &
            Merge( passed, FAILED, Maxval( Abs( C - tmp ) ) < tol )
    End If
    Call ks_matrix_finalise
    
  End Subroutine test_ks_matrix_matmul_complex_TT

  ! KS_array tests

  ! Extract
  
  Subroutine test_ks_array_extract_transpose()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, ks_eval_storage
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: B
    
    Type( ks_array ) :: base_k

    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s
    Integer :: n
    Integer :: r1, r2, c1, c2
    Integer :: error

    r1 = 1
    r2 = n_block + 1
    c1 = 1
    c2 = n_block + 2

    n = Max( m, r1, r2, c1, c2 )

    Allocate( A_r( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:n ) )
    Allocate( tmp_c( 1:n, 1:n ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )
    
    If( me == 0 ) Then

       Call Random_number( A_r )

       Do s = 1, ns
          Do k = 1, nk
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = tmp_r
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
          End Do
       End Do
       
       Do k = 1, nk
          k_points( :, k ) = [ k - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do

    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( n, n, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
          Else
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
          End If
       End Do
    End Do

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If

    Deallocate( tmp_r )
    Deallocate( tmp_c )
    Allocate( tmp_r( r1:r2, c1:c2 ) )
    Allocate( tmp_c( r1:r2, c1:c2 ) )

    A_split = .Dagger. A_split
    B = A_split%extract( r1, r2, c1, c2 )
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call B%get_by_global( k_points( :, k ), s, 1, r2 - r1 + 1, 1, c2 - c1 + 1, tmp_r )
             tmp_r = Abs( tmp_r - Transpose( A_r( c1:c2, r1:r2, k, s ) ) )
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, r2 - r1 + 1, 1, c2 - c1 + 1, tmp_c )
             tmp_r = Abs( tmp_c - Transpose( Conjg( A_c( c1:c2, r1:r2, k, s ) ) ) )
          End If
          max_diff = Max( max_diff, Maxval( tmp_r ) )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split extract transpose ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If
    
    Call ks_array_finalise

  End Subroutine test_ks_array_extract_transpose

  Subroutine test_ks_array_extract()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, ks_eval_storage
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: B
    
    Type( ks_array ) :: base_k

    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s
    Integer :: n
    Integer :: r1, r2, c1, c2
    Integer :: error

    r1 = 2
    r2 = n_block + 1
    c1 = 3
    c2 = n_block + 2

    n = Max( m, r1, r2, c1, c2 )

    Allocate( A_r( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:n ) )
    Allocate( tmp_c( 1:n, 1:n ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )
    
    If( me == 0 ) Then

       Call Random_number( A_r )

       Do s = 1, ns
          Do k = 1, nk
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = tmp_r
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
          End Do
       End Do
       
       Do k = 1, nk
          k_points( :, k ) = [ k - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do

    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( n, n, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
          Else
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
          End If
       End Do
    End Do

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If

    Deallocate( tmp_r )
    Deallocate( tmp_c )
    Allocate( tmp_r( r1:r2, c1:c2 ) )
    Allocate( tmp_c( r1:r2, c1:c2 ) )
    
    B = A_split%extract( r1, r2, c1, c2 )
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call B%get_by_global( k_points( :, k ), s, 1, r2 - r1 + 1, 1, c2 - c1 + 1, tmp_r )
             tmp_r = Abs( tmp_r - A_r( r1:r2, c1:c2, k, s ) )
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, r2 - r1 + 1, 1, c2 - c1 + 1, tmp_c )
             tmp_r = Abs( tmp_c - A_c( r1:r2, c1:c2, k, s ) )
          End If
          max_diff = Max( max_diff, Maxval( tmp_r ) )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split extract ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If
    
    Call ks_array_finalise

  End Subroutine test_ks_array_extract

  ! Unary plus/minus
  
  Subroutine test_ks_array_plus()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, ks_eval_storage
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: B
    
    Type( ks_array ) :: base_k

    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s
    Integer :: n
    Integer :: error

    n = m

    Allocate( A_r( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:n ) )
    Allocate( tmp_c( 1:n, 1:n ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )
    
    If( me == 0 ) Then

       Call Random_number( A_r )
       Do s = 1, ns
          Do k = 1, nk
             ! Make sure matrix is positive definite
             A_r( :, :, k, s ) = Matmul( A_r( :, :, k, s ), Transpose( A_r( :, :, k, s  ) ) )
          End Do
       End Do

       Do s = 1, ns
          Do k = 1, nk
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = tmp_r
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
             ! Make sure matrix is positive definite
             A_c( :, :, k, s ) = Matmul( A_c( :, :, k, s ), Conjg( Transpose( A_c( :, :, k, s  ) ) ) )
          End Do
       End Do
       
       Do k = 1, nk
          k_points( :, k ) = [ k - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do

    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( n, n, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
          Else
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
          End If
       End Do
    End Do

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If

    B = + A_split
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
       Do k = 1, nk
          tmp_r = 0.0_wp
          If( k_types( k ) == K_POINT_REAL ) Then
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_r )
             tmp_r = Abs( tmp_r - A_r( :, :, k, s ) )
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_c )
             tmp_r = Abs( tmp_c - A_c( :, :, k, s ) )
          End If
          max_diff = Max( max_diff, Maxval( tmp_r ) )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split plus ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If
    
    Call ks_array_finalise

  End Subroutine test_ks_array_plus

  Subroutine test_ks_array_minus()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, ks_eval_storage
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: B
    
    Type( ks_array ) :: base_k

    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s
    Integer :: n
    Integer :: error

    n = m

    Allocate( A_r( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:n ) )
    Allocate( tmp_c( 1:n, 1:n ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )
    
    If( me == 0 ) Then

       Call Random_number( A_r )
       Do s = 1, ns
          Do k = 1, nk
             ! Make sure matrix is positive definite
             A_r( :, :, k, s ) = Matmul( A_r( :, :, k, s ), Transpose( A_r( :, :, k, s  ) ) )
          End Do
       End Do

       Do s = 1, ns
          Do k = 1, nk
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = tmp_r
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
             ! Make sure matrix is positive definite
             A_c( :, :, k, s ) = Matmul( A_c( :, :, k, s ), Conjg( Transpose( A_c( :, :, k, s  ) ) ) )
          End Do
       End Do
       
       Do k = 1, nk
          k_points( :, k ) = [ k - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do

    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( n, n, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
          Else
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
          End If
       End Do
    End Do

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If

    B = - A_split
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
       Do k = 1, nk
          tmp_r = 0.0_wp
          If( k_types( k ) == K_POINT_REAL ) Then
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_r )
             tmp_r = Abs( tmp_r + A_r( :, :, k, s ) )
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_c )
             tmp_r = Abs( tmp_c + A_c( :, :, k, s ) )
          End If
          max_diff = Max( max_diff, Maxval( tmp_r ) )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split minus ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If
    
    Call ks_array_finalise

  End Subroutine test_ks_array_minus

  ! Add tests
  Subroutine test_ks_split_add_NN

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

    Allocate( A_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( B_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( B_c( 1:m, 1:n, 1:nk, 1:ns ) )
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
                C_r( :, :, kpoint, spin ) = A_r( :, :, kpoint, spin ) + B_r( :, :, kpoint, spin )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = A_c( :, :, kpoint, spin ) + B_c( :, :, kpoint, spin )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, n, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( m, n, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am + Bm
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split add NN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_add_NN

  Subroutine test_ks_split_add_TN

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

    Allocate( A_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( B_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( B_c( 1:m, 1:n, 1:nk, 1:ns ) )
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
                     Transpose( A_r( :, :, kpoint, spin ) ) + B_r( :, :, kpoint, spin )
             Else
                ! Complex
                Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Transpose( Conjg( A_c( :, :, kpoint, spin ) ) ) + B_c( :, :, kpoint, spin )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( n, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( m, n, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    Cm = AmT + Bm
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split add TN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_add_TN

  Subroutine test_ks_split_add_NT

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

    Allocate( A_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( B_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( B_c( 1:n, 1:m, 1:nk, 1:ns ) )
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
                C_r( :, :, kpoint, spin ) = A_r( :, :, kpoint, spin ) + &
                     Transpose( B_r( :, :, kpoint, spin ) )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = A_c( :, :, kpoint, spin ) + &
                     Transpose( Conjg( B_c( :, :, kpoint, spin ) ) )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, n, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( n, m, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    BmT = .Dagger. Bm
    Cm = Am + BmT
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split add NT ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_add_NT

  Subroutine test_ks_split_add_TT

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

    Allocate( A_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( B_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( B_c( 1:n, 1:m, 1:nk, 1:ns ) )
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
                     Transpose( A_r( :, :, kpoint, spin ) ) + Transpose( B_r( :, :, kpoint, spin ) )
             Else
                ! Complex
                Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Transpose( Conjg( A_c( :, :, kpoint, spin ) ) ) + &
                     Transpose( Conjg( B_c( :, :, kpoint, spin ) ) )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( n, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( n, m, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT + BmT
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split add TT ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_add_TT

  Subroutine test_ks_split_post_add_diagonal

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ), Dimension( : ), Allocatable :: d
    
    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Cm
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin
    Integer :: i

    Allocate( A_r( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:m, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    Allocate( d( 1:m ) )
    If( me == 0 ) Then
       Call Random_number( d )
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
                C_r( :, :, kpoint, spin ) = A_r( :, :, kpoint, spin )
                Do i = 1, m
                   C_r( i, i, kpoint, spin ) = A_r( i, i, kpoint, spin ) + d( i )
                End Do
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:m ), rand2( 1:m, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = A_c( :, :, kpoint, spin )
                Do i = 1, m
                   C_c( i, i, kpoint, spin ) = A_c( i, i, kpoint, spin ) + d( i )
                End Do
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Call mpi_bcast( d, Size( d ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp_r( 1:m, 1:m ) )
    Allocate( tmp_c( 1:m, 1:m ) )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, A_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, A_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am + d
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split post-add ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_post_add_diagonal

  Subroutine test_ks_split_pre_add_diagonal

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ), Dimension( : ), Allocatable :: d
    
    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Cm
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin
    Integer :: i

    Allocate( A_r( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:m, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    Allocate( d( 1:m ) )
    If( me == 0 ) Then
       Call Random_number( d )
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
                C_r( :, :, kpoint, spin ) = A_r( :, :, kpoint, spin )
                Do i = 1, m
                   C_r( i, i, kpoint, spin ) = A_r( i, i, kpoint, spin ) + d( i )
                End Do
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:m ), rand2( 1:m, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = A_c( :, :, kpoint, spin )
                Do i = 1, m
                   C_c( i, i, kpoint, spin ) = A_c( i, i, kpoint, spin ) + d( i )
                End Do
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Call mpi_bcast( d, Size( d ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp_r( 1:m, 1:m ) )
    Allocate( tmp_c( 1:m, 1:m ) )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, A_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, A_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = d + Am
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split pre-add ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_pre_add_diagonal

  ! Subtract tests
  Subroutine test_ks_split_subtract_NN

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

    Allocate( A_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( B_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( B_c( 1:m, 1:n, 1:nk, 1:ns ) )
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
                C_r( :, :, kpoint, spin ) = A_r( :, :, kpoint, spin ) - B_r( :, :, kpoint, spin )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = A_c( :, :, kpoint, spin ) - B_c( :, :, kpoint, spin )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, n, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( m, n, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am - Bm
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split subtract NN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_subtract_NN

  Subroutine test_ks_split_subtract_TN

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

    Allocate( A_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( B_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( B_c( 1:m, 1:n, 1:nk, 1:ns ) )
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
                     Transpose( A_r( :, :, kpoint, spin ) ) - B_r( :, :, kpoint, spin )
             Else
                ! Complex
                Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Transpose( Conjg( A_c( :, :, kpoint, spin ) ) ) - B_c( :, :, kpoint, spin )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( n, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( m, n, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    Cm = AmT - Bm
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split subtract TN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_subtract_TN

  Subroutine test_ks_split_subtract_NT

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

    Allocate( A_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( B_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( B_c( 1:n, 1:m, 1:nk, 1:ns ) )
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
                C_r( :, :, kpoint, spin ) = A_r( :, :, kpoint, spin ) - &
                     Transpose( B_r( :, :, kpoint, spin ) )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = A_c( :, :, kpoint, spin ) - &
                     Transpose( Conjg( B_c( :, :, kpoint, spin ) ) )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, n, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( n, m, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    BmT = .Dagger. Bm
    Cm = Am - BmT
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split subtract NT ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_subtract_NT

  Subroutine test_ks_split_subtract_TT

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

    Allocate( A_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( B_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( B_c( 1:n, 1:m, 1:nk, 1:ns ) )
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
                     Transpose( A_r( :, :, kpoint, spin ) ) - Transpose( B_r( :, :, kpoint, spin ) )
             Else
                ! Complex
                Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                Allocate( rand1( 1:n, 1:m ), rand2( 1:n, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                B_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = &
                     Transpose( Conjg( A_c( :, :, kpoint, spin ) ) ) - &
                     Transpose( Conjg( B_c( :, :, kpoint, spin ) ) )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( n, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( n, m, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, m, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT - BmT
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split subtract TT ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_subtract_TT

  Subroutine test_ks_split_post_subtract_diagonal

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ), Dimension( : ), Allocatable :: d
    
    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Cm
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin
    Integer :: i

    Allocate( A_r( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:m, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    Allocate( d( 1:m ) )
    If( me == 0 ) Then
       Call Random_number( d )
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
                C_r( :, :, kpoint, spin ) = A_r( :, :, kpoint, spin )
                Do i = 1, m
                   C_r( i, i, kpoint, spin ) = A_r( i, i, kpoint, spin ) - d( i )
                End Do
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:m ), rand2( 1:m, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = A_c( :, :, kpoint, spin )
                Do i = 1, m
                   C_c( i, i, kpoint, spin ) = A_c( i, i, kpoint, spin ) - d( i )
                End Do
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Call mpi_bcast( d, Size( d ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp_r( 1:m, 1:m ) )
    Allocate( tmp_c( 1:m, 1:m ) )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, A_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, A_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am - d
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split post-subtract ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_post_subtract_diagonal

  Subroutine test_ks_split_pre_subtract_diagonal

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ), Dimension( : ), Allocatable :: d
    
    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Cm
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin
    Integer :: i

    Allocate( A_r( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:m, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:m, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
    C_c = Huge( 1.0_wp )
    Allocate( d( 1:m ) )
    If( me == 0 ) Then
       Call Random_number( d )
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
                C_r( :, :, kpoint, spin ) = - A_r( :, :, kpoint, spin )
                Do i = 1, m
                   C_r( i, i, kpoint, spin ) = d( i ) - A_r( i, i, kpoint, spin )
                End Do
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:m ), rand2( 1:m, 1:m ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = - A_c( :, :, kpoint, spin )
                Do i = 1, m
                   C_c( i, i, kpoint, spin ) = d( i ) - A_c( i, i, kpoint, spin )
                End Do
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Call mpi_bcast( d, Size( d ), mpi_double_precision, 0, mpi_comm_world, error )

    Allocate( tmp_r( 1:m, 1:m ) )
    Allocate( tmp_c( 1:m, 1:m ) )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, A_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, A_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = d - Am
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, m, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split pre-subtract ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_pre_subtract_diagonal

  ! multiply tests
  
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( m, k, base )
    Call Bm%create( k, n, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am * Bm

    If( verbose ) Then
       Call Cm%print_info( 'Cm - the result matrix', 9999 )
    End If
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_array matmul NN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( k, m, base )
    Call Bm%create( k, n, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    Cm = AmT * Bm

    If( verbose ) Then
       Call Cm%print_info( 'Cm - the result matrix', 9999 )
    End If
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_array matmul TN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( m, k, base )
    Call Bm%create( n, k, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, k, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, k, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    BmT = .Dagger. Bm
    Cm = Am * BmT

    If( verbose ) Then
       Call Cm%print_info( 'Cm - the result matrix', 9999 )
    End If
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_array matmul NN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( k, m, base )
    Call Bm%create( n, k, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, k, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, k, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT * BmT

    If( verbose ) Then
       Call Cm%print_info( 'Cm - the result matrix', 9999 )
    End If
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_array matmul TT ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_array_matmul_TT

  Subroutine test_ks_rejoin_matmul_NN

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
    Type( ks_array ) :: Am
    Type( ks_array ) :: Ams, Bms, Cms
    Type( ks_array ) :: Cmj
    
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am%create( m, k, base )
    Call Am%split_ks( 2.0_wp, Ams )
    Call Bms%create( k, n, Ams   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Ams%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bms%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Ams%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bms%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cms = Ams * Bms

    Call Cms%join_ks( Cmj )

    If( verbose ) Then
       Call Cmj%print_info( 'Cmj - the result matrix', 9999 )
    End If
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cmj%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cmj%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_rejoin matmul NN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_rejoin_matmul_NN

  Subroutine test_ks_split_real_pre_scale

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Cm
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
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
                C_r( :, :, kpoint, spin ) = 3.0_wp * A_r( :, :, kpoint, spin )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = 3.0_wp * A_c( :, :, kpoint, spin )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, n, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = 3.0_wp * Am
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split pre-scale ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_real_pre_scale

  Subroutine test_ks_split_real_post_scale

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
         

    Real   ( wp ), Dimension( :, :, :, : ), Allocatable :: A_r, C_r
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c, C_c

    Real   ( wp ), Dimension( :, : ), Allocatable :: tmp_r
    Complex( wp ), Dimension( :, : ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, : ), Allocatable :: rand1, rand2

    Real( wp ) :: max_diff, this_diff
    Real( wp ) :: rand

    Type( ks_array ) :: base
    Type( ks_array ) :: Am, Cm
    Type( ks_array ) :: Am_base
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: kpoint, spin

    Allocate( A_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( C_r( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:m, 1:n, 1:nk, 1:ns ) )
    Allocate( C_c( 1:m, 1:n, 1:nk, 1:ns ) )
    A_r = Huge( 1.0_wp )
    C_r = Huge( 1.0_wp )
    A_c = Huge( 1.0_wp )
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
                C_r( :, :, kpoint, spin ) = 3.0_wp * A_r( :, :, kpoint, spin )
             Else
                ! Complex
                Allocate( rand1( 1:m, 1:n ), rand2( 1:m, 1:n ) )
                Call Random_number( rand1 ); Call Random_number( rand2 ) 
                A_c( :, :, kpoint, spin ) = Cmplx( rand1, rand2, wp )
                Deallocate( rand1, rand2 )
                C_c( :, :, kpoint, spin ) = 3.0_wp * A_c( :, :, kpoint, spin )
             End If
          End Do
       End Do
    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )
    Call mpi_bcast( C_r, Size( C_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    Call mpi_bcast( C_c, Size( C_c ), mpi_double_complex  , 0, mpi_comm_world, error )
    
    Allocate( tmp_r( 1:m, 1:n ) )
    Allocate( tmp_c( 1:m, 1:n ) )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, n, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, A_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am * 3.0_wp
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split post-scale ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_real_post_scale

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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, k, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    If( verbose ) Then
       Call Am%print_info( 'Am - the split matrix', 9999 )
    End If

    Call Bm%create( k, n, Am )
    If( verbose ) Then
       Call Bm%print_info( 'Bm - the derived matrix', 9999 )
    End If
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    Cm = Am * Bm
    If( verbose ) Then
       Call Cm%print_info( 'Cm_split - the result matrix', 9999 )
    End If

    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split matmul NN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( k, m, base )

    Call Am_base%split_ks( 2.0_wp, Am )
    Call Bm%create( k, n, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, n, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    Cm = AmT * Bm

    If( verbose ) Then
       Call Cm%print_info( 'Cm - the result matrix', 9999 )
    End If
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split matmul TN ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( m, k, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    Call Bm%create( n, k, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, k, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, m, 1, k, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, k, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    BmT = .Dagger. Bm
    Cm = Am * BmT

    If( verbose ) Then
       Call Cm%print_info( 'Cm - the result matrix', 9999 )
    End If
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split matmul NT ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
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

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    Call Am_base%create( k, m, base )
    Call Am_base%split_ks( 2.0_wp, Am )
    Call Bm%create( n, k, Am   )
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, m, A_r( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, k, B_r( :, :, kpoint, spin ) )
          Else
             Call Am%set_by_global( k_points( :, kpoint ), spin, 1, k, 1, m, A_c( :, :, kpoint, spin ) )
             Call Bm%set_by_global( k_points( :, kpoint ), spin, 1, n, 1, k, B_c( :, :, kpoint, spin ) )
          End If
       End Do
    End Do

    AmT = .Dagger. Am
    BmT = .Dagger. Bm
    Cm = AmT * BmT

    If( verbose ) Then
       Call Cm%print_info( 'Cm - the result matrix', 9999 )
    End If
    max_diff = -1.0_wp
    Do spin = 1, ns
       Do kpoint = 1, nk
          If( k_types( kpoint ) == K_POINT_REAL ) Then
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_r )
             this_diff = Maxval( Abs( C_r( :, :, kpoint, spin ) - tmp_r ) )
          Else
             Call Cm%get_by_global( k_points( :, kpoint ), spin, 1, m, 1, n, tmp_c ) 
             this_diff = Maxval( Abs( C_c( :, :, kpoint, spin ) - tmp_c ) )
          End If
          max_diff = Max( this_diff, max_diff )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split matmul TT ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_split_matmul_TT

  Subroutine test_ks_array_diag()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, ks_eval_storage
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: Q
    Type( ks_array ) :: QT
    Type( ks_array ) :: B
    
    Type( ks_eval_storage ), Dimension( 1:nk * ns ) :: E
    
    Type( ks_array ) :: base_k

    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s
    Integer :: i
    Integer :: n
    Integer :: error

    n = m

    Allocate( A_r( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:n ) )
    Allocate( tmp_c( 1:n, 1:n ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )
    
    If( me == 0 ) Then

       Call Random_number( A_r )
       Do s = 1, ns
          Do k = 1, nk
             A_r( :, :, k, s ) = A_r( :, :, k, s ) + Transpose( A_r( :, :, k, s  ) )
          End Do
       End Do

       Do s = 1, ns
          Do k = 1, nk
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = tmp_r
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Conjg( Transpose( A_c( :, :, k, s  ) ) )
          End Do
       End Do
       
       Do k = 1, nk
          k_points( :, k ) = [ k - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do

    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( n, n, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
          Else
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
          End If
       End Do
    End Do

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If

    Call A_split%diag( Q, E )
    QT = .Dagger. Q
    B = QT * A_split * Q
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
       Do k = 1, nk
          tmp_r = 0.0_wp
          If( k_types( k ) == K_POINT_REAL ) Then
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_r )
             Do i = 1, n
                tmp_r( i, i ) = Abs( tmp_r( i, i ) - E( k + ( s - 1 ) * nk )%evals( i ) )
             End Do
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_c )
             Do i = 1, n
                tmp_r( i, i ) = Abs( tmp_c( i, i ) - E( k + ( s - 1 ) * nk )%evals( i ) )
             End Do
          End If
          max_diff = Max( max_diff, Maxval( tmp_r ) )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split diag ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If
    
    Call ks_array_finalise

  End Subroutine test_ks_array_diag

  Subroutine test_ks_array_choleski()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, ks_eval_storage
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: L
    Type( ks_array ) :: LT
    Type( ks_array ) :: B
    
    Type( ks_array ) :: base_k

    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s
    Integer :: n
    Integer :: error

    n = m

    Allocate( A_r( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:n ) )
    Allocate( tmp_c( 1:n, 1:n ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )
    
    If( me == 0 ) Then

       Call Random_number( A_r )
       Do s = 1, ns
          Do k = 1, nk
             ! Make sure matrix is positive definite
             A_r( :, :, k, s ) = Matmul( A_r( :, :, k, s ), Transpose( A_r( :, :, k, s  ) ) )
          End Do
       End Do

       Do s = 1, ns
          Do k = 1, nk
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = tmp_r
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
             ! Make sure matrix is positive definite
             A_c( :, :, k, s ) = Matmul( A_c( :, :, k, s ), Conjg( Transpose( A_c( :, :, k, s  ) ) ) )
          End Do
       End Do
       
       Do k = 1, nk
          k_points( :, k ) = [ k - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do

    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( n, n, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
          Else
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
          End If
       End Do
    End Do

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If

    L = .Choleski. A_split
    LT = .Dagger. L
    B = L * LT - A_split
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
       Do k = 1, nk
          tmp_r = 0.0_wp
          If( k_types( k ) == K_POINT_REAL ) Then
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_r )
             tmp_r = Abs( tmp_r )
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_c )
             tmp_r = Abs( tmp_c )
          End If
          max_diff = Max( max_diff, Maxval( tmp_r ) )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split choleski ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If
    
    Call ks_array_finalise

  End Subroutine test_ks_array_choleski

  Subroutine test_ks_array_tr_inv()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, ks_eval_storage
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: L
    Type( ks_array ) :: L_inv
    Type( ks_array ) :: B
    
    Type( ks_array ) :: base_k

    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s
    Integer :: i
    Integer :: n
    Integer :: error

    n = m

    Allocate( A_r( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:n ) )
    Allocate( tmp_c( 1:n, 1:n ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )
    
    If( me == 0 ) Then

       Call Random_number( A_r )
       Do s = 1, ns
          Do k = 1, nk
             ! Make sure matrix is positive definite
             A_r( :, :, k, s ) = Matmul( A_r( :, :, k, s ), Transpose( A_r( :, :, k, s  ) ) )
          End Do
       End Do

       Do s = 1, ns
          Do k = 1, nk
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = tmp_r
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
             ! Make sure matrix is positive definite
             A_c( :, :, k, s ) = Matmul( A_c( :, :, k, s ), Conjg( Transpose( A_c( :, :, k, s  ) ) ) )
          End Do
       End Do
       
       Do k = 1, nk
          k_points( :, k ) = [ k - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do

    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( n, n, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
    Do s = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
          Else
             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
          End If
       End Do
    End Do

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If

    ! Generate a traingular array by Choleski Decomposition
    L = .Choleski. A_split
    L_inv = .TrInv. L
    B = L * L_inv
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
       Do k = 1, nk
          tmp_r = 0.0_wp
          If( k_types( k ) == K_POINT_REAL ) Then
             ! B should be the unit matrix
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_r )
             tmp_r = Abs( tmp_r )
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_c )
             tmp_r = Abs( tmp_c )
          End If
          Do i = 1, n
             tmp_r( i, i ) = tmp_r( i, i ) - 1.0_wp
          End Do
          max_diff = Max( max_diff, Maxval( tmp_r ) )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split tr_inv ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If
    
    Call ks_array_finalise

  End Subroutine test_ks_array_tr_inv

  Subroutine test_ks_array_tr_inv_with_iterator()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_point_info, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, K_POINT_NOT_EXIST, ks_eval_storage
    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: L
    Type( ks_array ) :: L_inv
    Type( ks_array ) :: B
    
    Type( ks_array ) :: base_k

    Type( ks_point_info ), Allocatable :: info
    
    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff
    
    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s
    Integer :: i
    Integer :: n
    Integer :: error

    n = m

    Allocate( A_r( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:n, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:n ) )
    Allocate( tmp_c( 1:n, 1:n ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )
    
    If( me == 0 ) Then

       Call Random_number( A_r )
       Do s = 1, ns
          Do k = 1, nk
             ! Make sure matrix is positive definite
             A_r( :, :, k, s ) = Matmul( A_r( :, :, k, s ), Transpose( A_r( :, :, k, s  ) ) )
          End Do
       End Do

       Do s = 1, ns
          Do k = 1, nk
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = tmp_r
             Call Random_number( tmp_r )
             A_c( :, :, k, s ) = A_c( :, :, k, s ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
             ! Make sure matrix is positive definite
             A_c( :, :, k, s ) = Matmul( A_c( :, :, k, s ), Conjg( Transpose( A_c( :, :, k, s  ) ) ) )
          End Do
       End Do
       
       Do k = 1, nk
          k_points( :, k ) = [ k - 1, 0, 0 ]
          Call Random_number( rand )
          k_types( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
       End Do

    End If

    Call mpi_bcast( k_points, Size( k_points ), mpi_integer, 0, mpi_comm_world, error )
    Call mpi_bcast( k_types , Size( k_types  ), mpi_integer, 0, mpi_comm_world, error )
    
    Call mpi_bcast( A_r, Size( A_r ), mpi_double_precision, 0, mpi_comm_world, error )

    Call mpi_bcast( A_c, Size( A_c ), mpi_double_complex  , 0, mpi_comm_world, error )

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( n, n, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
!!$    Do s = 1, ns
!!$       Do k = 1, nk
!!$          If( k_types( k ) == K_POINT_REAL ) Then
!!$             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
!!$          Else
!!$             Call A%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
!!$          End If
!!$       End Do
!!$    End Do
    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If

    Call A_split%iterator_init()
    Do
       info = A_split%iterator_next()
!!$       Write( *, * ) info%k_type
       If( info%k_type == K_POINT_NOT_EXIST ) Exit
       Outer: Do s = 1, ns
          Do k = 1, nk
             If( All( info%k_indices == k_points( :, k ) .And. info%spin == s ) ) Then
                If( info%k_type == K_POINT_REAL ) Then
                   Call A_split%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_r( :, :, k, s ) )
                Else
                   Call A_split%set_by_global( k_points( :, k ), s, 1, n, 1, n, A_c( :, :, k, s ) )
                End If
                Exit Outer
             End If
          End Do
       End Do Outer
    End Do
    Call A_split%iterator_reset()

    ! Generate a traingular array by Choleski Decomposition
    L = .Choleski. A_split
    L_inv = .TrInv. L
    B = L * L_inv
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
       Do k = 1, nk
          tmp_r = 0.0_wp
          If( k_types( k ) == K_POINT_REAL ) Then
             ! B should be the unit matrix
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_r )
             tmp_r = Abs( tmp_r )
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_c )
             tmp_r = Abs( tmp_c )
          End If
          Do i = 1, n
             tmp_r( i, i ) = tmp_r( i, i ) - 1.0_wp
          End Do
          max_diff = Max( max_diff, Maxval( tmp_r ) )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split tr_inv iterator ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If
    
    Call ks_array_finalise

  End Subroutine test_ks_array_tr_inv_with_iterator

End Program test_distributed_matrix
