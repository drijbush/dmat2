Program test_distributed_matrix

  Use numbers_module, Only : wp
  Use mpi, Only : mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, mpi_bcast, &
    mpi_comm_world, mpi_integer, mpi_abort
  use add_tests
  use diag_tests
  use multiply_tests
  use subtract_tests
  use ks_add_tests
  use ks_diag_tests
  use ks_multiply_tests
  use ks_subtract_tests
  use ks_extract_tests
  use ks_unary_tests
  use test_params

  Implicit None

  character(len=256), dimension(:), allocatable :: args
  integer, parameter :: num_test_sets = 11
  character(len=11), dimension(num_test_sets) :: test_sets = [ &
    "add        ", &
    "subtract   ", &
    "multiply   ", &
    "diag       ", &
    "ks_multiply", &
    "ks_extract ", &
    "ks_unary   ", &
    "ks_add     ", &
    "ks_subtract", &
    "ks_multiply", &
    "ks_diag    " ]
  integer :: start_arg
  integer :: num_args
  integer :: file_handle
  integer :: i

  Call mpi_init( error )

  Call mpi_comm_size( mpi_comm_world, nproc, error ) 
  Call mpi_comm_rank( mpi_comm_world,    me, error )

  num_args = command_argument_count()

  start_arg = 1
  if (num_args > 0) then
    allocate(args(num_args), stat=error)
    if (error .ne. 0) then
      write(*,*) "Failed to allocate args"
      call mpi_abort(mpi_comm_world, 1, error)
    end if
    do i = 1, num_args
      call get_command_argument(i, args(i))
    end do
  else
    allocate(args(num_test_sets), stat=error)
    if (error .ne. 0) then
      write(*,*) "Failed to allocate args"
      call mpi_abort(mpi_comm_world, 1, error)
    end if
    args = test_sets
  end if

  If( me == 0 .and. args(1) == "interactive" ) Then
    start_arg = 2
    Write( *, * ) 'm, n, k ?'
    Read ( *, * )  m, n, k
    Write( *, * ) 'n_block ?'
    Read ( *, * )  n_block
    Write( *, * ) 'ns, nk ?'
    Read ( *, * )  ns, nk
  Else If ( me == 0 ) Then
    m = m_def
    n = n_def
    k = k_def
    n_block = n_block_def
    ns = ns_def
    nk = nk_def
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
  do i = start_arg, num_args
    select case(trim(args(i)))
    case("add")
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
    case("subtract")
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
    case("multiply")
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
    case("diag")
      ! Diagonalisation tests
      If( me == 0 ) Then
        Write( *, title_format ) 'Diagonalisations'
      End If
      Call test_diag_real
      Call test_diag_complex
    case("ks_multiply")
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
        Write( *, title_format ) 'All distribution multiplies'
      End If
      Call test_ks_array_matmul_NN
      Call test_ks_array_matmul_TN
      Call test_ks_array_matmul_NT
      Call test_ks_array_matmul_TT
    case("ks_extract")
      ! Extract
      If( me == 0 ) Then
        Write( *, title_format ) 'Split distribution extract'
      End If
      Call test_ks_array_extract
      Call test_ks_array_extract_transpose
      ! Unary +/-
    case("ks_unary")
      If( me == 0 ) Then
        Write( *, title_format ) 'Split distribution unary +/-'
      End If
      Call test_ks_array_plus
      Call test_ks_array_minus
    case("ks_add")
      ! Adds
      Call test_ks_split_add_NN
      Call test_ks_split_add_TN
      Call test_ks_split_add_NT
      Call test_ks_split_add_TT
      Call test_ks_split_pre_add_diagonal
      Call test_ks_split_post_add_diagonal
    case("ks_subtract")
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
    case("ks_diag")
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
    case default
      write( *,* ) "Unknown test: ", trim(args(i))
    end select

  end do
  Call mpi_finalize( error )

Contains

  ! DISTRIBUTED_MATRIX tests


  ! multiply tests

End Program test_distributed_matrix
