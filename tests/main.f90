Program test_distributed_matrix

  Use mpi, Only : mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, mpi_bcast, &
       mpi_comm_world, mpi_integer, mpi_abort
  Use add_tests
  Use diag_tests
  Use multiply_tests
  Use subtract_tests
  Use ks_add_tests
  Use ks_diag_tests
  Use ks_multiply_tests
  Use ks_subtract_tests
  Use ks_extract_tests
  Use ks_unary_tests
  Use test_params

  Implicit None

  Character(len=256), Dimension(:), Allocatable :: args
  Integer, Parameter :: num_test_sets = 11
  Character(len=11), Dimension(num_test_sets) :: test_sets = [ &
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
  Integer :: start_arg
  Integer :: num_args
  Integer :: i

  Logical :: interactive, supplied
  
  Call mpi_init( error )

  Call mpi_comm_size( mpi_comm_world, nproc, error ) 
  Call mpi_comm_rank( mpi_comm_world,    me, error )

  num_args = command_argument_Count()
  Allocate(args(num_args), stat=error)
  If (error .Ne. 0) Then
     Write(*,*) "Failed to allocate args"
     Call mpi_abort(mpi_comm_world, 1, error)
  End If
  Do i = 1, num_args
     Call get_command_Argument(i, args(i))
  End Do

  If( num_args > 0 ) Then
     interactive = args( 1 ) == 'interactive'
  Else
     interactive = .False.
  End If

  supplied = ( num_args > 0 .And. .Not. interactive ) .Or. ( num_args > 1 .And. interactive )

  start_arg = Merge( 1, 2, .Not. interactive )
  If( supplied ) Then
  Else
     args = test_sets
     num_args = Size( test_sets )
  End If

  If( me == 0 .And. interactive ) Then
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
  Do i = start_arg, num_args
     Select Case(Trim(args(i)))
     Case("add")
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
     Case("subtract")
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
     Case("multiply")
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
     Case("diag")
        ! Diagonalisation tests
        If( me == 0 ) Then
           Write( *, title_format ) 'Diagonalisations'
        End If
        Call test_diag_real
        Call test_diag_complex
     Case("ks_multiply")
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
     Case("ks_extract")
        ! Extract
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution extract'
        End If
        Call test_ks_array_extract
        Call test_ks_array_extract_transpose
        ! Unary +/-
     Case("ks_unary")
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution unary +/-'
        End If
        Call test_ks_array_plus
        Call test_ks_array_minus
     Case("ks_add")
        ! Adds
        Call test_ks_split_add_NN
        Call test_ks_split_add_TN
        Call test_ks_split_add_NT
        Call test_ks_split_add_TT
        Call test_ks_split_pre_add_diagonal
        Call test_ks_split_post_add_diagonal
     Case("ks_subtract")
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
     Case("ks_diag")
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
     Case default
        Write( *,* ) "Unknown test: ", Trim(args(i))
     End Select

  End Do
  Call mpi_finalize( error )

End Program test_distributed_matrix
