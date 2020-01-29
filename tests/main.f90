Program test_distributed_matrix

  Use, intrinsic :: iso_fortran_env, Only : output_unit

!!$  Use mpi, Only : mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, mpi_bcast, &
!!$       mpi_comm_world, mpi_integer, mpi_abort
  Use mpi, Only : mpi_init, mpi_finalize, mpi_comm_rank, mpi_comm_size, &
       mpi_comm_world, mpi_integer, mpi_abort
  Use add_tests, Only : test_real_add_NN, test_real_add_TN, test_real_add_NT, test_real_add_TT, &
       test_real_post_add_diagonal, test_complex_post_add_diagonal, test_complex_add_NN, test_complex_add_TN, &
       test_complex_add_NT, test_complex_add_NT, test_complex_add_TT
  Use diag_tests       , Only : test_diag_real, test_diag_complex
  Use multiply_tests   , Only : test_real_pre_scale_real, &
       test_real_post_scale_real, &
       test_complex_pre_scale_real, &
       test_complex_post_scale_real, &
       test_real_matmul_NN, &
       test_real_matmul_TN, &
       test_real_matmul_NT, &
       test_real_matmul_TT, &
       test_complex_matmul_NN, &
       test_complex_matmul_TN, &
       test_complex_matmul_NT, &
       test_complex_matmul_TT
  Use subtract_tests, Only : test_real_subtract_NN, &
       test_real_subtract_TN, &
       test_real_subtract_NT, &
       test_real_subtract_TT, &
       test_complex_subtract_NN, &
       test_complex_subtract_TN, &
       test_complex_subtract_NT, &
       test_complex_subtract_TT
  Use ks_add_tests
  Use ks_diag_tests
  Use ks_multiply_tests
  Use ks_subtract_tests
  Use ks_extract_tests
  Use ks_unary_tests
  Use test_params

  Implicit None

  Character( len = * ), Dimension( * ), Parameter :: test_sets = [ &
       "add        ", &
       "subtract   ", &
       "multiply   ", &
       "misc       ", &
       "ks_multiply", &
       "ks_extract ", &
       "ks_unary   ", &
       "ks_add     ", &
       "ks_subtract", &
       "ks_misc    " ]


  Integer :: start_arg
  Integer :: num_args
  Integer :: i

  Logical :: interactive, supplied

  Character( len = 256 ), Dimension( : ), Allocatable :: args
  Character( len = :   ), Dimension( : ), Allocatable :: tests
  
  Call mpi_init( error )

  Call mpi_comm_size( mpi_comm_world, nproc, error ) 
  Call mpi_comm_rank( mpi_comm_world,    me, error )

  If( me == 0 ) Then
     Write( *, * ) 'Running on ', nproc, ' processes'
  End If

  num_args = command_argument_Count()
  Allocate( args( 1:num_args ), stat=error )
  If (error .Ne. 0) Then
     Write(*,*) "Failed to allocate args"
     Call mpi_abort( mpi_comm_world, 1, error )
     Stop ! Compiler doesn't know MPI_Abort kills the program, so this stops warnings
          ! About potentially uninited arrays array later on
  End If
  Do i = 1, num_args
     Call get_command_argument( i, args( i ) )
  End Do

  If( num_args > 0 ) Then
     interactive = args( 1 ) == 'interactive'
  Else
     interactive = .False.
  End If

  supplied = ( num_args > 0 .And. .Not. interactive ) .Or. ( num_args > 1 .And. interactive )

  If( supplied ) Then
     start_arg = Merge( 1, 2, .Not. interactive )
     Allocate( tests, Source = args( start_arg:num_args ) )
  Else
     Allocate( tests, Source = test_sets )
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
     Flush( output_unit )
  End If
  ! Adds
  Do i = 1, Size( tests )
     Select Case( Trim( Adjustl( tests( i ) ) ) )
     Case("add")
        ! Additions
        If( me == 0 ) Then
           Write( *, title_format ) 'Additions'
           Flush( output_unit )
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
           Flush( output_unit )
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
           Flush( output_unit )
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
        
     Case("misc")
        ! Diagonalisation and other misc tests
        If( me == 0 ) Then
           Write( *, title_format ) 'Diagonalisations'
           Flush( output_unit )
        End If
        Call test_diag_real
        Call test_diag_complex
        
     Case("ks_multiply")
        ! Multiplies
        If( me == 0 ) Then
           Write( *, title_format ) 'Multiplies'
           Flush( output_unit )
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
           Flush( output_unit )
        End If
        Call test_ks_array_matmul_NN
        Call test_ks_array_matmul_TN
        Call test_ks_array_matmul_NT
        Call test_ks_array_matmul_TT
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution multiplies'
           Flush( output_unit )
        End If
        Call test_ks_split_real_pre_scale
        Call test_ks_split_real_post_scale
        Call test_ks_split_real_pre_multiply_diagonal
        Call test_ks_split_real_post_multiply_diagonal
        Call test_ks_split_real_pre_multiply_diagonal_vary
        Call test_ks_split_real_post_multiply_diagonal_vary
        Call test_ks_split_matmul_NN
        Call test_ks_split_matmul_TN
        Call test_ks_split_matmul_NT
        Call test_ks_split_matmul_TT
        If( me == 0 ) Then
           Write( *, title_format ) 'Split and re-join multiplies'
           Flush( output_unit )
        End If
        Call test_ks_rejoin_matmul_NN
        
        If( me == 0 ) Then
           Write( *, title_format ) 'Split double dot products'
           Flush( output_unit )
        End If
        Call test_ks_split_double_dot
        Call test_ks_split_double_dot_T
        
     Case("ks_extract")
        ! Extract
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution extract'
           Flush( output_unit )
        End If
        Call test_ks_array_extract
        Call test_ks_array_extract_transpose
        Call test_ks_array_extract_vary
        
     Case("ks_unary")
        ! Unary +/-
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution unary +/-'
           Flush( output_unit )
        End If
        Call test_ks_array_plus
        Call test_ks_array_minus
        Call test_ks_assign_real
        
     Case("ks_add")
        ! Adds
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution adds'
           Flush( output_unit )
        End If
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
           Flush( output_unit )
        End If
        Call test_ks_split_subtract_NN
        Call test_ks_split_subtract_TN
        Call test_ks_split_subtract_NT
        Call test_ks_split_subtract_TT
        Call test_ks_split_pre_subtract_diagonal
        Call test_ks_split_post_subtract_diagonal

     Case("ks_misc")
        ! Diags and other misc
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution diags'
           Flush( output_unit )
        End If
        Call test_ks_array_diag
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution Choleskis'
           Flush( output_unit )
        End If
        Call test_ks_array_choleski
        If( me == 0 ) Then
           Write( *, title_format ) 'Split distribution Triangular inverts'
           Flush( output_unit )
        End If
        Call test_ks_array_tr_inv
        Call test_ks_array_tr_inv_with_iterator
     Case default
        Write( *,* ) "Unknown test: ", Trim( tests( i ) )
        Flush( output_unit )
     End Select

  End Do
  Call mpi_finalize( error )

End Program test_distributed_matrix
