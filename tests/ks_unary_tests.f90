module ks_unary_tests
  use test_params
  use mpi
  implicit none

contains

  ! Assignment to a real

  Subroutine test_ks_assign_real()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
      K_POINT_REAL, K_POINT_COMPLEX
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision

    Real( wp ), Parameter:: test_val = 1.23456789_wp
    
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
    Integer :: error

    Allocate( A_r( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( A_c( 1:n, 1:m, 1:nk, 1:ns ) )
    Allocate( tmp_r( 1:n, 1:m ) )
    Allocate( tmp_c( 1:n, 1:m ) )

    A_r = Huge( A_r )
    A_c = Huge( Real( A_c, Kind( A_c ) ) )

    If( me == 0 ) Then

       A_r = test_val
       A_c = test_val
       
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

    Call A%create( n, m, base_k )
    If( verbose ) Then
      Call A%print_info( 'A', 200 )
    End If

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
      Call A_split%print_info( 'A_split', 100 )
    End If

    A_split = test_val

    B = + A_split
    If( verbose ) Then
      Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    Do s = 1, ns
      Do k = 1, nk
        tmp_r = 0.0_wp
        If( k_types( k ) == K_POINT_REAL ) Then
          Call B%get_by_global( k_points( :, k ), s, 1, n, 1, m, tmp_r )
          tmp_r = Abs( tmp_r - A_r( :, :, k, s ) )
        Else
          Call B%get_by_global( k_points( :, k ), s, 1, n, 1, m, tmp_c )
          tmp_r = Abs( tmp_c - A_c( :, :, k, s ) )
        End If
        max_diff = Max( max_diff, Maxval( tmp_r ) )
      End Do
    End Do
    If( me == 0 ) Then
      Write( *, error_format ) 'Error in ks_split assign ', max_diff, &
        Merge( passed, FAILED, Abs( max_diff ) < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_assign_real

  ! Unary plus/minus

  Subroutine test_ks_array_plus()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
      K_POINT_REAL, K_POINT_COMPLEX
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision

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
      K_POINT_REAL, K_POINT_COMPLEX
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision

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
end module ks_unary_tests
