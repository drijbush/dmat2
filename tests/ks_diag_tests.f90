module ks_diag_tests
  use test_params
  use mpi
  implicit none
contains

  Subroutine test_ks_array_diag()

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX, ks_array_replicated_1D
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: Q
    Type( ks_array ) :: QT
    Type( ks_array ) :: B

    Type( ks_array_replicated_1D ), Dimension( 1:nk * ns ) :: E

    Type( ks_array ) :: base_k

    Complex( wp ), Dimension( :, :, :, : ), Allocatable :: A_c
    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: A_r
    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ), Dimension( : ), Allocatable :: evals_ks

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
          evals_ks = E( k + ( s - 1 ) * nk )
          If( k_types( k ) == K_POINT_REAL ) Then
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_r )
             Do i = 1, n
                tmp_r( i, i ) = Abs( tmp_r( i, i ) - evals_ks( i ) )
             End Do
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, n, 1, n, tmp_c )
             Do i = 1, n
                tmp_r( i, i ) = Abs( tmp_c( i, i ) - evals_ks( i ) )
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
         K_POINT_REAL, K_POINT_COMPLEX
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision

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
         K_POINT_REAL, K_POINT_COMPLEX
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision

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
         K_POINT_REAL, K_POINT_COMPLEX, K_POINT_NOT_EXIST
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision

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
end module ks_diag_tests

