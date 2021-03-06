module ks_extract_tests
  use test_params
  use mpi
  implicit none
contains

  Subroutine test_ks_array_extract_transpose()

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

  Subroutine test_ks_array_extract_vary()

    ! checks on varying size creates and extracts

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision

    Integer, Parameter :: ns = 2
    Integer, Parameter :: nk = 3

    Type( ks_array ) :: A, A_split
    Type( ks_array ) :: B, C

    Type( ks_array ) :: base_k

    Type replicated_data
       Real   ( wp ), Dimension( :, : ), Allocatable :: A_r
       Complex( wp ), Dimension( :, : ), Allocatable :: A_c
    End type replicated_data

    Type( replicated_data ), Dimension( :, : ), Allocatable :: A_rep

    Complex( wp ), Dimension( :, :    ), Allocatable :: tmp_c

    Real( wp ), Dimension( :, :    ), Allocatable :: tmp_r

    Real( wp ) :: rand
    Real( wp ) :: max_diff

    Integer, Dimension( :, :    ), Allocatable :: sizes
    Integer, Dimension( :, :, : ), Allocatable :: shapes

    Integer, Dimension( 1:3, 1:nk ) :: k_points
    Integer, Dimension(      1:nk ) :: k_types

    Integer :: k, s, ks
    Integer :: n
    Integer :: r1, r2, c1, c2
    Integer :: error

    r1 = 2
    r2 = n_block + 1
    c1 = 3
    c2 = n_block + 2

    n = Max( m, r1, r2, c1, c2 )

    Allocate( A_rep( 1:nk, 1:ns ) )
    Allocate( sizes( 1:2, 1:nk * ns ) )
    Allocate( shapes( 1:2, 1:2, 1:nk * ns ) )
    
    ks = 0
    Do s = 1, ns
       Do k = 1, nk
          ks = ks + 1
          
          n = m + s * k * 10
          
          sizes( :, ks ) = n

          shapes( 1, 1, ks ) = 1
          shapes( 1, 2, ks ) = 1
          shapes( 2, 1, ks ) = Min( m, m / 2 + ks * 2 )
          shapes( 2, 2, ks ) = Min( m, m / 3 + ks * 4 )
          
          Allocate( A_rep( k, s )%A_r( 1:n, 1:n ) )           
          Allocate( A_rep( k, s )%A_c( 1:n, 1:n ) ) 
          
       End Do
    End Do

    If( me == 0 ) Then

       ks = 0
       Do s = 1, ns
          Do k = 1, nk
             ks = ks + 1


             Call Random_number( A_rep( k, s )%A_r )
             
             n = Size( A_rep( k, s )%A_c, Dim = 1 )
             Allocate( tmp_r( 1:n, 1:n ) )
             Call Random_number( tmp_r )
             A_rep( k, s )%A_c( :, : ) = tmp_r
             Call Random_number( tmp_r )
             A_rep( k, s )%A_c( :, : ) = A_rep( k, s )%A_c( :, : ) + Cmplx( 0.0_wp, tmp_r, Kind = wp )
             Deallocate( tmp_r )
             
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

    Do s = 1, ns
       Do k = 1, nk
          Call mpi_bcast( A_rep( k, s )%A_r, Size( A_rep( k, s )%A_r ), mpi_double_precision, 0, mpi_comm_world, error )
          Call mpi_bcast( A_rep( k, s )%A_c, Size( A_rep( k, s )%A_c ), mpi_double_complex  , 0, mpi_comm_world, error )
       End Do
    End Do

    Call ks_array_init( n_block )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base_k )

    Call A%create( sizes, base_k )
    If( verbose ) Then
       Call A%print_info( 'A', 200 )
    End If
    ks = 0
    Do s = 1, ns
       Do k = 1, nk
          ks = ks + 1
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%set_by_global( k_points( :, k ), s, 1, sizes( 1, ks ), 1, sizes( 1, ks ), A_rep( k, s )%A_r )
          Else
             Call A%set_by_global( k_points( :, k ), s, 1, sizes( 1, ks ), 1, sizes( 1, ks ), A_rep( k, s )%A_c )
          End If
       End Do
    End Do

    Call A%split_ks( 2.0_wp, A_split )
    If( verbose ) Then
       Call A_split%print_info( 'A_split', 100 )
    End If


    B = A_split%extract( shapes )
    If( verbose ) Then
       Call B%print_info( 'B', 100 )
    End If

    max_diff = -1.0_wp
    ks = 0
    Do s = 1, ns
       Do k = 1, nk
          ks = ks + 1
          r1 = shapes( 1, 1, ks )
          r2 = shapes( 2, 1, ks )
          c1 = shapes( 1, 2, ks )
          c2 = shapes( 2, 2, ks )
          Allocate( tmp_r( r1:r2, c1:c2 ) )
          Allocate( tmp_c( r1:r2, c1:c2 ) )
          If( k_types( k ) == K_POINT_REAL ) Then
             Call B%get_by_global( k_points( :, k ), s, 1, r2 - r1 + 1, 1, c2 - c1 + 1, tmp_r )
             tmp_r = Abs( tmp_r - A_rep( k, s )%A_r( r1:r2, c1:c2 ) )
          Else
             Call B%get_by_global( k_points( :, k ), s, 1, r2 - r1 + 1, 1, c2 - c1 + 1, tmp_c )
             tmp_r = Abs( tmp_c - A_rep( k, s )%A_c( r1:r2, c1:c2 ) )
          End If
          max_diff = Max( max_diff, Maxval( tmp_r ) )
          Deallocate( tmp_r )
          Deallocate( tmp_c )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split extract vary ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    ! Check multiplies on the varying size arrays work
    C = .Dagger. B
    A = B * C
    If( verbose ) Then
       Call A%print_info( 'A', 100 )
    End If

    max_diff = -1.0_wp
    ks = 0
    Do s = 1, ns
       Do k = 1, nk
          ks = ks + 1
          r1 = shapes( 1, 1, ks )
          r2 = shapes( 2, 1, ks )
          c1 = shapes( 1, 2, ks )
          c2 = shapes( 2, 2, ks )
          Allocate( tmp_r( r1:r2, r1:r2 ) )
          Allocate( tmp_c( r1:r2, r1:r2 ) )
          If( k_types( k ) == K_POINT_REAL ) Then
             Call A%get_by_global( k_points( :, k ), s, 1, r2 - r1 + 1, 1, r2 - r1 + 1, tmp_r )
             tmp_r = Abs( tmp_r - &
                  Matmul( A_rep( k, s )%A_r( r1:r2, c1:c2 ), Transpose( A_rep( k, s )%A_r( r1:r2, c1:c2 ) ) ) )
          Else
             Call A%get_by_global( k_points( :, k ), s, 1, r2 - r1 + 1, 1, r2 - r1 + 1, tmp_c )
             tmp_r = Abs( tmp_c - &
                  Matmul( A_rep( k, s )%A_c( r1:r2, c1:c2 ), Transpose( Conjg( A_rep( k, s )%A_c( r1:r2, c1:c2 ) ) ) ) )
          End If
          max_diff = Max( max_diff, Maxval( tmp_r ) )
          Deallocate( tmp_r )
          Deallocate( tmp_c )
       End Do
    End Do
    If( me == 0 ) Then
       Write( *, error_format ) 'Error in ks_split multiply vary ', max_diff, &
            Merge( passed, FAILED, max_diff < tol )
    End If

    Call ks_array_finalise

  End Subroutine test_ks_array_extract_vary

end module ks_extract_tests

