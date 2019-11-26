module ks_add_tests
  use test_params
  use mpi
  implicit none
contains
  ! Add tests
  Subroutine test_ks_split_add_NN

    Use numbers_module , Only : wp
    Use ks_array_module, Only : ks_array, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
         K_POINT_REAL, K_POINT_COMPLEX
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision


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
       If( nk /= 1 ) Then
          Do While( All( k_types == K_POINT_REAL ) .Or. All( k_types == K_POINT_COMPLEX ) )
             Do kpoint = 1, nk
                k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
                Call Random_number( rand )
                k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
             End Do
          End Do
       End If
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
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision


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
       If( nk /= 1 ) Then
          Do While( All( k_types == K_POINT_REAL ) .Or. All( k_types == K_POINT_COMPLEX ) )
             Do kpoint = 1, nk
                k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
                Call Random_number( rand )
                k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
             End Do
          End Do
       End If
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
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision


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
       If( nk /= 1 ) Then
          Do While( All( k_types == K_POINT_REAL ) .Or. All( k_types == K_POINT_COMPLEX ) )
             Do kpoint = 1, nk
                k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
                Call Random_number( rand )
                k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
             End Do
          End Do
       End If
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
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision


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
       If( nk /= 1 ) Then
          Do While( All( k_types == K_POINT_REAL ) .Or. All( k_types == K_POINT_COMPLEX ) )
             Do kpoint = 1, nk
                k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
                Call Random_number( rand )
                k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
             End Do
          End Do
       End If
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
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision


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
       If( nk /= 1 ) Then
          Do While( All( k_types == K_POINT_REAL ) .Or. All( k_types == K_POINT_COMPLEX ) )
             Do kpoint = 1, nk
                k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
                Call Random_number( rand )
                k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
             End Do
          End Do
       End If
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
!!$    Use mpi            , Only : mpi_bcast, mpi_comm_world, mpi_double_complex, mpi_double_precision
    Use mpi            , Only : mpi_comm_world, mpi_double_complex, mpi_double_precision


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
       If( nk /= 1 ) Then
          Do While( All( k_types == K_POINT_REAL ) .Or. All( k_types == K_POINT_COMPLEX ) )
             Do kpoint = 1, nk
                k_points( :, kpoint ) = [ kpoint - 1, 0, 0 ]
                Call Random_number( rand )
                k_types( kpoint ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, rand > 0.5_wp )
             End Do
          End Do
       End If
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
end module ks_add_tests

