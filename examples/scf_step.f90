
Module scf_step_module

  Use numbers_module , Only : wp
  Use mpi            , Only : MPI_COMM_WORLD
  Use ks_array_module, Only : ks_array, ks_eval_storage, ks_array_init, ks_array_comm_to_base, ks_array_finalise, &
      K_POINT_REAL, K_POINT_COMPLEX

  Implicit None

  Type repl_matrix
     Real   ( wp ), Dimension( :, : ), Allocatable :: r_data
     Complex( wp ), Dimension( :, : ), Allocatable :: c_data
  End type repl_matrix

Contains

  Subroutine scf_step( ns, k_types, k_points, H, S, ne, P )

    Integer            ,                    Intent( In    ) :: ns
    Integer            , Dimension(    : ), Intent( In    ) :: k_types
    Integer            , Dimension( :, : ), Intent( In    ) :: k_points
    Type( repl_matrix ), Dimension( :, : ), Intent( In    ) :: H
    Type( repl_matrix ), Dimension( :, : ), Intent( In    ) :: S
    Integer                               , Intent( In    ) :: ne
    Type( repl_matrix ), Dimension( :, : ), Intent( InOut ) :: P

    Type( ks_array ) :: base
    Type( ks_array ) :: Hd
    Type( ks_array ) :: Hds, Sds, Pds
    Type( ks_array ) :: Lds, Lds_inv, Lds_inv_T
    Type( ks_array ) :: Qds

    Type( ks_eval_storage ), Dimension( : ), Allocatable :: E
    
    Integer :: n, nk
    Integer :: spin, k

    If( Allocated( H( 1, 1 )%r_data ) ) Then
       n = Size( H( 1, 1 )%r_data, Dim = 1 )
    Else
       n = Size( H( 1, 1 )%c_data, Dim = 1 )
    End If

    nk = Size( k_types )

    ! Initialise the system and get a base matrix to act as a template for the distribution
    ! Note hte base matrix is in "All K-points" mode
    Call ks_array_init( nb = 32 )
    Call ks_array_comm_to_base( MPI_COMM_WORLD, ns, k_types, k_points, base )

    ! Using the base matrix to describe the distribution create a new Hamiltonian matrix
    Call Hd%create( n, n, base )
    ! Split this matrix so that different k points and spins are held by different processes,
    ! so allowing dual level paralleism
    Call Hd%split_ks( 2.0_wp, Hds )
    ! Alternatively just take a copy so all that follows now happens in "All K-point"
!!$    Hds = Hd

    ! Create  the overlap matrix with the same distribution as the Hamiltonian matrix
    Call Sds%create( n, n, Hds )

    ! Load up the overlap and hamiltonian matrices
    Do spin = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call Hds%set_by_global( k_points( :, k ), spin, 1, n, 1, n, H( k, spin )%r_data )
             Call Sds%set_by_global( k_points( :, k ), spin, 1, n, 1, n, S( k, spin )%r_data )
          Else
             Call Hds%set_by_global( k_points( :, k ), spin, 1, n, 1, n, H( k, spin )%c_data )
             Call Sds%set_by_global( k_points( :, k ), spin, 1, n, 1, n, S( k, spin )%c_data )
          End If
       End Do
    End Do

    ! Generate the orthogonalising transform
    Lds       = .Choleski. Sds
    Lds_inv   = .TrInv.    Lds
    Lds_inv_T = .Dagger.   Lds_inv

    ! Transfrom AO -> MO
    Hds = Lds_inv * Hds * Lds_inv_T

    ! Diagonalise in the MO basis
    Allocate( E( 1:nk * ns ) )
    Call Hds%diag( Qds, E )

    ! Transform the Evecs MO -> AO
    Qds = Lds_inv_T * Qds

    ! Extract the evecs for the occupied states into a new distributed matrix
    Qds = Qds%extract( 1, n, 1, ne )

    ! Form the density matrix - note dealing with rectangular matrices
    Pds = Qds * .Dagger. Qds

    ! Report some information about the layout of the resulting density matrix
    Call Pds%print_info( 'Pds ', 9999 )

    ! Re-replicate into the result array
    Do spin = 1, ns
       Do k = 1, nk
          If( k_types( k ) == K_POINT_REAL ) Then
             Call Pds%get_by_global( k_points( :, k ), spin, 1, n, 1, n, P( k, spin )%r_data )
          Else
             Call Pds%get_by_global( k_points( :, k ), spin, 1, n, 1, n, P( k, spin )%c_data )
          End If
       End Do
    End Do

    ! Close down the system
    Call ks_array_finalise
    
  End Subroutine scf_step
  
End Module scf_step_module

Program scf_example

  Use numbers_module , Only : wp
  Use mpi
  Use scf_step_module, Only : scf_step, repl_matrix
  Use ks_array_module, Only : K_POINT_REAL, K_POINT_COMPLEX
  
  Implicit None

  Type( repl_matrix ), Dimension( :, : ), Allocatable :: H
  Type( repl_matrix ), Dimension( :, : ), Allocatable :: S
  Type( repl_matrix ), Dimension( :, : ), Allocatable :: P
  Type( repl_matrix ), Dimension( :, : ), Allocatable :: P_lapack

  Complex( wp ), Dimension( :, : ), Allocatable :: Qc

  Complex( wp ), Dimension( : ), Allocatable :: cwork
  
  Real( wp ), Dimension( :, : ), Allocatable :: Qr
  Real( wp ), Dimension( :, : ), Allocatable :: rtmp1, rtmp2

  Real( wp ), Dimension( : ), Allocatable :: E
  Real( wp ), Dimension( : ), Allocatable :: rwork
  Real( wp ), Dimension( : ), Allocatable :: rwork2
  
  Real( wp ) :: scal_tmp

  Integer, Dimension(    : ), Allocatable :: k_types
  Integer, Dimension( :, : ), Allocatable :: k_points

  Integer :: me, nproc
  Integer :: n, ns, nk, ne
  Integer :: spin, k
  Integer :: i
  Integer :: error

  Call MPI_Init( error )
  Call MPI_Comm_size( MPI_COMM_WORLD, nproc, error )
  Call MPI_Comm_rank( MPI_COMM_WORLD, me   , error )

  ! Read the input
  If( me == 0 ) Then
     Write( *, * ) 'n, ns, nk, ne'
     Read ( *, * )  n, ns, nk, ne
  End If
  Call MPI_Bcast( n , 1, MPI_Integer, 0, MPI_COMM_World, error )
  Call MPI_Bcast( ns, 1, MPI_Integer, 0, MPI_COMM_World, error )
  Call MPI_Bcast( nk, 1, MPI_Integer, 0, MPI_COMM_World, error )
  Call MPI_Bcast( ne, 1, MPI_Integer, 0, MPI_COMM_World, error )

  Allocate( rtmp1( 1:n, 1:n ) )
  Allocate( rtmp2( 1:n, 1:n ) )

  ! Define some k_points and their types
  Allocate( k_types (      1:nk ) )
  Allocate( k_points( 1:3, 1:nk ) )
  If( me == 0 ) Then
     Do k = 1, nk
        k_points( :, k ) = [ k, 0, 0 ]
     End Do
     k_types = Merge( K_POINT_REAL, K_POINT_COMPLEX, scal_tmp > 0.6_wp )
     Do While( ( All( k_types == K_POINT_REAL ) .Or. All( k_types == K_POINT_COMPLEX ) ) .And. nk /= 1 )
        Do k = 1, nk
           Call Random_number( scal_tmp )
           k_types ( k ) = Merge( K_POINT_REAL, K_POINT_COMPLEX, scal_tmp > 0.6_wp )
        End Do
     End Do
  End If
  Call MPI_Bcast( k_types , Size( k_types  ), MPI_Integer, 0, MPI_COMM_World, error )
  Call MPI_Bcast( k_points, Size( k_points ), MPI_Integer, 0, MPI_COMM_World, error )
     

  ! Now set up the S and H matrices, and allocate the P matrix
  Allocate( H( 1:nk, 1:ns ) )
  Allocate( S( 1:nk, 1:ns ) )
  Allocate( P( 1:nk, 1:ns ) )
  Allocate( P_lapack( 1:nk, 1:ns ) )
  Do spin = 1, ns
     Do k = 1, nk
        If( k_types( k ) == K_POINT_REAL ) Then
           Allocate( H( k, spin )%r_data( 1:n, 1:n ) )
           Allocate( S( k, spin )%r_data( 1:n, 1:n ) )
           Allocate( P( k, spin )%r_data( 1:n, 1:n ) )
           Allocate( P_lapack( k, spin )%r_data( 1:n, 1:n ) )
           If( me == 0 ) Then
              ! H Matrix - symmetric and make it somewhat diagonal dominant
              Call Random_number( H( k, spin )%r_data )
              H( k, spin )%r_data = 0.05_wp * ( H( k, spin )%r_data + Transpose( H( k, spin )%r_data ) )
              Do i = 1, n
                 H( k, spin )%r_data( i, i ) = H( k, spin )%r_data( i, i ) * 100.0_wp
              End Do
              ! S Matrix - same for all s but needs to be symmetric positive definite and unit diagonal
              If( spin == 1 ) Then
                 Call Random_number( S( k, spin )%r_data )
                 S( k, spin )%r_data = 0.05_wp * ( S( k, spin )%r_data + Transpose( S( k, spin )%r_data ) )
                 Do i = 1, n
                    S( k, spin )%r_data( i, i ) = 1.0_wp
                 End Do
                 S( k, spin )%r_data = S( k, spin )%r_data * Transpose( S( k, spin )%r_data )
              Else
                 S( k, spin )%r_data = S( k, 1 )%r_data
              End If
           End If
           Call MPI_Bcast( H( k, spin )%r_data, Size( H( k, spin )%r_data ), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error )
           Call MPI_Bcast( S( k, spin )%r_data, Size( S( k, spin )%r_data ), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error )
        Else
           Allocate( H( k, spin )%c_data( 1:n, 1:n ) )
           Allocate( S( k, spin )%c_data( 1:n, 1:n ) )
           Allocate( P( k, spin )%c_data( 1:n, 1:n ) )
           Allocate( P_lapack( k, spin )%c_data( 1:n, 1:n ) )
           If( me == 0 ) Then
              ! H Matrix - symmetric and make it somewhat diagonal dominant
              Call Random_number( rtmp1 )
              Call Random_number( rtmp2 )
              H( k, spin )%c_data = Cmplx( rtmp1, rtmp2, Kind = wp )
              H( k, spin )%c_data = 0.05_wp * ( H( k, spin )%c_data + Conjg( Transpose( H( k, spin )%c_data ) ) )
              Do i = 1, n
                 H( k, spin )%c_data( i, i ) = H( k, spin )%c_data( i, i ) * 100.0_wp
              End Do
              ! S Matrix - same for all s but needs to be hermitian positive definite and unit diagonal
              If( spin == 1 ) Then
                 Call Random_number( rtmp1 )
                 Call Random_number( rtmp2 )
                 S( k, spin )%c_data = Cmplx( rtmp1, rtmp2, Kind = wp )
                 S( k, spin )%c_data = 0.05_wp * ( S( k, spin )%c_data + Conjg( Transpose( S( k, spin )%c_data ) ) )
                 Do i = 1, n
                    S( k, spin )%c_data( i, i ) = 1.0_wp
                 End Do
                 S( k, spin )%c_data = S( k, spin )%c_data * Conjg( Transpose( S( k, spin )%c_data ) )
              Else
                 S( k, spin )%c_data = S( k, 1 )%c_data
              End If
           End If
           Call MPI_Bcast( H( k, spin )%c_data, Size( H( k, spin )%c_data ), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, error )
           Call MPI_Bcast( S( k, spin )%c_data, Size( S( k, spin )%c_data ), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, error )
        End If
     End Do
  End Do

  Call scf_step( ns, k_types, k_points, H, S, ne, P )

  Allocate( Qr( 1:n, 1:n ) )
  Allocate( Qc( 1:n, 1:n ) )
  Allocate( E ( 1:n ) )
  Allocate( rwork( 1:n * 100 ) )
  Allocate( cwork( 1:n * 100 ) )
  Allocate( rwork2( 1:n * 3 ) )
  ! Now generate P_lapack in serial via Scalapack routines
  Do spin = 1, ns
     Do k = 1, nk
        If( k_types( k ) == K_POINT_REAL ) Then
           Qr = H( k, spin )%r_data
           Call dsygv( 1, 'V', 'U', n, Qr, n, S( k, spin )%r_data, n, E, rwork, Size( rwork ), error )
           P_lapack( k, spin )%r_data = Matmul( Qr( :, 1:ne ), Transpose( Qr( :, 1:ne ) ) )
        Else
           Qc = H( k, spin )%c_data
           Call zhegv( 1, 'V', 'U', n, Qc, n, S( k, spin )%c_data, n, E, cwork, Size( rwork ), rwork2, error )
           P_lapack( k, spin )%c_data = Matmul( Qc( :, 1:ne ), Transpose( Conjg( Qc( :, 1:ne ) ) ) )
        End If
     End Do
  End Do

  ! Now check the results
  If( me == 0 ) Then
     Write( *, '( a )' ) 'Maximum errors in density matrices'
     Write( *, '( a, t10, a )' ) 'Spin', 'K-point'
     Do spin = 1, ns
        Do k = 1, nk
           If( k_types( k ) == K_POINT_REAL ) Then
              Write( *, '( i1, t10, i3, t20, g20.10 )' ) spin, k, Maxval( Abs( P( k, spin )%r_data - P_lapack( k, spin )%r_data ) )
           Else
              Write( *, '( i1, t10, i3, t20, g20.10 )' ) spin, k, Maxval( Abs( P( k, spin )%c_data - P_lapack( k, spin )%c_data ) )
           End If
        End Do
     End Do
  End If

  ! Bye bye!
  Call MPI_Finalize( error )

End Program scf_example
