Module ks_matrix_module

  !! Wrapper module for distributed matrix classes so we can implement an opaque array - 
  !! At this level we only consider 1 k/spin point 

  Use numbers_module           , Only : wp
  Use distributed_matrix_module, Only : distributed_matrix

  Implicit None

  Type, Public :: ks_matrix
     !! Wrapper type for distributed matices, i.e. either real or complex
     Class( distributed_matrix ), Allocatable, Private :: matrix
   Contains
     ! Public methods
     ! Methods at all levels
     Procedure, Public  :: create               => ks_matrix_create                     !! Create a ks_matrix
     Generic  , Public  :: Operator( .Dagger. ) => dagger                               !! Dagger a ks_matrix
     Generic  , Public  :: Operator( * )        => multiply                             !! multiply 2 ks_matrix's
     Generic  , Public  :: set_by_global        => set_global_real, set_global_complex  !! Set elements by global indices
     Generic  , Public  :: get_by_global        => get_global_real, get_global_complex  !! Get elements using global indices
     ! Private implementations
     Procedure, Private :: dagger               => ks_matrix_dagger
     Procedure, Private :: multiply             => ks_matrix_mult
     Procedure, Private :: set_global_real      => ks_matrix_set_global_real
     Procedure, Private :: set_global_complex   => ks_matrix_set_global_complex
     Procedure, Private :: get_global_real      => ks_matrix_get_global_real
     Procedure, Private :: get_global_complex   => ks_matrix_get_global_complex
  End Type ks_matrix

  Public :: ks_matrix_init
  Public :: ks_matrix_comm_to_base
  Public :: ks_matrix_finalise

Contains

  !##########################################################################################################
  ! Non-type bound procedures

  Subroutine ks_matrix_init( nb )

    !! Initialise the matrix system and optionally sets the default blocking factor

    Use distributed_matrix_module, Only : distributed_matrix_init, &
         distributed_matrix_set_default_blocking

    Integer, Intent( In ), Optional :: nb !! Set a default blocking factor

    Call distributed_matrix_init
    
    If( Present( nb ) ) Then
       Call distributed_matrix_set_default_blocking( nb )
    End If

  End Subroutine ks_matrix_init

  Subroutine ks_matrix_comm_to_base( comm, base_matrix )

    !! Converts an MPI communicator into the data structures
    !! required to describe a matrix mapped onto it

    Use distributed_matrix_module, Only : distributed_matrix_comm_to_base, &
         real_distributed_matrix

    Integer             , Intent( In    ) :: comm
    Type   ( ks_matrix ), Intent(   Out ) :: base_matrix

    Allocate( real_distributed_matrix:: base_matrix%matrix )

    Call distributed_matrix_comm_to_base( comm, base_matrix%matrix )
    
  End Subroutine ks_matrix_comm_to_base

  Subroutine ks_matrix_finalise

    !! Finalise the matrix system

    Use distributed_matrix_module, Only : distributed_matrix_finalise

    Call distributed_matrix_finalise
    
  End Subroutine ks_matrix_finalise

  !##########################################################################################################
  ! Type bound procedures

  Subroutine ks_matrix_create( A, is_complex, m, n, source_matrix )

    !! Create a distributed matrix

    Use distributed_matrix_module, Only : real_distributed_matrix, &
         complex_distributed_matrix

    Class  ( ks_matrix ), Intent(   Out ) :: A
    Logical             , Intent( In    ) :: is_complex
    Integer             , Intent( In    ) :: m
    Integer             , Intent( In    ) :: n
    Type   ( ks_matrix ), Intent( In    ) :: source_matrix

    If( is_complex ) Then
      Allocate( complex_distributed_matrix :: A%matrix ) 
    Else
      Allocate( real_distributed_matrix    :: A%matrix ) 
    End If

    Call A%matrix%create( m, n, source_matrix%matrix )
    
  End Subroutine ks_matrix_create

  Function ks_matrix_dagger( A ) Result( tA )
    
    Type( ks_matrix ) :: tA

    Class( ks_matrix ), Intent( In ) :: A 

    tA%matrix = .Dagger. A%matrix

  End Function ks_matrix_dagger

  Function ks_matrix_mult( A, B ) Result( C )
    
    Type( ks_matrix ) :: C

    Class( ks_matrix ), Intent( In ) :: A
    Type ( ks_matrix ), Intent( In ) :: B

    C%matrix = A%matrix * B%matrix

  End Function ks_matrix_mult

  Subroutine ks_matrix_set_global_real( A, m, n, p, q, data )
    
    Class( ks_matrix )             , Intent( InOut ) :: A
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Integer                        , Intent( In    ) :: p
    Integer                        , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Call A%matrix%set_by_global( m, n, p, q, data )

  End Subroutine ks_matrix_set_global_real

  Subroutine ks_matrix_set_global_complex( A, m, n, p, q, data )
    
    Class( ks_matrix )                , Intent( InOut ) :: A
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent( In    ) :: data

    Call A%matrix%set_by_global( m, n, p, q, data )

  End Subroutine ks_matrix_set_global_complex

  Subroutine ks_matrix_get_global_real( A, m, n, p, q, data )
    
    Class( ks_matrix )             , Intent( In    ) :: A
    Integer                        , Intent( In    ) :: m
    Integer                        , Intent( In    ) :: n
    Integer                        , Intent( In    ) :: p
    Integer                        , Intent( In    ) :: q
    Real( wp ), Dimension( m:, p: ), Intent(   Out ) :: data

    Call A%matrix%get_by_global( m, n, p, q, data )

  End Subroutine ks_matrix_get_global_real

  Subroutine ks_matrix_get_global_complex( A, m, n, p, q, data )
    
    Class( ks_matrix )                , Intent( In    ) :: A
    Integer                           , Intent( In    ) :: m
    Integer                           , Intent( In    ) :: n
    Integer                           , Intent( In    ) :: p
    Integer                           , Intent( In    ) :: q
    Complex( wp ), Dimension( m:, p: ), Intent(   Out ) :: data

    Call A%matrix%get_by_global( m, n, p, q, data )

  End Subroutine ks_matrix_get_global_complex

End Module ks_matrix_module
