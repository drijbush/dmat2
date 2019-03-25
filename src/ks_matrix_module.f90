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
     Procedure, Public :: create               => ks_matrix_create                     !! Create a ks_matrix
     Generic  , Public :: Operator( .Dagger. ) => dagger                               !! Dagger a ks_matrix
     Generic  , Public :: Operator( * )        => multiply                             !! multiply 2 ks_matrix's
     Procedure, Public :: size                 => ks_matrix_size                       !! Get the dimensions of the matrix
     Generic  , Public :: set_by_global        => set_global_real, set_global_complex  !! Set elements by global indices
     Generic  , Public :: get_by_global        => get_global_real, get_global_complex  !! Get elements using global indices
     Procedure, Public :: get_comm             => ks_matrix_communicator               !! Get the communicator containing the processes holding the matrix
     Procedure, Public :: global_to_local      => ks_matrix_global_to_local            !! Get an array for mapping global indices to local  ones
     Procedure, Public :: local_to_global      => ks_matrix_local_to_global            !! Get an array for mapping local  indices to global ones
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
  Public :: ks_matrix_remap_data
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

  Subroutine ks_matrix_remap_data( A, parent_communicator, B )

    ! Issues here because either in the ource or remapped
    ! matrix a process may not actually hold any part of it
    ! and so the unallocated actual argument doesn't contain
    ! any information about what it is

    Use distributed_matrix_module, Only : distributed_matrix_remap_data

    Type   ( ks_matrix ), Allocatable, Intent( In    ) :: A
    Integer             ,              Intent( In    ) :: parent_communicator
    Type   ( ks_matrix ), Allocatable, Intent( InOut ) :: B

    Type( ks_matrix ), Allocatable :: dummy_A, dummy_B

    Logical :: p_A, p_B

    p_A = Allocated( A )
    p_B = Allocated( A )

    ! One of A or B must be present to define the data type of what is being redistributed
    If( .Not. p_A .And. .Not. p_B ) Then
       Stop "Must specify one of the matrices in ks_matrix_remap_data"
    End If

    ! Both matrices MUST be of the same type. Thus if one is not specified
    ! create a dummy with the same type as one that is present
    If( .Not. p_A ) Then
       Allocate( dummy_A, Source = B )
    End If
    If( .Not. p_B ) Then
       Allocate( dummy_B, Source = A )
    End If

    ! Now can call remap data routines, carefully indicating which
    ! arguments are dummies because this process

    If     (       p_A .And. p_B ) Then
       Call distributed_matrix_remap_data(       A%matrix, .False., parent_communicator,        B%matrix, .False. )

    Else If( .Not. p_A .And. p_B ) Then
       Call distributed_matrix_remap_data( dummy_A%matrix, .True. , parent_communicator,        B%matrix, .False. )

    Else If(       p_A .And. .Not. p_B ) Then
       Call distributed_matrix_remap_data(       A%matrix, .False. , parent_communicator, dummy_B%matrix, .True.  )

    Else
       Stop "Must specify one of the matrices in ks_matrix_remap_data"
    End If


    
  End Subroutine ks_matrix_remap_data

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

  Function ks_matrix_size( A, dim ) Result( n )

    Integer :: n
    
    Class( ks_matrix ), Intent( In ) :: A
    Integer           , Intent( In ) :: dim

    n = A%matrix%size( dim )
    
  End Function ks_matrix_size

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

  Function ks_matrix_communicator( A ) Result( c )

    !! Get the communicator containing the processes holding the matrix

    Integer :: c

    Class( ks_matrix ), Intent( In ) :: A

    c = A%matrix%get_comm()
    
  End Function ks_matrix_communicator

  Function ks_matrix_global_to_local( A, what ) Result( gl_indexing )

    !! For the given matrix get an array for mapping global indices to local  ones
    !! Get row mapping arrays if WHAT='C' or 'c', row ones if it is equal to 'R' or 'r'

    Integer, Dimension( : ), Allocatable :: gl_indexing
    
    Class( ks_matrix   ), Intent( In ) :: A
    Character( Len = * ), Intent( In ) :: what

    gl_indexing = A%matrix%global_to_local( what )

  End Function ks_matrix_global_to_local
  
  Function ks_matrix_local_to_global( A, what ) Result( lg_indexing )

    !! For the given matrix get an array for mapping local indices to global  ones
    !! Get row mapping arrays if WHAT='C' or 'c', row ones if it is equal to 'R' or 'r'

    Integer, Dimension( : ), Allocatable :: lg_indexing
    
    Class( ks_matrix   ), Intent( In ) :: A
    Character( Len = * ), Intent( In ) :: what

    lg_indexing = A%matrix%local_to_global( what )

  End Function ks_matrix_local_to_global

End Module ks_matrix_module
