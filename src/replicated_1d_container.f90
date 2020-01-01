Module replicated_1D_container_module

  !! Implements a simple container which complex or real numbers can be stored in
  !! This module deals with one dimensional arrays
  !! The primary use in dmat2 is to form array of relicated inputs or outputs to
  !! routines where a mix of real and complex values in the arrays are required
  
  Use numbers_module, Only : wp

  Implicit None

  Type, Abstract, Private :: replicated_1D
   Contains
     Generic, Public :: Assignment( = ) => store_real_data, store_complex_data
     Procedure( store_real_data    ), Pass( A ), Deferred, Private :: store_real_data
     Procedure( store_complex_data ), Pass( A ), Deferred, Private :: store_complex_data
     Generic, Public :: Assignment( = ) => get_real_data, get_complex_data
     Procedure(   get_real_data    ), Pass( A ), Deferred, Private :: get_real_data
     Procedure(   get_complex_data ), Pass( A ), Deferred, Private :: get_complex_data
  End type replicated_1D

  Type, Private, Extends( replicated_1D ) :: real_replicated_1D
     Real( wp ), Dimension( : ), Allocatable, Private :: data
   Contains
     Procedure, Pass( A ), Private :: store_real_data    => store_real_data_into_real
     Procedure, Pass( A ), Private :: store_complex_data => store_complex_data_into_real
     Procedure, Pass( A ), Private :: get_real_data      => get_real_data_from_real
     Procedure, Pass( A ), Private :: get_complex_data   => get_complex_data_from_real
  End type real_replicated_1D

  Type, Private, Extends( replicated_1D ) :: complex_replicated_1D
     Complex( wp ), Dimension( : ), Allocatable, Private :: data
   Contains
     Procedure, Pass( A ), Private :: store_real_data    => store_real_data_into_complex
     Procedure, Pass( A ), Private :: store_complex_data => store_complex_data_into_complex
     Procedure, Pass( A ), Private :: get_real_data      => get_real_data_from_complex
     Procedure, Pass( A ), Private :: get_complex_data   => get_complex_data_from_complex
  End type complex_replicated_1D

  Type, Public :: replicated_1D_container
     Class( replicated_1D ), Allocatable, Private :: data
   Contains
     Generic, Public :: Assignment( = ) => store_real, store_complex
     Procedure, Pass( A ), Private :: store_real
     Procedure, Pass( A ), Private :: store_complex
     Generic, Public :: Assignment( = ) => get_real, get_complex
     Procedure, Pass( A ), Private :: get_real
     Procedure, Pass( A ), Private :: get_complex
  End type replicated_1D_container

  Private

  Abstract Interface

     Subroutine store_real_data( A, data )
       Import :: wp
       Import :: replicated_1D
       Implicit None
       Class( replicated_1D )                , Intent( InOut ) :: A
       Real ( wp )           , Dimension( : ), Intent( In    ) :: data
     End Subroutine store_real_data
     
     Subroutine store_complex_data( A, data )
       Import :: wp
       Import :: replicated_1D
       Implicit None
       Class  ( replicated_1D ),                 Intent( InOut ) :: A
       Complex( wp )           , Dimension( : ), Intent( In    ) :: data
     End Subroutine store_complex_data
     
     Subroutine get_real_data( data, A )
       Import :: wp
       Import :: replicated_1D
       Implicit None
       Real( wp )            , Dimension( : ), Allocatable, Intent(   Out ) :: data
       Class( replicated_1D ),                              Intent( In    ) :: A
     End Subroutine get_real_data
     
     Subroutine get_complex_data( data, A )
       Import :: wp
       Import :: replicated_1D
       Implicit None
       Complex( wp )         , Dimension( : ), Allocatable, Intent(   Out ) :: data
       Class( replicated_1D ),                              Intent( In    ) :: A
     End Subroutine get_complex_data
     
  End Interface

Contains

  ! Top level store routines
  
  Subroutine store_real( A, data )

    Class( replicated_1D_container ),                 Intent( InOut ) :: A
    Real( wp )                      , Dimension( : ), Intent( In    ) :: data

    If( Allocated( A%data ) ) Deallocate( A%data )
    Allocate( real_replicated_1D :: A%data )

    A%data = data
    
  End Subroutine store_real
  
  Subroutine store_complex( A, data )

    Class( replicated_1D_container ),                 Intent( InOut ) :: A
    Complex( wp )                   , Dimension( : ), Intent( In    ) :: data

    If( Allocated( A%data ) ) Deallocate( A%data )
    Allocate( complex_replicated_1D :: A%data )
    
    A%data = data

  End Subroutine store_complex

  ! Storing into a real 1D
  
  Subroutine store_real_data_into_real( A, data )

    Class( real_replicated_1D ),                 Intent( InOut ) :: A
    Real( wp )                 , Dimension( : ), Intent( In    ) :: data

    A%data = data
    
  End Subroutine store_real_data_into_real
  
  Subroutine store_complex_data_into_real( A, data )

    Class( real_replicated_1D ),                 Intent( InOut ) :: A
    Complex( wp )              , Dimension( : ), Intent( In    ) :: data

    A%data = Real( data, Kind( A%data ) )
    
  End Subroutine store_complex_data_into_real
  
  ! Storing into a complex 1D

  Subroutine store_real_data_into_complex( A, data )

    Class( complex_replicated_1D ),                 Intent( InOut ) :: A
    Real( wp )                    , Dimension( : ), Intent( In    ) :: data

    A%data = data
    
  End Subroutine store_real_data_into_complex
  
  Subroutine store_complex_data_into_complex( A, data )

    Class( complex_replicated_1D ),                 Intent( InOut ) :: A
    Complex( wp )                 , Dimension( : ), Intent( In    ) :: data

    A%data = data
    
  End Subroutine store_complex_data_into_complex

  ! Top level get routines
  
  Subroutine get_real( data, A )

    Real( wp )                      , Dimension( : ), Allocatable, Intent(   Out ) :: data
    Class( replicated_1D_container )                             , Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_real

  Subroutine get_complex( data, A )

    Complex( wp )                   , Dimension( : ), Allocatable, Intent(   Out ) :: data
    Class( replicated_1D_container ),                              Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_complex
  
  ! Getting into a real 1D
  
  Subroutine get_real_data_from_real( data, A )

    Real( wp )                 , Dimension( : ), Allocatable, Intent(   Out ) :: data
    Class( real_replicated_1D ),                              Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_real_data_from_real
  
  Subroutine get_real_data_from_complex( data, A )

    Real( wp )                    , Dimension( : ), Allocatable, Intent(   Out ) :: data
    Class( complex_replicated_1D ),                              Intent( In    ) :: A

    data = Real( A%data, Kind( data ) )
    
  End Subroutine get_real_data_from_complex

  ! Getting into a complex 1D
  
  Subroutine get_complex_data_from_real( data, A )

    Complex( wp )              , Dimension( : ), Allocatable, Intent(   Out ) :: data
    Class( real_replicated_1D ),                              Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_complex_data_from_real
  
  Subroutine get_complex_data_from_complex( data, A )

    Complex( wp )                 , Dimension( : ), Allocatable, Intent(   Out ) :: data
    Class( complex_replicated_1D ),                              Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_complex_data_from_complex
    
End Module replicated_1D_container_module
