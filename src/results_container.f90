Module replicated_result_container_module

  Use numbers_module, Only : wp

  Implicit None

  Type, Abstract, Private :: replicated_result
   Contains
     Generic, Public :: Assignment( = ) => store_real_data, store_complex_data
     Procedure( store_real_data    ), Pass( A ), Deferred, Private :: store_real_data
     Procedure( store_complex_data ), Pass( A ), Deferred, Private :: store_complex_data
     Generic, Public :: Assignment( = ) => get_real_data, get_complex_data
     Procedure(   get_real_data    ), Pass( A ), Deferred, Private :: get_real_data
     Procedure(   get_complex_data ), Pass( A ), Deferred, Private :: get_complex_data
  End type replicated_result

  Type, Private, Extends( replicated_result ) :: real_replicated_result
     Real( wp ), Private :: data
   Contains
     Procedure, Pass( A ), Private :: store_real_data    => store_real_data_into_real
     Procedure, Pass( A ), Private :: store_complex_data => store_complex_data_into_real
     Procedure, Pass( A ), Private :: get_real_data      => get_real_data_from_real
     Procedure, Pass( A ), Private :: get_complex_data   => get_complex_data_from_real
  End type real_replicated_result

  Type, Private, Extends( replicated_result ) :: complex_replicated_result
     Complex( wp ), Private :: data
   Contains
     Procedure, Pass( A ), Private :: store_real_data    => store_real_data_into_complex
     Procedure, Pass( A ), Private :: store_complex_data => store_complex_data_into_complex
     Procedure, Pass( A ), Private :: get_real_data      => get_real_data_from_complex
     Procedure, Pass( A ), Private :: get_complex_data   => get_complex_data_from_complex
  End type complex_replicated_result

  Type, Public :: replicated_result_container
     Class( replicated_result ), Allocatable, Private :: data
   Contains
     Generic, Public :: Assignment( = ) => store_real, store_complex
     Procedure, Pass( A ), Private :: store_real
     Procedure, Pass( A ), Private :: store_complex
     Generic, Public :: Assignment( = ) => get_real, get_complex
     Procedure, Pass( A ), Private :: get_real
     Procedure, Pass( A ), Private :: get_complex
  End type replicated_result_container

  Private

  Abstract Interface
     
     Subroutine store_real_data( A, data )
       Import :: wp
       Import :: replicated_result
       Implicit None
       Class( replicated_result ), Intent(   Out ) :: A
       Real( wp )                , Intent( In    ) :: data
     End Subroutine store_real_data
     
     Subroutine store_complex_data( A, data )
       Import :: wp
       Import :: replicated_result
       Implicit None
       Class( replicated_result ), Intent(   Out ) :: A
       Complex( wp )             , Intent( In    ) :: data
     End Subroutine store_complex_data
     
     Subroutine get_real_data( data, A )
       Import :: wp
       Import :: replicated_result
       Implicit None
       Real( wp )                , Intent(   Out ) :: data
       Class( replicated_result ), Intent( In    ) :: A
     End Subroutine get_real_data
     
     Subroutine get_complex_data( data, A )
       Import :: wp
       Import :: replicated_result
       Implicit None
       Complex( wp )             , Intent(   Out ) :: data
       Class( replicated_result ), Intent( In    ) :: A
     End Subroutine get_complex_data
     
  End Interface

Contains

  Subroutine store_real( A, data )

    Class( replicated_result_container ), Intent(   Out ) :: A
    Real( wp )                          , Intent( In    ) :: data

    Allocate( real_replicated_result :: A%data )

    A%data = data
    
  End Subroutine store_real
  
  Subroutine store_complex( A, data )

    Class( replicated_result_container ), Intent(   Out ) :: A
    Complex( wp )                       , Intent( In    ) :: data

    Allocate( complex_replicated_result :: A%data )
    
    A%data = data

  End Subroutine store_complex
  
  Subroutine store_real_data_into_real( A, data )

    Class( real_replicated_result ), Intent(   Out ) :: A
    Real( wp )                     , Intent( In    ) :: data

    A%data = data
    
  End Subroutine store_real_data_into_real
  
  Subroutine store_complex_data_into_real( A, data )

    Class( real_replicated_result ), Intent(   Out ) :: A
    Complex( wp )                  , Intent( In    ) :: data

    A%data = Real( data, Kind( A%data ) )
    
  End Subroutine store_complex_data_into_real
  
  Subroutine store_real_data_into_complex( A, data )

    Class( complex_replicated_result ), Intent(   Out ) :: A
    Real( wp )                        , Intent( In    ) :: data

    A%data = data
    
  End Subroutine store_real_data_into_complex
  
  Subroutine store_complex_data_into_complex( A, data )

    Class( complex_replicated_result ), Intent(   Out ) :: A
    Complex( wp )                     , Intent( In    ) :: data

    A%data = data
    
  End Subroutine store_complex_data_into_complex
  
  Subroutine get_real( data, A )

    Real( wp )                          , Intent(   Out ) :: data
    Class( replicated_result_container ), Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_real
  
  Subroutine get_complex( data, A )

    Complex( wp )                       , Intent(   Out ) :: data
    Class( replicated_result_container ), Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_complex
  
  Subroutine get_real_data_from_real( data, A )

    Real( wp )                     , Intent(   Out ) :: data
    Class( real_replicated_result ), Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_real_data_from_real
  
  Subroutine get_real_data_from_complex( data, A )

    Real( wp )                        , Intent(   Out ) :: data
    Class( complex_replicated_result ), Intent( In    ) :: A

    data = Real( A%data, Kind( data ) )
    
  End Subroutine get_real_data_from_complex
    
  Subroutine get_complex_data_from_real( data, A )

    Complex( wp )                  , Intent(   Out ) :: data
    Class( real_replicated_result ), Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_complex_data_from_real
  
  Subroutine get_complex_data_from_complex( data, A )

    Complex( wp )                     , Intent(   Out ) :: data
    Class( complex_replicated_result ), Intent( In    ) :: A

    data = A%data
    
  End Subroutine get_complex_data_from_complex
    
End Module replicated_result_container_module
