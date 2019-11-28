Module replicated_result_container_module

  Use numbers_module, Only : wp

  Implicit None

  Type, Abstract, Public :: replicated_result
   Contains
  End type replicated_result

  Type, Public, Extends( replicated_result ) :: real_replicated_result
     Real( wp ) :: data
  End type real_replicated_result

  Type, Public, Extends( replicated_result ) :: complex_replicated_result
     Real( wp ) :: complex
  End type complex_replicated_result
  
End Module replicated_result_container_module
