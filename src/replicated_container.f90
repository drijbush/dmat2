Module replicated_container_module

  !! Module to act as union of the various replicated data containers
  
  Use replicated_scalar_container_module, Only : replicated_scalar_container
  Use replicated_1D_container_module    , Only : replicated_1D_container
  Use replicated_2D_container_module    , Only : replicated_2D_container

  Implicit None

  Public :: replicated_scalar_container
  Public :: replicated_1D_container
  Public :: replicated_2D_container
  
  Private
  
End Module replicated_container_module
