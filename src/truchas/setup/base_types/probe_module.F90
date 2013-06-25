MODULE probe_module 
  !=======================================================================
  ! Purpose: 
  !    Define the probe variable structure. 
  !    Probes are defined in the user input by their name and physical positions
  !    Field variables at the probe positions are outputted for diagnostic purposes.
  ! 
  ! Author(s): Sharen Cummins (scummins@lanl.gov)
  !======================================================================= 
 
  use kind_module,      only: int_kind, real_kind 
#ifdef USE_TBROOK
  use tbrook_module,    only: brook
#endif
  use parameter_module, only: ndim, string_len
  use constants_module, only: preset
#ifdef USE_DANU
  use,intrinsic :: iso_c_binding, only: c_ptr
#endif
 
  implicit none 
 
  private 
 
  ! public procedures and structures
 
  public :: probes, probe_node, probe_cell, probe_scalarfield, &
       probe_vectorfield, probe_tensorfield, probe_variable
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
 
  type probe_node
     integer(KIND=int_kind)                :: index
     real(KIND=real_kind), dimension(ndim) :: coords
  end type probe_node
  
  type probe_cell
     integer(KIND=int_kind)                :: index
     real(KIND=real_kind), dimension(ndim) :: coords
  end type probe_cell

  type probe_scalarfield
     real(KIND=real_kind)                 :: field
     character(LEN=string_len)            :: meshspace
  end type probe_scalarfield

  type probe_vectorfield
     real(KIND=real_kind),dimension(ndim) :: field
     character(LEN=string_len)            :: meshspace
  end type probe_vectorfield

  type probe_tensorfield
     real(KIND=real_kind),dimension(6)    :: field
     character(LEN=string_len)            :: meshspace
  end type probe_tensorfield

  type probe_variable  
     character(LEN=string_len)                        :: name         = 'Unnamed' !probe name defined in user input
     real(KIND=real_kind), dimension(ndim)            :: coords                   !coords of the probe location defined in user input
     character(LEN=string_len)                        :: description  = 'None'    !optional description of probe defined by user input
     real(KIND=real_kind)                             :: coords_scale = preset

     type(probe_node)                                 :: node                !nearest node (index,coords) to probe
     type(probe_cell)                                 :: cell                !nearest cell (index,coords) to probe

#ifdef USE_DANU
     type(c_ptr), dimension(:), pointer               :: pid        ! Danu probe data set id
#endif
#ifdef USE_TBROOK
     type(Brook), dimension(:), pointer               ::BrookLU     !list containing all brook pointers to the probe's field output files
#endif
     character(LEN=string_len), dimension(:), pointer ::NameLU      !list containing names of all probe's fields
     type(probe_scalarfield), dimension(:), pointer   ::ScalarVarLU !list containing the probe's scalar fields
     type(probe_vectorfield), dimension(:), pointer   ::VectorVarLU !list containing the probe's vector fields
     type(probe_tensorfield), dimension(:), pointer   ::TensorVarLU !list containing the probe's tensor fields
  end type probe_variable
 
  type(probe_variable), dimension(:), pointer :: probes
                                           
 
END MODULE probe_module
