!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE probe_module 
  !=======================================================================
  ! Purpose: 
  !    Define the probe variable structure. 
  !    Probes are defined in the user input by their name and physical positions
  !    Field variables at the probe positions are outputted for diagnostic purposes.
  ! 
  ! Author(s): Sharen Cummins (scummins@lanl.gov)
  !======================================================================= 
 
  use kinds, only: r8
  use parameter_module, only: string_len
  use legacy_mesh_api,  only: ndim
  use truchas_h5_outfile, only: th5_probe
  use,intrinsic :: iso_c_binding, only: c_ptr
  implicit none 
  private 
 
  ! public procedures and structures
 
  public :: probes, probe_node, probe_cell, probe_scalarfield, &
       probe_vectorfield, probe_tensorfield, probe_variable
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
 
  type probe_node
     integer :: index
     real(r8), dimension(ndim) :: coords
  end type probe_node
  
  type probe_cell
     integer :: index
     real(r8), dimension(ndim) :: coords
  end type probe_cell

  type probe_scalarfield
     real(r8) :: field
     character(LEN=string_len) :: meshspace
  end type probe_scalarfield

  type probe_vectorfield
     real(r8), dimension(ndim) :: field
     character(LEN=string_len) :: meshspace
  end type probe_vectorfield

  type probe_tensorfield
     real(r8), dimension(6)    :: field
     character(LEN=string_len) :: meshspace
  end type probe_tensorfield

  type probe_variable  
     character(LEN=string_len) :: name         = 'Unnamed' !probe name defined in user input
     real(r8), dimension(ndim) :: coords                   !coords of the probe location defined in user input
     character(LEN=string_len) :: description  = 'None'    !optional description of probe defined by user input
     real(r8) :: coords_scale

     type(probe_node) :: node                !nearest node (index,coords) to probe
     type(probe_cell) :: cell                !nearest cell (index,coords) to probe

     class(th5_probe), dimension(:), pointer :: probe        ! Truchas IO probe data
     character(LEN=string_len), dimension(:), pointer ::NameLU      !list containing names of all probe's fields
     type(probe_scalarfield), dimension(:), pointer   ::ScalarVarLU !list containing the probe's scalar fields
     type(probe_vectorfield), dimension(:), pointer   ::VectorVarLU !list containing the probe's vector fields
     type(probe_tensorfield), dimension(:), pointer   ::TensorVarLU !list containing the probe's tensor fields
  end type probe_variable
 
  type(probe_variable), dimension(:), pointer :: probes
 
END MODULE probe_module
