!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PARAMETER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all integer and scalar parameters.
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe, LANL (dbk@lanl.gov)
  !
  !=======================================================================
  implicit none
  public
  
  save

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Number of physical dimensions
  integer, parameter :: ndim = 3

  ! Number of rotation axes
  integer, parameter :: nrot = 2*ndim - 3

  ! Maximum number of outputs
  integer, parameter :: mops = 21

  ! Maximum output character string array size
  integer, parameter :: string_dim = 30, string_len = 256

  ! Maximum number of BC namelists and surfaces
  integer, parameter :: mbc_surfaces = 100
  integer, parameter :: mbcsrf       = 5
  integer, parameter :: mbc_nodes    = 50
  !  The radiation BC requires max_bc_dof be at least 4 greater than the number of time,reference-temperature pairs
  integer, parameter :: max_bc_dof   = 24 ! Maximum number of degrees of freedom any kind of BC
  

  ! Number of different BC variable character string forms. 
  integer, parameter :: bc_forms = 30

  ! Number of interface topology character strings allowed
  integer, parameter :: max_topology_models = 5

  ! Current number of BC types and variables
  integer, parameter :: nbcs = 32, nvar = 10

  ! Interface (body) initialization parameters
  integer, parameter :: mtype = 10
  integer, parameter :: msurf = 16
  integer, parameter :: mbody = 50
  integer, parameter :: mtab  = 50
  integer, parameter :: mcoef =  3
  integer, parameter :: mphi  =  5

  integer, parameter :: mregion =  50

  ! Maximum number of material constants          
  integer, parameter :: maxmat = 64
  integer, parameter :: maxcon = 10
  integer, parameter :: maxform = 8

  ! Maximum number of mesh domains
  integer, parameter :: max_domains = 4196

  ! Number of stress/strain components
  integer, Parameter :: ncomps = (ndim*(ndim + 1))/2

  ! Array sizes
  integer, dimension(ndim) :: Nx, Mx, Nx_tot, Mx_tot
  integer :: nicells, nicells_tot, nfaces, boundary_faces, boundary_faces_tot

  ! maximum number of probes allowed in input file
  integer, parameter :: MAX_PROBES = 20
  integer            :: nprobes

  ! Maximum time intervals for electromagnetics

  integer, parameter, public :: MAXSV = 32

END MODULE PARAMETER_MODULE
