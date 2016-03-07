!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PROJECTION_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define projection solution solver variables. Note: This module is
  !   declared as public and all the variables are saved.
  !
  ! Contains: None
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: string_len
  implicit none
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Linear solution name
  character(LEN = string_len), public, save :: projection_linear_solution

  ! Ubik_user element number to use for projection
  integer, public, save :: UBIK_PRESSURE

  ! Current number of iterations to convergence.
  integer, public, save :: mac_projection_iterations
  integer, public, save :: prelim_projection_iterations
  integer, public, save :: projection_iterations

  ! Current number of preconditioning iterations.
  integer, public, save :: mac_projection_precond_iter
  integer, public, save :: projection_precond_iterations

  ! Dirichlet BC flag
  logical, public, save :: dirichlet_pressure

  ! Rank-1 real coefficient arrays.
  real(r8), pointer, public, save, dimension(:) :: Solution_Vtx

  ! Rank-2 real coefficient arrays.
  real(r8), allocatable, public, save :: Coeff(:,:)
  real(r8), pointer, public, save, dimension(:,:) :: Solution_Ngbr, Solution_Vertex

  ! Face densities and pressure gradients
  real(r8), pointer, public, save, dimension(:,:)   :: Face_Density
  real(r8), pointer, public, save, dimension(:,:,:) :: dt_gradP_over_Rho
  real(r8), public, allocatable, dimension(:,:,:)   :: dtRhoG_over_Rho

  ! Array related to boundary conditions.
  integer, pointer, public, save, dimension(:,:) :: Boundary_Flag

  ! Compressibility Array
  real(r8), allocatable, public, dimension(:) :: Vol_over_RhoCsqDt

  ! Volume Change Rate Array
  real(r8), allocatable, public, dimension(:) :: DVol_by_Dt_over_Vol

  ! pressure head arrays...
  real(r8), public, allocatable, dimension(:,:) :: ghc
  real(r8), public, allocatable, dimension(:,:) :: ghn

  !-mf
  real(r8), public, allocatable, dimension(:,:,:) :: dtCsf_over_Rho
  real(r8), public, allocatable, dimension(:,:,:) :: Fcsf_new

END MODULE PROJECTION_DATA_MODULE
