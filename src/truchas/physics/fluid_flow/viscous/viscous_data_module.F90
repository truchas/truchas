!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE VISCOUS_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define viscous solution solver variables. Note: This module is
  !   declared as public and all the variables are saved.
  !
  ! Contains: None
  !
  ! Author(s): Jim Sicilian (sicilian@lanl.gov)
  !            Ed Dendy     (dendy@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: string_len
  implicit none
  private

  ! Linear solution name (NOTE: ADD to numerics_input_module -> NUMERICS_INPUT_PARALLEL)
  character(string_len), public, save :: viscous_linear_solution
 
  integer, public, save :: ubik_viscous

  ! Current number of iterations to convergence.
  integer, public, save :: viscous_iterations
  integer, public, save :: prelim_viscous_iterations

  ! Current number of preconditioning iterations.
  integer, public, save :: viscous_precond_interations

  ! Physics Namelist Variables
  logical, save, public :: inviscid
  logical, save, public :: stokes

  ! Numerics Namelist Variables
  real(r8), save, public :: viscous_implicitness

  ! Orthogonal mesh matrix and preconditioner
  real(r8), pointer, public, save, dimension(:,:) :: A_Ortho

  ! Face densities and pressure gradients
  real(r8), pointer, public, save, dimension(:,:) :: Face_Viscosity

  real(r8), pointer, public, save, dimension(:,:) :: Stress_Grad_BC

  ! things used by the stress gradient method...
  logical,  public, dimension(:),     allocatable :: Mask
  real(r8), public, dimension(:,:,:), allocatable :: Grad
  real(r8), public, dimension(:,:),   allocatable :: Mu_Face
  real(r8), public, dimension(:),     allocatable :: Normal, Face_Velocity

END MODULE VISCOUS_DATA_MODULE
