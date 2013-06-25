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
  use kind_module,        only: int_kind, log_kind, real_kind
  use parameter_module,   only: string_len

  implicit none

  ! Private Module
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Linear solution name (NOTE: ADD to numerics_input_module -> NUMERICS_INPUT_PARALLEL)
  character(LEN = string_len),  public, save :: viscous_linear_solution

  integer(int_kind), public, save :: ubik_viscous

  ! Current number of iterations to convergence.
  integer(int_kind), public, save :: viscous_iterations
  integer(int_kind), public, save :: prelim_viscous_iterations

  ! Current number of preconditioning iterations.
  integer(int_kind), public, save :: viscous_precond_interations

  ! Physics Namelist Variables
  logical(KIND = log_kind), save, public :: inviscid

  logical(KIND = log_kind), save, public :: stokes

  ! Numerics Namelist Variables
  real(KIND = real_kind),   save, public :: viscous_implicitness

  ! Orthogonal mesh matrix and preconditioner
  real(real_kind),   pointer, public, save, dimension(:,:)   :: A_Ortho

  ! Face densities and pressure gradients
  real(real_kind),   pointer, public, save, dimension(:,:)   :: Face_Viscosity



  real(real_kind),   pointer, public, save, dimension(:,:)   :: Stress_Grad_BC



  ! things used by the stress gradient method...

  logical(KIND = log_kind), public, dimension(:),     allocatable :: Mask
  real(KIND = real_kind), public, dimension(:,:,:), allocatable :: Grad
  real(KIND = real_kind), public, dimension(:,:),   allocatable :: Mu_Face
  real(KIND = real_kind), public, dimension(:),     allocatable :: Normal, Face_Velocity


END MODULE VISCOUS_DATA_MODULE
