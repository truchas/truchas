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
  use kind_module,        only: int_kind, log_kind, real_kind
  use parameter_module,   only: string_len

  implicit none

  ! Private Module
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Linear solution name
  character(LEN = string_len), public, save :: projection_linear_solution

  ! Ubik_user element number to use for projection
  integer(int_kind),           public, save :: UBIK_PRESSURE

  ! Current number of iterations to convergence.
  integer(int_kind),           public, save :: mac_projection_iterations
  integer(int_kind),           public, save :: prelim_projection_iterations
  integer(int_kind),           public, save :: projection_iterations

  ! Current number of preconditioning iterations.
  integer(int_kind),           public, save :: mac_projection_precond_iter
  integer(int_kind),           public, save :: projection_precond_iterations

  ! Dirichlet BC flag
  logical(log_kind),           public, save :: dirichlet_pressure


  ! Orthogonal mesh matrix and preconditioner
  real(real_kind),   pointer, public, save, dimension(:,:)   :: A_Ortho

  ! Rank-1 real coefficient arrays.
  real(real_kind),   pointer, public, save, dimension(:)     :: Solution_Vtx

  ! Rank-2 real coefficient arrays.
  real(real_kind),   pointer, public, save, dimension(:,:)   :: Coeff, Solution_Ngbr, Solution_Vertex

  ! Face densities and pressure gradients
  real(real_kind),   pointer, public, save, dimension(:,:)   :: Face_Density
  real(real_kind),   pointer, public, save, dimension(:,:,:) :: dt_gradP_over_Rho
  real(real_kind),   public, allocatable, dimension(:,:,:)   :: dtRhoG_over_Rho

  ! Array related to boundary conditions.
  integer(int_kind), pointer, public, save, dimension(:,:)   :: Boundary_Flag

  ! Compressibility Array
  real(real_kind),   allocatable, public,   dimension(:)     :: Vol_over_RhoCsqDt

  ! Volume Change Rate Array
  real(real_kind),   allocatable, public,   dimension(:)     :: DVol_by_Dt_over_Vol

  ! pressure head arrays...
  real(KIND = real_kind),   public, allocatable, dimension(:,:)          :: ghc
  real(KIND = real_kind),   public, allocatable, dimension(:,:)          :: ghn

  !-mf
  real(real_kind), public, allocatable, dimension(:,:,:)      :: dtCsf_over_Rho
  real(real_kind), public, allocatable, dimension(:,:,:)      :: Fcsf_new

END MODULE PROJECTION_DATA_MODULE
