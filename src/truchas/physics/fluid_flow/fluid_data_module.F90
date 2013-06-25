MODULE FLUID_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables needed to solve the Navier-Stokes equations.
  !
  ! Contains: (none)
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Jim Sicilian, LANL CCS-2 (sicilian@lanl.gov)   Jan 2003
  !            M. A. Christon, LANL CCS-2 (christon@lanl.gov) Sep 2007
  !
  !=======================================================================
  use kind_module,      only: log_kind, real_kind, int_kind
  use parameter_module, only: maxmat, ndim

  implicit none

  ! Private Module
  private

  public read_flow_data

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Flag for enabling/disabling fluid flow.
  logical(KIND = log_kind),        public, save                    :: fluid_flow
  logical(KIND = log_kind),        public, save                    :: fluid_to_move
  logical(KIND = log_kind),        public, save                    :: mass_limiter
  real(KIND = real_kind),          public, save                    :: mass_limiter_cutoff

  ! Flag for applying a prescribed velocity field.
  logical(KIND = log_kind),        public, save                    :: applyflow

  ! Flag that defines whether a material flows.
  logical(KIND = log_kind),        public, save, dimension(maxmat) :: isImmobile

  ! Flag that defines the existance of a void material in this calculation and its index
  logical(log_kind),               public, save                    :: Void_Material_Exists
  integer(int_kind),               public, save, dimension(maxmat) :: Void_Material_Index
  integer(int_kind),               public, save                    :: Void_material_Count

  ! Value of pressure in the void.
  real(KIND = real_kind),          public, save                    :: void_pressure

  ! Boussinesq approximation
  logical(KIND = log_kind),        public, save                    :: boussinesq_approximation

  ! Total volume flowing into and out of the domain.
  real(KIND = real_kind),          public, save                    :: qin, qout

  ! Volume Fraction Cutoff Value  (no flow solution if Vof < fluid_cutoff)
  real(KIND = real_kind),          public, parameter               :: fluid_cutoff = 0.01

  ! Time-weighting for treatment of momentum deposition due to solidification
  real(real_kind),                 public, save        :: momentum_solidify_implicitness

  ! Flow arrays.
  real(KIND = real_kind), pointer, public, save, dimension(:,:)     :: Fluxing_Velocity
  real(KIND = real_kind), pointer, public, save, dimension(:,:)     :: Face_Interpolation_Factor
  real(KIND = real_kind), pointer, public, save, dimension(:,:)     :: Momentum_by_Volume
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: fluidVof, fluidVof_n
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: fluidRho
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: fluidRho_n
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: fluidDeltaRho
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: realfluidVof
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: cutRho

! Volume-fraction averaged density at n, n+1 for body force terms
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: avgRho
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: avgRho_n

  real(KIND = real_kind), allocatable, public, save, dimension(:,:) :: Drag_Coefficient
  real(KIND = real_kind), pointer, public, save, dimension(:)       :: courant

  ! Arrays related to partially/totally solidified cells.
  logical(log_kind),      pointer, public, save, dimension(:)      :: isPureImmobile
  logical(log_kind),      pointer, public, save, dimension(:,:)    :: Solid_Face
 
  ! Non-Void Cell list
  logical(KIND = log_kind), pointer, public, save, dimension(:)   :: Cell_isnt_Void
  logical(KIND = log_kind), pointer, public, save, dimension(:,:) :: Ngbr_isnt_Void

  ! Projection / Predictor Communication Arrays
  real(KIND = real_kind),   public, allocatable, save, dimension(:,:) :: Centered_GradP_Dynamic
  real(KIND = real_kind),   public, allocatable, save, dimension(:,:) :: Rho_Face, Rho_Face_n
  real(real_kind), public, save                                       :: MinFluidRho, MinFaceFraction

  ! Storage moved to a single-allocation per run, rather than once per time-step (MAC)
  real(real_kind), public, allocatable, save, dimension(:,:)  :: Mom_Delta
  real(real_kind), public, allocatable, save, dimension(:)    :: pc_delta_rho
  real(real_kind),          public, save, dimension(ndim) :: minVel, maxVel   

 contains 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_FLOW_DATA
 !!
 !! Jim Sicilian <sicilian@lanl.gov>
 !! 20 Apr 2006
 !!
 !! This subroutine reads the fluid flow data from a restart file opened (and pre-
 !! positioned) on UNIT, and initializes the Fluxing_Velocity with this
 !! data (properly distributed and permuted).  VERSION is the version number
 !! of the restart file format.
 !!

  subroutine read_flow_data (unit, version)

    use mesh_module, only: pcell => unpermute_mesh_vector
    use restart_utilities, only: read_dist_array

    integer, intent(in) :: unit, version


    call read_dist_array (unit, Fluxing_Velocity, pcell, 'READ_FLOW_DATA: error reading Fluxing_Velocity records')

  end subroutine read_flow_data

END MODULE FLUID_DATA_MODULE
