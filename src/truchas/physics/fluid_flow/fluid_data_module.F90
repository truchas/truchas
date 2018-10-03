!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use kinds, only: r8
  use parameter_module, only: maxmat
  use legacy_mesh_api, only: ndim
  implicit none
  private

  public read_flow_data

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Flag for enabling/disabling fluid flow.
  logical,  public, save :: fluid_flow
  logical,  public, save :: fluid_to_move
  logical,  public, save :: mass_limiter
  real(r8), public, save :: mass_limiter_cutoff

  ! Flag for applying a prescribed velocity field.
  logical, public, save :: applyflow

  ! Flag that defines the existance of a void material in this calculation and its index
  logical, public, save :: Void_Material_Exists
  integer, public, save, dimension(maxmat) :: Void_Material_Index
  integer, public, save :: Void_material_Count

  ! Value of pressure in the void.
  real(r8), public, save :: void_pressure

  ! Boussinesq approximation
  logical, public, save :: boussinesq_approximation

  ! Total volume flowing into and out of the domain.
  real(r8), public, save :: qin = 0.0_r8, qout = 0.0_r8

  ! Volume Fraction Cutoff Value  (no flow solution if Vof < fluid_cutoff)
  real(r8), public, parameter :: fluid_cutoff = 0.01

  ! Time-weighting for treatment of momentum deposition due to solidification
  real(r8), public, save :: momentum_solidify_implicitness

  ! Flow arrays.
  real(r8), pointer, public, save, dimension(:,:) :: Fluxing_Velocity => null()
  real(r8), pointer, public, save, dimension(:,:) :: Face_Interpolation_Factor => null()
  real(r8), pointer, public, save, dimension(:,:) :: Momentum_by_Volume => null()
  real(r8), pointer, public, save, dimension(:)   :: fluidVof => null(), fluidVof_n => null()
  real(r8), pointer, public, save, dimension(:)   :: fluidRho => null()
  real(r8), pointer, public, save, dimension(:)   :: fluidRho_n => null()
  real(r8), pointer, public, save, dimension(:)   :: fluidDeltaRho => null()
  real(r8), pointer, public, save, dimension(:)   :: realfluidVof => null()
  real(r8), pointer, public, save, dimension(:)   :: cutRho => null()

! Volume-fraction averaged density at n, n+1 for body force terms
  real(r8), pointer, public, save, dimension(:) :: avgRho => null()
  real(r8), pointer, public, save, dimension(:) :: avgRho_n => null()

  real(r8), allocatable, public, save, dimension(:,:) :: Drag_Coefficient
  real(r8), pointer, public, save, dimension(:)       :: courant => null()

  ! Arrays related to partially/totally solidified cells.
  logical, pointer, public, save, dimension(:)   :: isPureImmobile => null()
  logical, pointer, public, save, dimension(:,:) :: Solid_Face => null()
 
  ! Non-Void Cell list
  logical, pointer, public, save, dimension(:)   :: Cell_isnt_Void => null()
  logical, pointer, public, save, dimension(:,:) :: Ngbr_isnt_Void => null()

  ! Projection / Predictor Communication Arrays
  real(r8), public, allocatable, save, dimension(:,:) :: Centered_GradP_Dynamic
  real(r8), public, allocatable, save, dimension(:,:) :: Rho_Face, Rho_Face_n
  real(r8), public, save                              :: MinFluidRho, MinFaceFraction

  ! Storage moved to a single-allocation per run, rather than once per time-step (MAC)
  real(r8), public, allocatable, save, dimension(:,:) :: Mom_Delta
  real(r8), public, allocatable, save, dimension(:)   :: pc_delta_rho
  real(r8), public, save, dimension(ndim) :: minVel, maxVel   

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

    use legacy_mesh_api, only: pcell => unpermute_mesh_vector
    use restart_utilities, only: read_dist_array

    integer, intent(in) :: unit, version

    call read_dist_array (unit, Fluxing_Velocity, pcell, 'READ_FLOW_DATA: error reading Fluxing_Velocity records')

  end subroutine read_flow_data

END MODULE FLUID_DATA_MODULE
