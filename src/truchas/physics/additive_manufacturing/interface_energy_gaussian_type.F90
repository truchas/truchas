!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module defines the derived type INTERFACE_ENERGY_GAUSSIAN, which
!! inherits from the abstract type INTERFACE_ENERGY. Objects of type 
!! INTERFACE_ENERGY_GAUSSIAN model the energy deposition process using an
!! energy source that has a Gaussian profile; energy is deposited at the
!! material-void interface, which consists of cells that have no void
!! and are under cells with some void (this avoids the delicate issues
!! associated with depositing energy in cells with tiny fluid fragments,
!! which can result in unphysical temperature values).
!!
!! PROGRAMMING INTERFACE
!!
!!  The class INTERFACE_ENERGY has the following type bound procedures:
!!
!!     ALLOC (THIS) 
!!      Allocates THIS%interface_cells(ncells)
!!
!!     INIT (THIS, PARAMS) 
!!      Initializes THIS%laser_absorp, THIS%laser_powder, etc., from parameters read in from the input file
!!      (encapsulated in the object parmas of type PARAMETER_LIST).
!!
!!   DEPOSIT_ENERGY (THIS, STATE, AM_GEOM, DT)
!!    Updates the state data (e.g. temperatures and enthalpies) stored in the object STATE 
!!    (of type AM_STATE). This procedure uses as input the object AM_GEOM
!!    (of some type derived from class AM_COORD_SYSTEM) and the time step DT
!!
!!   T_AFTER_ENERGY_DEPOS(THIS, STATE, AM_GEOM, NCELL, DT)
!!    Returns the updated temperature in cell number NCELL after the energy deposition process takes
!!    place over a time step of DT.
!!
!!   LOCATE_INTERFACE_CELLS (THIS, STATE)
!!    Finds all cells that have no void and are directly underneath a cell with at least some void. (Initially,
!!    interface cells were defined as those with at least some void, but depositing energy in such cells
!!    led to numerical difficulties if the cell was mostly void). This assumes that the mesh cells
!!    are hexagons. 
!!

module interface_energy_gaussian_type

  use kinds, only: r8
  use parameter_list_type
  use interface_energy_class
  use am_coord_system_class
  use am_state_type
  implicit none

  public :: print_obj, energy_source_gaussian
  private

  type, extends(interface_energy), public :: interface_energy_gaussian
    real(r8) :: stefan_boltzmann, abszero
    real(r8) :: laser_absorp, laser_sigma, laser_power
    real(r8) :: vof_cutoff
    real(r8) :: source_params(9)
    logical, allocatable :: interface_cells(:)
    type(root_finder) :: root_solver
    ! root_solver is used to solve for H(T) - H_0 = dt  * (Q_source  - Q_sink(T)), 
    ! where T = temp and H = enthalpy density
    contains
       procedure :: init => init_energy_depos
       procedure :: alloc => alloc_energy_depos
       procedure :: deposit_energy
       procedure :: T_after_energy_depos => T_after_energy_depos
       procedure :: locate_interface_cells => locate_interface_cells
       procedure :: print_obj => print_obj
  end type interface_energy_gaussian


  contains


    subroutine alloc_energy_depos(this)

      use parameter_module, only: ncells

      class(interface_energy_gaussian), intent(inout) :: this
 
      integer :: status

      ! Allocate array to hold fraction of heat to be deposited in cell

      allocate( this%interface_cells(ncells), STAT = status )

      if (status /= 0) call TLS_fatal ('alloc_gaussian: allocation failed')      

    end subroutine alloc_energy_depos



    subroutine init_energy_depos(this, params)

      use parameter_list_type

      class(interface_energy_gaussian), intent(inout) :: this
      type(parameter_list), intent(inout) :: params

      call params%get('laser_absorp', this%laser_absorp)
      call params%get('laser_power', this%laser_power)
      call params%get('laser_sigma', this%laser_sigma)
      call params%get('stefan_boltzmann', this%stefan_boltzmann)
      call params%get('abszero', this%abszero)

      this%source_params(1) = this%laser_absorp
      this%source_params(2) = this%laser_power
      this%source_params(3) = this%laser_sigma
    !  this%source_params(4) = this%laser_depth
      this%source_params(5) = this%stefan_boltzmann ! Stefan-Boltzmann coeffient
      this%source_params(6) = this%abszero ! absolute zero
      !!!!!! CHANGE !!!!!!
      this%vof_cutoff = 10.0_r8**(-12)
      !!!!!! CHANGE !!!!!!

      ! Initialize the root solver object (this is used for
      ! computing the enthalpy density from the temperature and 
      ! volume fractions
      !!!!!! TODO: CHANGE---READ THIS IN FROM INPUT !!!!!!!!
      call this%root_solver%init(100, 100, 100, 10.0d0**(-12))
      !!!!!! CHANGE !!!!!!!!    

    end subroutine init_energy_depos



    subroutine locate_interface_cells(this, state)

      use gs_module,     only: EE_GATHER
      use mesh_module,      only: Cell

      class(interface_energy_gaussian), intent(inout) :: this
      type(am_state), intent(in) :: state

      integer :: n, f, f_up(1), ncells, nfc, status, ndim
      real(r8) :: upward_unit_vec(1:3), volume_void_neighbor, &
                  vof_cutoff
      real(r8), allocatable :: normals_dotted_with_up_vec(:), &
                               void_vol_fracs_nghbrs(:,:)

      nfc = state%nfaces
      ncells = state%ncells
      vof_cutoff = this%vof_cutoff
      ndim = 3

      allocate( void_vol_fracs_nghbrs(1:nfc,1:ncells), &
                normals_dotted_with_up_vec(1:nfc), &
                STAT = status )
      if (status /= 0) call TLS_panic ('locate_interface_cells_vof: allocation failed') 

      call EE_GATHER(void_vol_fracs_nghbrs, state%void_vol_fracs)
    
      upward_unit_vec = (/0,0,1/) 
      this%interface_cells = .false.
      do n = 1, ncells
        ! Find all cells without any void that are directly under cells with some void
        ! NOTE: THIS ASSUMES THAT THE CELLS ARE HEXES
        do f = 1, nfc
          normals_dotted_with_up_vec(f) = &
                 sum( Cell(n)%Face_Normal(1:ndim,f) * upward_unit_vec(1:ndim) )
        end do
        f_up(1:1) = maxloc( normals_dotted_with_up_vec(1:nfc) )
        ! Find the void volume of the cell above
        volume_void_neighbor = void_vol_fracs_nghbrs(f_up(1),n)
        ! If the cell has no void and the cell above has some void, then this
        ! in in the interface layer
        this%interface_cells(n) = &
              (volume_void_neighbor > state%am_tol) .and. &
              (state%void_vol_fracs(n) < vof_cutoff)
      end do

      deallocate(normals_dotted_with_up_vec, void_vol_fracs_nghbrs)


    end subroutine locate_interface_cells


    function energy_source_gaussian(absorp, power, sigma, dt, dx, dy) result(q)

      real(r8), intent(in) :: absorp, power, sigma, dt, dx, dy
      real(r8) :: q

      real(r8) :: r, two_sigma2
      real(r8), parameter :: PI = 3.141592653589793_r8

      two_sigma2 = sigma**2
      r = 2 * (dx**2 + dy**2) / two_sigma2
      q = dt * absorp * (2*power) * exp(-r) / (PI*two_sigma2)

    end function energy_source_gaussian


    function energy_residual_gaussian(T, residual_params) result(residual)

      !!!!! CHANGE !!!!!
      use parameter_module,     only: nmat
      !!!!! CHANGE !!!!!
      use material_interop,     only: ds_enthalpy_density

      real(r8), intent(in) :: T, residual_params(:) 
      real(r8) :: residual, dt, A_over_V, stefan_boltzmann, abszero

      integer :: m
      real(r8) :: state(1), H_old, vol_fracs(1:nmat), H, q

      ! extract parameters for evaluation of the residual
      q = residual_params(1)
      H_old = residual_params(2)
      vol_fracs(1:nmat) = residual_params(3:nmat+2)
      dt = residual_params(nmat+3)
      A_over_V = residual_params(nmat+4)
      stefan_boltzmann = residual_params(nmat+5)
      abszero = residual_params(nmat+6)

      ! compute the enthalpy density rho(Vof) * H(Vof,T) 
      state(1) = T
      H = 0
      do m = 1, nmat
        H = H + vol_fracs(m) * ds_enthalpy_density(m,state)
      end do

      ! remove energy due to radiation
      q = q - dt*A_over_V * stefan_boltzmann * (T**4 - abszero**4)

      residual = (H - H_old - q)/max(H_old,1.0_r8)


    end function energy_residual_gaussian



    function T_after_energy_depos(this, state, am_geom, ncell, dt) result(T_new)

      use mesh_module,          only: Cell

      class(interface_energy_gaussian), intent(inout) :: this
      type(am_state), intent(in) :: state
      class(am_coord_system), intent(in) :: am_geom
      integer, intent(in) :: ncell
      real(r8), intent(in) :: dt
      real(r8) :: T_new

      integer :: nmat, nfc, f, status, f_up(1)
      real(r8) :: q, Tmin, Tmax, T, upward_vec(1:3), A, V, dx, dy
      real(r8), allocatable :: params(:), proj_onto_up_vec(:)

      nmat = state%nmat
      nfc = state%nfaces
      allocate(params(1:nmat+6), proj_onto_up_vec(1:nfc), STAT = status )
      if (status /= 0) call TLS_fatal ('T_after_energy_depos: allocation failed')

      dx = Cell(ncell)%Centroid(1) - am_geom%am_coords(1)
      dy = Cell(ncell)%Centroid(2) - am_geom%am_coords(2)
      T = state%temps(ncell)
      q = energy_source_gaussian(this%laser_absorp, this%laser_power, this%laser_sigma, dt, dx, dy)

      ! Multiply by (cell volume)/(surface area), where (surface area) denotes
      ! the area of the surface that the energy is being fluxed through
      ! NOTE: THIS ASSUMES THAT THE CELLS ARE HEXES
      upward_vec(1:3) = (/0,0,1/)
      do f = 1, nfc
        proj_onto_up_vec(f) = &
            sum( Cell(ncell)%Face_Normal(:,f) * upward_vec(:) )
      end do
      f_up(1:1) = maxloc( proj_onto_up_vec(1:nfc) )
      A = Cell(ncell)%Face_Area(f_up(1))
      V = Cell(ncell)%Volume
      ! compute the enthalpy density change in the cell
      q = q*A/V

      if ( q < 10.0_r8**(-10) ) then
        this%root_solver%is_bracketed = .false.
        this%root_solver%found_root = .false.
        return
      end if

      ! Compute the temperature such that 
      ! energy_residual_gaussian(T_new) = 0
      params(1) = q 
      params(2) = state%enthalpies(ncell)
      params(3:nmat+2) = state%vol_fracs(1:nmat,ncell)
      params(nmat+3) = dt
      params(nmat+4) = A/V
      params(nmat+5) = this%stefan_boltzmann
      params(nmat+6) = this%abszero

      this%root_solver%found_root = .false.
      this%root_solver%is_bracketed = .false.
      T_new = T
      call this%root_solver%bracket_root(energy_residual_gaussian, params, Tmin, Tmax)
      if ( this%root_solver%is_bracketed ) then
        call this%root_solver%compute_root(energy_residual_gaussian, params, Tmin, Tmax)
        if ( this%root_solver%found_root ) then
          T_new = this%root_solver%root
        else
          print *, "T_after_energy_depos: couldn't find root: "
        end if
      else
        print *, "T_after_energy_depos: couldn't bracket root: "
      end if

      deallocate(params, proj_onto_up_vec)
      

    end function T_after_energy_depos



    subroutine deposit_energy(this, state, am_geom, dt)

      class(interface_energy_gaussian), intent(inout) :: this
      type(am_state), intent(inout) :: state
      class(am_coord_system), intent(in) :: am_geom
      real(r8), intent(in) :: dt

      integer ::n
      real(r8) :: T, H

      ! store time step for use in function "energy_source_gaussian"
     ! this%source_params(9) = dt

      ! compute the interface cells for energy deposition
      call this%locate_interface_cells(state)

      ! deposit energy in each cell
      do n = 1, state%ncells
        if ( .not. this%interface_cells(n)  ) cycle
        ! compute the temperature, T, after energy is deposited in the cell
        ! this involves solving some equation F(H(T),T, dt) = 0 in each cell
        T = this%T_after_energy_depos(state, am_geom, n, dt)
        ! if the value for T was computed correctly, update the state
        ! with the new temp. and the new corresponding enthalpy density;
        ! if the computation of the new T failed (invariably because a cell
        ! contains a tiny volume fraction of fluid on the order of roundoff),
        ! then keep the old temp. and enthalpy
        if ( this%root_solver%found_root .and. T > 0  ) then
          state%temps(n) = T
          call state%compute_enthalpy_density(n, H)
          state%enthalpies(n) = H
        end if
      end do


    end subroutine deposit_energy



    subroutine print_obj(this)

      class(interface_energy_gaussian), intent(inout) :: this

      print *, "Stefan-Boltzmann constant: ", this%stefan_boltzmann
      print *, "abszero constant: ", this%abszero
      print *, "laser_absorp: ", this%laser_absorp
      print *, "laser_sigma: ", this%laser_sigma
      print *, "laser_power: ", this%laser_power

    end subroutine print_obj






end module interface_energy_gaussian_type

