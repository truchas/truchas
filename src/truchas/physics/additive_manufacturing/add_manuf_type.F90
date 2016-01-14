!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module defines the main additive manufacturing type ADD_MANUF.
!! The type add_manuf_type contains an object, state, of type AM_STATE
!! (which is responsible for holding and manipulating state data such as temperatures,
!! enthalpies, etc.), as well as objects of abstract type AM_COORD_SYSTEM (for
!! representing and advancing the laser-powder nozzle system), INTERFACE_ENERGY
!! (for depositing energy at the material-void interface), and INTERFACE_MASS
!! (for depositing  mass at the material-void interface)
!!
!! PROGRAMMING INTERFACE
!!
!!  The derived type ADD_MANUF has the following type bound procedures:
!!
!!     INIT (THIS, PARAMS_AM, VOL_FRACS, TEMPS, &
!!           DENSITIES, VELOCITIES_X, &
!!           VELOCITIES_Y, VELOCITIES_Z) 
!!      Initializes the ADD_MANUF object THIS from parameters read in from the input file 
!!      (encapsulated in parmas_am), as well as global state data for 
!!      volume fractions, temperatures, densities, and velocities. Ultimately, the state data
!!      stored in the object STATE of type AM_STATE gets initialized appropriately
!!
!!     SET_STATE_VARIABLES (THIS, VOL_FRACS, TEMPS, &
!!                          DENSITIES, VELOCITIES_X, &
!!                          VELOCITIES_Y, VELOCITIES_Z) 
!!      Updates the ADD_MANUF object THIS from global state data for 
!!      volume fractions, temperatures, densities, and velocities (the state data
!!      stored in the object STATE of type AM_STATE gets set appropriately)
!!
!!    UPDATE_STATE_FROM_ADVECTION (THIS)
!!      Updates the state data (temperatures, enthalpies, etc.) from the advection step
!!      (note that the advection step is a separate physics component of truchas and 
!!      uses the volume-of-fluid method). The subroutine computes the energy increment, dQ,
!!      that results fluid being fluxed into and out of cells from the advection step; the
!!      change in energy density for each cell is computed from dQ and the total enthalpy
!!      density and temperature of each cell is adjusted accordingly. Note that, for cells
!!      with tiny amounts of fluid, the energy increment (computed from COMPUTE_ADVECTED_ENTHALPY,
!!      which is part of the advection physics module) can, in some cases, be such that the updated energy
!!      density is negative or non-physical; this special care is taken to fixup such 
!!      pathological cases by setting the temperature of such cases to the average
!!      temperature of its neighboring cells; since the non-void volume fraction such cells 
!!      is invariably very small, the hope is that the (arbitrary) manner in which this energy fixup 
!!      is performed should not impact the final solution too much.
!! 
!!    UPDATE_STATE_FROM_MASS_DEPOS (THIS)
!!      Updates the temperatures and enthalpies (stored in the object STATE) from the mass deposition step.
!!      (changes in temperatures and enthalpies can result from volume fraction changes from depositing mass
!!      into cells).



module add_manuf_type

  use kinds, only: r8
  use parameter_list_type
  use truchas_logging_services

  use am_state_type
  use am_coord_system_class
  use interface_energy_class
  use interface_mass_class

  use am_coord_system_factory, only : alloc_am_coord_system
  use interface_energy_factory, only : alloc_interface_energy
  use interface_mass_factory, only:   alloc_interface_mass

  implicit none

  type, public :: add_manuf
    type(am_state) :: state
    class(am_coord_system), allocatable :: am_geom
    class(interface_energy), allocatable :: energy_depos
    class(interface_mass), allocatable :: mass_depos
  contains
    procedure :: init => init_add_manuf
    procedure :: set_state_variables
    procedure :: update_state_from_advection
    procedure :: update_state_from_mass_depos
  end type add_manuf


  contains


    subroutine init_add_manuf(this, params_am, vol_fracs, temps, &
                              densities, velocities_x, &
                              velocities_y, velocities_z)
  
      class(add_manuf), intent(inout) :: this
      type(parameter_list) :: params_am
      real(r8), intent(in) :: vol_fracs(:,:), temps(:), &
                              densities(:), velocities_x(:), &
                              velocities_y(:), velocities_z(:)

      type(parameter_list), pointer :: p_params



      ! allocate and initialize am_state object
      call this%state%alloc_and_init(vol_fracs, temps, &
                                     densities, velocities_x, &
                                     velocities_y, velocities_z)

      ! allocate and initialize the laser-powder nozzle object
      if ( params_am%is_sublist('am_coord_system')) then
        p_params => params_am%sublist('am_coord_system')
        call alloc_am_coord_system(this%am_geom, p_params)
        call this%am_geom%init(p_params)
      else
        call TLS_fatal ('init_add_manuf: params list for am_coord_system is not initialized')
      end if


      ! allocate and initialize the interface energy object
      if ( params_am%is_sublist('interface_energy')) then
        p_params => params_am%sublist('interface_energy')
        call alloc_interface_energy(this%energy_depos, p_params)
        call this%energy_depos%init(p_params)
      else
        call TLS_fatal ('init_add_manuf: params list for interface_energy is not initialized')
      end if

      ! allocate and initialize the interface mass object
      if ( params_am%is_sublist('interface_mass')) then
        p_params => params_am%sublist('interface_mass')
        call alloc_interface_mass(this%mass_depos, p_params)
        call this%mass_depos%init(p_params)
      else
        call TLS_fatal ('init_add_manuf: params list for interface_mass is not initialized')
      end if


    end subroutine init_add_manuf



    subroutine set_state_variables(this, vol_fracs, temps, &
                                   densities, velocities_x, &
                                   velocities_y, velocities_z)
  
      class(add_manuf), intent(inout) :: this
      real(r8), intent(in) :: vol_fracs(:,:), temps(:), &
                              densities(:), velocities_x(:), &
                              velocities_y(:), velocities_z(:)

      call this%state%set_state(vol_fracs, temps, &
                                densities, velocities_x, &
                                velocities_y, velocities_z)

    end subroutine set_state_variables




    subroutine update_state_from_advection(this)
 
      use mesh_module, only: Cell
      use advection_module, only: compute_advected_enthalpy

      class(add_manuf), intent(inout) :: this

      integer :: status, n, ncells
      real(r8) :: T, H, one_tol, am_tol, void_frac, void_frac_old
      real(r8), allocatable :: dQ(:)
      logical, allocatable :: is_a_fixup_cell(:)

      one_tol = this%state%one_tol
      am_tol = this%state%am_tol
      ncells = this%state%ncells


      allocate( dQ(ncells), is_a_fixup_cell(ncells), &
                STAT = status )
      if (status /= 0) call TLS_panic ('update_sate_from_advection: allocation failed')

 
      ! Compute the total enthalpy change, dQ, from advection
      call compute_advected_enthalpy (this%state%temps, dQ)
      
      is_a_fixup_cell = .false.
      do n = 1, ncells
        void_frac = this%state%void_vol_fracs(n)
        void_frac_old = this%state%void_vol_fracs_old(n)
        ! Don't deposit energy into void cells
        ! This is to avoid issues with depositing energy into cells with 
        if ( void_frac > 1.0_r8 - one_tol ) cycle
        ! Update formerly void cells; here it is assumed that the energy of the cell
        ! is nonzero only if energy has fluxed in from the advection step
        if ( void_frac_old > (1 - one_tol) ) then
          ! This is a formerly void cell and previously had no enthalpy
          this%state%enthalpies_old(n) = 0
          ! The cell's current enthalpy comes from what was advected into it
          this%state%enthalpies(n) = dQ(n)/Cell(n)%Volume
        else ! the cell was not void on the previous time step
          ! Update the corresponding enthalpy density
          this%state%enthalpies(n) = this%state%enthalpies(n) + dQ(n)/Cell(n)%Volume
        end if

        ! Update the corresponding temperature to be consistent with the enthalpy density,
        ! if the cell is not nearly void (NOTE: can be nearly void and totally void on the
        ! previous time step, in which case need to set the enthalpy density via a call to
        ! "fixup_energy_in_nearly_void_cells" below
        H = this%state%enthalpies(n)
        if ( void_frac > (1 - am_tol) .or. H < 0  ) then
          is_a_fixup_cell(n) = .true.
        else
          ! Update the temperature to be consistent with the new enthalpy density
          call this%state%compute_temp(n, H, T)
          if ( this%state%root_solver%found_root ) then ! can save the computed enthalpy density
            this%state%temps(n) = T
            is_a_fixup_cell(n) = .false.
          else
            is_a_fixup_cell(n) = .true.
          end if
        end if

      end do
 
      ! Set the temperature and enthalpy to be that of its
      ! neighbor with the smallest amount of void
      call this%state%fixup_energy_in_nearly_void_cells(is_a_fixup_cell)

      deallocate(dQ, is_a_fixup_cell)

     
    end subroutine update_state_from_advection


    subroutine update_state_from_mass_depos(this)

      class(add_manuf), intent(inout) :: this

      integer :: n, ncells
      real(r8) :: H, one_tol

      ! Compute the enthalpy density change corresponding to changes 
      ! in material volume fractions from mass deposition and redistribution;
      ! note that mass redistribution occurs only if a cell overflows

      one_tol = this%state%one_tol
      ncells = this%state%ncells

      do n = 1, ncells
        if ( this%state%void_vol_fracs(n) > (1 - one_tol) .or. this%state%ignore_cell(n)  ) cycle
        ! Update enthalpy from mass deposition
      !  if ( this%state%temps(n) <= 0 ) then
      !    print *, 'update_enthalpy_density_and_temp_from_mass_depos: non-positive temp'
      !    call TLS_fatal ('update_enthalpy_density_and_temp_from_mass_depos: non-positive temp')
      !  end if
        call this%state%compute_enthalpy_density(n, H)
        this%state%enthalpies(n) = H
        ! The volume fraction in a cell may have changed due to 
        ! the redistribution stage, where fluid "overflowed" from a neighboring
        ! cell. In order to form a source for the diffusion solver that 
        ! captures the associated enthalpy change, we set the old enthalpy of
        ! the cell to be zero
        if ( this%state%void_vol_fracs_old(n) > (1- one_tol) ) then
          this%state%enthalpies_old(n) = 0
        end if

      end do    

    end subroutine update_state_from_mass_depos






end module add_manuf_type
