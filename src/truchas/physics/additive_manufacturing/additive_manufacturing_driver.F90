!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module additive_manufacturing_driver

  use kinds, only: r8
  use parameter_list_type
  use truchas_logging_services
  use add_manuf_type
  use additive_manufacturing_data, only: params_am 
  implicit none

  logical, private, save :: am_initialized = .false.
  type(add_manuf), private, save :: am

  public :: deposit_energy_and_mass, compute_add_manuf_source_for_ds
  
  contains


    subroutine deposit_energy_and_mass()

      use parameter_module, only: nmat
      use legacy_mesh_api,  only: ncells
      use zone_module,      only: Zone
      use matl_utilities,       only: MATL_GET_VOF, MATL_SET_VOF
      use time_step_module,     only: dt

      integer :: status
      real(r8), allocatable :: Vof(:,:)

      allocate( Vof(nmat,ncells), &
                STAT = status )
      if (status /= 0) call TLS_fatal ('deposit_energy_and_mass: allocation failed')

      ! Grab a copy of the volume fractions from Matl
      call MATL_GET_VOF (Vof)

      if (.not. am_initialized) then
        call am%init(params_am, Vof, Zone%Temp, Zone%Rho, &
                     Zone%Vc(1), Zone%Vc(2), Zone%Vc(3))
        am_initialized = .true.
      else
        call am%set_state_variables(Vof, Zone%Temp, Zone%Rho, &
                                    Zone%Vc(1), Zone%Vc(2), Zone%Vc(3))
      end if

      !!!!! REMOVE !!!!!
    !  print *, "max temp: ", maxval(am%state%temps)
    !  print *, "min temp: ", minval(am%state%temps)
      !!!!! REMOVE !!!!!

      ! save enthalpy values in order to form a source term; temperature
      ! values also need to be saved for the energy fixup routine
      ! for the diffusion solver
      am%state%enthalpies_old = am%state%enthalpies
      am%state%temps_old = am%state%temps
      ! update state (temps, enthalpies, etc.) from advection
      call am%update_state_from_advection()
      ! deposit powder to the melt pool surface
      call am%mass_depos%deposit_mass(am%state, am%am_geom, dt)
      ! update enthalpies from mass deposition (to be consistent 
      ! with the updated volume fractions)
      call am%update_state_from_mass_depos()
      ! deposit energy at free interface from laser source
      call am%energy_depos%deposit_energy(am%state, am%am_geom, dt)
      ! advance the laser-powder nozzle system one time step
      call am%am_geom%advance_am_nozzle(dt)

      ! save old void volume fractions
      am%state%void_vol_fracs_old = am%state%void_vol_fracs

      ! Update global Truchas state
      Zone%Temp = am%state%temps
      Zone%Rho = am%state%densities
      Zone%Vc(1) = am%state%velocities_x
      Zone%Vc(2) = am%state%velocities_y
      Zone%Vc(3) = am%state%velocities_z
      ! Pass the new volume fractions back into the material structure 
      call MATL_SET_VOF (am%state%vol_fracs)

      !!!! REMOVE !!!!!
      print *, "nozzle y: ", am%am_geom%am_coords(2)
      !!!! REMOVE !!!!!

      deallocate(Vof)


    end subroutine deposit_energy_and_mass


    subroutine compute_add_manuf_source_for_ds(ncell_onP, cell_ip, dQ_ds, Tmin, Tmax)

      use legacy_mesh_api, only: Cell
      use index_partitioning, only: gather_boundary, ip_desc
      use pgslib_module,        only: PGSLIB_GLOBAL_MAXVAL
    
      integer, intent(inout) :: ncell_onP
      type(ip_desc), intent(in) :: cell_ip
      real(r8), intent(out) :: dQ_ds(:), Tmin(:), Tmax(:)

      integer :: n, ncells, status
      real(r8) :: max_temp, max_temp_global, one_tol, updated_enthalpy
      real(r8), allocatable :: dQ_t(:), Tmin_t(:), Tmax_t(:)

      ncells = am%state%ncells
      one_tol = am%state%one_tol

      allocate(dQ_t(ncells), Tmin_t(ncells), Tmax_t(ncells), &
                STAT = status )
      if (status /= 0) call TLS_fatal ('compute_add_manuf_source_for_ds: allocation failed')

      ! Get the maximum temperature (global)
      max_temp = 0
      do n = 1, ncells
        if ( am%state%void_vol_fracs(n) > (1 - one_tol) ) cycle
        max_temp = max(max_temp, am%state%temps(n))
      end do
      max_temp_global = PGSLIB_GLOBAL_MAXVAL(max_temp)

      ! Set Tmin_t and Tmax_t, which is used in the diffusion solver 
      ! (this is an estimate of max/min. temp. change
      ! from advection)
      ! TODO: CHANGE THE WAY TEMPERATURE BOUNDS ARE COMPUTED
      do n = 1, ncells
        if ( am%state%temps(n) <= 0.0_r8 ) then
          Tmin_t(n) = 0
          Tmax_t(n) = max_temp_global
        else
          Tmin_t(n) = (50.0_r8/100) * am%state%temps(n)
          Tmax_t(n) = (100.0_r8/50) * am%state%temps(n)
        end if
      end do

      ! Compute the source term from the change in the enthalpy density
      ! due to fluid advection, mass deposition, and energy deposition
      dQ_t = 0
      do n = 1, ncells
        if ( am%state%void_vol_fracs(n) > (1 - one_tol) ) cycle
        ! If a cell had a negative temperature from the diffusion solve, ignore
        if ( am%state%ignore_cell(n) ) cycle
        dQ_t(n) = Cell(n)%Volume * (am%state%enthalpies(n) - am%state%enthalpies_old(n))
        ! Check that the change in energy does not result in a negative value
        updated_enthalpy = am%state%enthalpies_old(n) + dQ_t(n)/Cell(n)%Volume
      end do

      ! Put into arrays consistent with the diffusion ordering
      dQ_ds(:ncell_onP) = dQ_t(:ncell_onP)
      Tmin(:ncell_onP) = Tmin_t(:ncell_onP)
      Tmax(:ncell_onP) = Tmax_t(:ncell_onP)

      deallocate(dQ_t, Tmin_t, Tmax_t)


    end subroutine compute_add_manuf_source_for_ds




end module additive_manufacturing_driver
