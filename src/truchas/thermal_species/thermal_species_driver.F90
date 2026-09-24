!!
!! THERMAL_SPECIES_DRIVER
!!
!! This module provides the legacy package-level procedural interface for the
!! thermal/species transport physics. It translates namelist input into the
!! parameter-list structure used by the refactored solvers, owns the package
!! object, and forwards cycle-driver and coupling calls to that object.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module thermal_species_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use thermal_species_physics_type
  use parameter_list_type
  use rad_problem_type, only: vf_event
  implicit none
  private

  public :: thermal_species_init, thermal_species_step, thermal_species_commit_step
  public :: thermal_species_read_namelists
  public :: thermal_species_enabled
  public :: thermal_species_mesh_name
  public :: thermal_species_set_initial_state, thermal_species_restart
  public :: thermal_species_delete
  public :: thermal_species_update_moving_vf, thermal_species_add_moving_vf_events, vf_event
  public :: thermal_species_have_heat_transport, thermal_species_have_species_transport, thermal_species_num_species

  !! These return cell-centered field values.
  public :: thermal_species_get_temp, thermal_species_get_enthalpy, thermal_species_get_phi
  public :: thermal_species_get_temp_grad

  !! Face temperature is a solver coupling variable, not a general physics field.
  public :: thermal_species_get_face_temp
  public :: thermal_species_get_face_temp_view

  !! TARGET is needed because solver/model/property objects retain references
  !! to package-owned components of this object.
  type(thermal_species_physics), target :: physics
  type(parameter_list), pointer :: thermal_species_params => null() ! owned object

  integer, parameter :: TS_SPECIES_SYS = 1
  integer, parameter :: TS_THERMAL_SYS = 2
  integer, parameter :: TS_THERMAL_SPECIES_SYS = 3
  character(*), parameter :: mesh_name = 'MAIN'

  logical :: package_initialized = .false.

contains

  logical function thermal_species_enabled()
    thermal_species_enabled = associated(thermal_species_params) .or. package_initialized
  end function thermal_species_enabled

  character(4) function thermal_species_mesh_name()
    thermal_species_mesh_name = mesh_name
  end function thermal_species_mesh_name

  subroutine thermal_species_read_namelists(lun, heat_transport, species_transport, number_of_species)
    integer, intent(in) :: lun
    logical, intent(in) :: heat_transport, species_transport
    integer, intent(in) :: number_of_species
    call read_params(lun, heat_transport, species_transport, number_of_species)
  end subroutine

  subroutine read_params(lun, heat_transport, species_transport, number_of_species)

    use diffusion_solver_namelist
    use physical_constants, only: stefan_boltzmann, absolute_zero

    integer, intent(in) :: lun
    logical, intent(in) :: heat_transport, species_transport
    integer, intent(in) :: number_of_species
    type(parameter_list), pointer :: solvers_params, ht_params, sd_params, htsd_params
    type(parameter_list), pointer :: heat_model_params, species_model_params
    type(parameter_list), pointer :: htsd_model_params, solver_params
    integer :: system_type
    real(r8) :: void_temperature

    if (associated(thermal_species_params)) deallocate(thermal_species_params)
    allocate(thermal_species_params)
    solvers_params => thermal_species_params%sublist('solvers')

    system_type = thermal_species_system_type(heat_transport, species_transport, number_of_species)
    select case (system_type)
    case (TS_SPECIES_SYS)
      call thermal_species_params%set('declared-solver', 'sd')
      sd_params => solvers_params%sublist('species')
      solver_params => sd_params%sublist('solver')
      species_model_params => sd_params%sublist('model')
    case (TS_THERMAL_SYS)
      call thermal_species_params%set('declared-solver', 'ht')
      ht_params => solvers_params%sublist('heat')
      solver_params => ht_params%sublist('solver')
      heat_model_params => ht_params%sublist('model')
    case (TS_THERMAL_SPECIES_SYS)
      call thermal_species_params%set('declared-solver', 'htsd')
      htsd_params => solvers_params%sublist('htsd')
      solver_params => htsd_params%sublist('solver')
      htsd_model_params => htsd_params%sublist('model')
      heat_model_params => htsd_model_params%sublist('heat')
      species_model_params => htsd_model_params%sublist('species')
    case default
      INSIST(.false.)
    end select

    call read_diffusion_solver_namelist(lun, system_type, number_of_species, solver_params)

    if (heat_transport) then
      call heat_model_params%set('stefan-boltzmann', stefan_boltzmann)
      call heat_model_params%set('absolute-zero', absolute_zero)
      call solver_params%get('void-temperature', void_temperature)
      call heat_model_params%set('void-temperature', void_temperature)
      call read_heat_model_namelists(lun, heat_model_params, solver_params)
    end if

    if (species_transport) then
      call species_model_params%set('num-species', number_of_species)
      call read_species_model_namelists(lun, species_model_params)
    end if

  end subroutine read_params

  integer function thermal_species_system_type(heat_transport, species_transport, number_of_species)
    logical, intent(in) :: heat_transport, species_transport
    integer, intent(in) :: number_of_species

    INSIST(heat_transport .or. species_transport)
    if (species_transport) then
      INSIST(number_of_species > 0)
    end if

    if (heat_transport .and. species_transport) then
      thermal_species_system_type = TS_THERMAL_SPECIES_SYS
    else if (heat_transport) then
      thermal_species_system_type = TS_THERMAL_SYS
    else
      thermal_species_system_type = TS_SPECIES_SYS
    end if
  end function thermal_species_system_type

  subroutine read_heat_model_namelists(lun, params, solver_params)

    use enclosure_radiation_namelist, only: read_enclosure_radiation_namelists
    use thermal_bc_namelist
    use thermal_source_namelist

    integer, intent(in) :: lun
    type(parameter_list), pointer, intent(in) :: params
    type(parameter_list), pointer, intent(in) :: solver_params

    integer :: vfr_precon_iter
    real(r8) :: vfr_solve_tol
    character(:), allocatable :: vfr_precon_method
    type(parameter_list), pointer :: bc_params, source_params, rad_params, solver_rad_params

    bc_params => params%sublist('bc')
    call read_thermal_bc_namelists(lun, bc_params)

    source_params => params%sublist('sources')
    call read_thermal_source_namelists(lun, source_params)

    rad_params => params%sublist('radiation')
    call read_enclosure_radiation_namelists(lun, rad_params)
    solver_rad_params => solver_params%sublist('radiation')
    call solver_rad_params%get('solve-tol', vfr_solve_tol)
    call solver_rad_params%get('precon-method', vfr_precon_method)
    call solver_rad_params%get('precon-iter', vfr_precon_iter)
    call rad_params%set('solve-tol', vfr_solve_tol)
    call rad_params%set('precon-method', vfr_precon_method)
    call rad_params%set('precon-iter', vfr_precon_iter)

  end subroutine read_heat_model_namelists

  subroutine read_species_model_namelists(lun, params)

    use species_bc_namelist
    use species_source_namelist

    integer, intent(in) :: lun
    type(parameter_list), pointer, intent(in) :: params

    type(parameter_list), pointer :: bc_params, source_params

    bc_params => params%sublist('bc')
    call read_species_bc_namelists(lun, bc_params)

    source_params => params%sublist('sources')
    call read_species_source_namelists(lun, source_params)

  end subroutine read_species_model_namelists

  subroutine thermal_species_step(t, hnext, errc)
    real(r8), intent(in)  :: t
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: errc
    call physics%step(t, hnext, errc)
  end subroutine

  subroutine thermal_species_commit_step
    call physics%commit_step
  end subroutine

  logical function thermal_species_have_heat_transport()
    thermal_species_have_heat_transport = physics%has_heat_transport()
  end function thermal_species_have_heat_transport

  logical function thermal_species_have_species_transport()
    thermal_species_have_species_transport = physics%has_species_transport()
  end function thermal_species_have_species_transport

  integer function thermal_species_num_species()
    thermal_species_num_species = physics%get_num_species()
  end function thermal_species_num_species

  subroutine thermal_species_get_temp (array)
    real(r8), intent(inout) :: array(:)
    call physics%get_temp(array)
  end subroutine

  subroutine thermal_species_get_enthalpy (array)
    real(r8), intent(inout) :: array(:)
    call physics%get_enthalpy(array)
  end subroutine

  subroutine thermal_species_get_phi(n, array)
    integer,  intent(in)  :: n
    real(r8), intent(inout) :: array(:)
    call physics%get_phi(n, array)
  end subroutine

  subroutine thermal_species_get_temp_grad(array)
    real(r8), intent(inout) :: array(:,:)
    call physics%get_temp_grad(array)
  end subroutine

  subroutine thermal_species_get_face_temp(array)
    real(r8), intent(inout) :: array(:)
    call physics%get_face_temp(array)
  end subroutine

  subroutine thermal_species_get_face_temp_view (view)
    real(r8), pointer :: view(:)
    call physics%get_face_temp_view(view)
  end subroutine

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! THERMAL_SPECIES_INIT
 !!
  subroutine thermal_species_init
    use mesh_manager, only: unstr_mesh_ptr
    use unstr_mesh_type

    type(unstr_mesh), pointer :: mesh

    mesh => unstr_mesh_ptr(mesh_name)
    INSIST(associated(mesh))
    INSIST(associated(thermal_species_params))
    call physics%init(mesh, thermal_species_params)
    package_initialized = .true.
  end subroutine

  subroutine thermal_species_delete
    if (allocated(physics%solver)) deallocate(physics%solver)
    physics%mesh => null()
    package_initialized = .false.
    if (associated(thermal_species_params)) then
      deallocate(thermal_species_params)
      thermal_species_params => null()
    end if
  end subroutine

  subroutine thermal_species_set_initial_state (t, dt, temp, conc)
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    call physics%set_initial_state(t, dt, temp, conc)
  end subroutine

  !! The effect of calling this subroutine is to restart or reset the solver so
  !! that its subsequent numerical behavior is as if it was starting integration
  !! from an initial state equal to the current state.  This mainly means
  !! dropping any previous solution history in the BDF2 solver and recomputing
  !! an approximation to the initial state time derivative.  The argument DT is
  !! used to compute that time derivative and is best chosen equal to the next
  !! time step size; however it has no bearing on the next step size used.

  subroutine thermal_species_restart(dt)
    real(r8), intent(in) :: dt
    call physics%restart(dt)
  end subroutine

  subroutine thermal_species_update_moving_vf
    call physics%update_moving_vf
  end subroutine

  subroutine thermal_species_add_moving_vf_events(eventq, rank)
    use sim_event_queue_type
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank
    call physics%add_moving_vf_events(eventq, rank)
  end subroutine

end module thermal_species_driver
