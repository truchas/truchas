!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

MODULE SETUP_MODULE
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures necessary to set up the problem.
  !
  !   Public Interface(s):
  !
  !     * call SETUP ()
  !
  !       Sets up the problem.
  !
  ! Contains: SETUP
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !            Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !=======================================================================
  implicit none

  ! Private module
  private

  ! public procedures
  public :: SETUP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE SETUP ()
    !=======================================================================
    ! Purpose:
    !
    !    set up the problem
    !=======================================================================
    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use zone_module, only: zone_init
    use legacy_matl_api, only: matl_init
    use restart_variables,      only: restart, ignore_t, ignore_dt, restart_t, restart_dt, &
                                      restart_cycle_number
    use restart_driver,         only: close_restart_file, skip_restart_mesh
    use time_step_module,       only: cycle_number, cycle_number_restart, &
                                      dt, t, t1, t2
    use init_module,            only: INITIAL
    use unstr_mesh_type
    use mesh_manager,           only: init_mesh_manager, unstr_mesh_ptr
    use truchas_danu_output,    only: TDO_write_mesh
    use truchas_timers

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    type(unstr_mesh), pointer :: mesh

    ! Start the initialization timer.
    call start_timer("Initialization")

    ! Skip over unused legacy mesh data in the restart file (TEMPORARY)
    if (restart) call skip_restart_mesh

    ! Instantiate the mesh objects.
    call init_mesh_manager

    ! Allocate the base types and set them to their defaults.
    mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(mesh))
    call zone_init(mesh%ncell_onP)
    call matl_init(mesh%ncell_onP)

    ! Write the primary truchas mesh.
    call TDO_write_mesh(mesh)

    ! NNC, Sep 2014.  This used to be done after the call to INITIAL, but is
    ! needed here because of time-dependent boundary conditions.  I left the
    ! rest (DT, CYCLE_NUMBER) where it was because I do not know the
    ! consequences of moving it here, especially CYCLE_NUMBER  as there are
    ! some "if cycle_number == 0" hacks in the flow code.
    if (restart) then
      if (.not.ignore_t)  t  = restart_t
      if (.not.ignore_dt) dt = restart_dt
    end if

    ! Initialize cell-centered fluid variables and thermodynamic quantities.
    call INITIAL (t, dt)

    ! Set the initial timestep value. If this is a restart, take what is in
    ! restart file; if not, then it was already set by the input file.
    if (restart) then
      if (.not.ignore_t) cycle_number = restart_cycle_number
    end if

    if (restart) call close_restart_file ()

    ! Set previous time, last time, and initial cycle number.
    t1 = MAX(t - dt, 0.0_r8)
    t2 = t
    cycle_number_restart = cycle_number

    ! Stop the initialization timer.
    call stop_timer("Initialization")

  END SUBROUTINE SETUP

END MODULE SETUP_MODULE
