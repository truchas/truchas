!!
!! THERMAL_SPECIES_SOLVER_CLASS
!!
!! This module defines the abstract solver interface used by the package-level
!! thermal/species physics object. Concrete solvers implement the active
!! thermal, species, or coupled thermal/species transport system while exposing
!! a common stepping, restart, field-query, and external-rate interface.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module thermal_species_solver_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  !! Solver result accessors are part of the solver/physics boundary.
  !! Prefer copy accessors for stable external use. View accessors expose live
  !! solver storage and should remain limited to synchronous coupling paths.
  !! Solver internals should access their own vector components directly.
  type, abstract, public :: thermal_species_solver
  contains
    procedure(step), deferred :: step
    procedure(commit_step), deferred :: commit_step
    procedure(restart), deferred :: restart
    procedure(set_initial_state), deferred :: set_initial_state
    procedure(last_time), deferred :: last_time
    procedure(log_step_stats), deferred :: log_step_stats
    !! Physics-facing hooks guarded by THERMAL_SPECIES_PHYSICS capability flags.
    !! Default implementations indicate an internal solver/capability
    !! mismatch.
    procedure :: get_cell_temp_copy => unsupported_get_cell_temp_copy
    procedure :: get_cell_heat_copy => unsupported_get_cell_heat_copy
    procedure :: get_face_temp_copy => unsupported_get_face_temp_copy
    procedure :: get_face_temp_view => unsupported_get_face_temp_view
    procedure :: get_cell_temp_grad => unsupported_get_cell_temp_grad
    procedure :: get_cell_conc_copy => unsupported_get_cell_conc_copy
    !! Set on-process, cell-integrated external rates for solver state variables.
    procedure :: set_ext_enthalpy_rate => unsupported_set_ext_enthalpy_rate
    procedure :: set_ext_species_rate => unsupported_set_ext_species_rate
  end type

  abstract interface
    subroutine step(this, t, hnext, errc)
      import thermal_species_solver, r8
      class(thermal_species_solver), intent(inout) :: this
      real(r8), intent(in)  :: t
      real(r8), intent(out) :: hnext
      integer,  intent(out) :: errc
    end subroutine
    subroutine commit_step(this)
      import thermal_species_solver
      class(thermal_species_solver), intent(inout) :: this
    end subroutine
    subroutine restart(this, dt)
      import thermal_species_solver, r8
      class(thermal_species_solver), intent(inout) :: this
      real(r8), intent(in) :: dt
    end subroutine
    subroutine set_initial_state(this, t, dt, temp, conc)
      import thermal_species_solver, r8
      class(thermal_species_solver), intent(inout) :: this
      real(r8), intent(in) :: t, dt
      real(r8), intent(in), optional :: temp(:), conc(:,:)
    end subroutine
    function last_time(this) result(t)
      import thermal_species_solver, r8
      class(thermal_species_solver), intent(in) :: this
      real(r8) :: t
    end function
    subroutine log_step_stats(this)
      import thermal_species_solver
      class(thermal_species_solver), intent(in) :: this
    end subroutine
  end interface

contains

  subroutine unsupported_get_cell_temp_copy(this, copy)
    class(thermal_species_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    call unsupported_solver_operation(this, 'get_cell_temp_copy', size(copy))
  end subroutine

  subroutine unsupported_get_cell_heat_copy(this, copy)
    class(thermal_species_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    call unsupported_solver_operation(this, 'get_cell_heat_copy', size(copy))
  end subroutine

  subroutine unsupported_get_face_temp_copy(this, copy)
    class(thermal_species_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    call unsupported_solver_operation(this, 'get_face_temp_copy', size(copy))
  end subroutine

  subroutine unsupported_get_face_temp_view(this, view)
    class(thermal_species_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    call unsupported_solver_operation(this, 'get_face_temp_view')
    view => null()
  end subroutine

  subroutine unsupported_get_cell_temp_grad(this, tgrad)
    class(thermal_species_solver), intent(inout) :: this
    real(r8), intent(out) :: tgrad(:,:)
    call unsupported_solver_operation(this, 'get_cell_temp_grad', size(tgrad))
  end subroutine

  subroutine unsupported_get_cell_conc_copy(this, n, copy)
    class(thermal_species_solver), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(inout) :: copy(:)
    call unsupported_solver_operation(this, 'get_cell_conc_copy', n + size(copy))
  end subroutine

  subroutine unsupported_set_ext_enthalpy_rate(this, enthalpy_rate)
    class(thermal_species_solver), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)
    call unsupported_solver_operation(this, 'set_ext_enthalpy_rate', size(enthalpy_rate))
  end subroutine

  subroutine unsupported_set_ext_species_rate(this, n, species_rate)
    class(thermal_species_solver), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: species_rate(:)
    call unsupported_solver_operation(this, 'set_ext_species_rate', n + size(species_rate))
  end subroutine

  subroutine unsupported_solver_operation(this, name, n)
    class(thermal_species_solver), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(in), optional :: n
    logical :: impossible
    impossible = .not.same_type_as(this, this) .or. len(name) == 0
    if (present(n)) impossible = impossible .or. n == huge(n)
    if (impossible) return
    INSIST(.false.)
  end subroutine

end module thermal_species_solver_class
