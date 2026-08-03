!!
!! HT_IDAESOL_MODEL_TYPE
!!
!! This module defines the IDAESOL adapter for the thermal transport model. It
!! bundles the thermal model, preconditioner, and norm objects behind the
!! abstract interface used by the implicit DAE integrator.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! IDAESOL operates on abstract vectors and only requires on-process entries
!! to be valid at callback boundaries. Thermal model and preconditioner
!! routines that need off-process input entries, including temperature and
!! enthalpy state data, are responsible for gathering them, and only on-process
!! output entries are significant on return to IDAESOL.
!!

#include "f90_assert.fpp"

module ht_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use new_idaesol_type, only: idaesol_model
  use vector_class
  use ht_vector_type
  use ht_model_type
  use ht_precon_type
  use ht_norm_type
  use truchas_timers, only: start_timer, stop_timer
  implicit none
  private

  type, extends(idaesol_model), public :: ht_idaesol_model
    type(ht_model),  pointer :: model =>  null() ! unowned reference
    type(ht_precon), pointer :: precon => null() ! unowned reference
    type(ht_norm),   pointer :: norm   => null() ! unowned reference
  contains
    procedure :: init
    !! Deferred procedures from IDAESOL_MODEL
    procedure :: alloc_vector
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
  end type

contains

  subroutine init(this, model, precon, norm)
    class(ht_idaesol_model), intent(out) :: this
    type(ht_model),  intent(in), target :: model
    type(ht_precon), intent(in), target :: precon
    type(ht_norm),   intent(in), target :: norm
    this%model => model
    this%precon => precon
    this%norm => norm
    ASSERT(associated(this%model, precon%model))
  end subroutine

  subroutine alloc_vector(this, vec)
    class(ht_idaesol_model), intent(in) :: this
    class(vector), allocatable, intent(out) :: vec
    type(ht_vector), allocatable :: tmp
    allocate(tmp)
    call this%model%init_vector(tmp)
    call move_alloc(tmp, vec)
  end subroutine

  subroutine compute_f(this, t, u, udot, f)
    class(ht_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u, udot ! data is intent(in)
    class(vector), intent(inout) :: f       ! data is intent(out)
    call start_timer('residual')
    select type (u)
    class is (ht_vector)
      select type (udot)
      class is (ht_vector)
        select type (f)
        class is (ht_vector)
          call this%model%residual(t, u, udot, f)
        end select
      end select
    end select
    call stop_timer('residual')
  end subroutine

  subroutine apply_precon(this, t, u, f)
    class(ht_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u  ! data is intent(in)
    class(vector), intent(inout) :: f  ! data is intent(inout)
    call start_timer('precon apply')
    select type (u)
    class is (ht_vector)
      select type (f)
      class is (ht_vector)
        call this%precon%apply(t, u, f)
      end select
    end select
    call stop_timer('precon apply')
  end subroutine

  subroutine compute_precon(this, t, u, udot, dt)
    class(ht_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    class(vector), intent(inout) :: u, udot
    call start_timer('precon compute')
    select type (u)
    class is (ht_vector)
      call this%precon%compute(t, u, dt)
    end select
    call stop_timer('precon compute')
  end subroutine

  subroutine du_norm(this, t, u, du, error)
    class(ht_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u, du
    real(r8), intent(out) :: error
    select type (u)
    class is (ht_vector)
      select type (du)
      class is (ht_vector)
        call this%norm%compute(t, u, du, error)
      end select
    end select
  end subroutine

end module ht_idaesol_model_type
