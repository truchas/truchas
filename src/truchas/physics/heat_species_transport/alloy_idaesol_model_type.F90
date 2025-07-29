!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use new_idaesol_type, only: idaesol_model
  use vector_class
  use alloy_vector_type
  use alloy_model_type
  use alloy_precon_type
  use alloy_norm_type
  implicit none
  private

  type, extends(idaesol_model), public :: alloy_idaesol_model
    type(alloy_model),  pointer :: model =>  null() ! reference only -- not owned
    type(alloy_precon), pointer :: precon => null() ! reference only -- not owned
    type(alloy_norm),   pointer :: norm   => null() ! reference only -- not owned
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
    class(alloy_idaesol_model), intent(out) :: this
    type(alloy_model),  intent(in), target :: model
    type(alloy_precon), intent(in), target :: precon
    type(alloy_norm),   intent(in), target :: norm
    this%model => model
    this%precon => precon
    this%norm => norm
    ASSERT(associated(this%model, precon%model))
  end subroutine

  subroutine alloc_vector(this, vec)
    class(alloy_idaesol_model), intent(in) :: this
    class(vector), allocatable, intent(out) :: vec
    type(alloy_vector), allocatable :: tmp
    allocate(tmp)
    call this%model%init_vector(tmp)
    call move_alloc(tmp, vec)
  end subroutine

  subroutine compute_f(this, t, u, udot, f)
    class(alloy_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u, udot ! data is intent(in)
    class(vector), intent(inout) :: f       ! data is intent(out)
    select type (u)
    class is (alloy_vector)
      select type (udot)
      class is (alloy_vector)
        select type (f)
        class is (alloy_vector)
          call this%model%compute_f(t, u, udot, f)
        end select
      end select
    end select
  end subroutine

  subroutine apply_precon(this, t, u, f)
    class(alloy_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u  ! data is intent(in)
    class(vector), intent(inout) :: f  ! data is intent(inout)
    select type (u)
    class is (alloy_vector)
      select type (f)
      class is (alloy_vector)
        call this%precon%apply(t, u, f)
      end select
    end select
  end subroutine

  subroutine compute_precon(this, t, u, udot, dt)
    class(alloy_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    class(vector), intent(inout) :: u, udot
    select type (u)
    class is (alloy_vector)
      select type (udot)
      class is (alloy_vector)
        call this%precon%compute(t, u, udot, dt)
      end select
    end select
  end subroutine

  subroutine du_norm(this, t, u, du, error)
    class(alloy_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u, du
    real(r8), intent(out) :: error
    select type (u)
    class is (alloy_vector)
      select type (du)
      class is (alloy_vector)
        call this%norm%compute(t, u, du, error)
      end select
    end select
  end subroutine

end module alloy_idaesol_model_type
