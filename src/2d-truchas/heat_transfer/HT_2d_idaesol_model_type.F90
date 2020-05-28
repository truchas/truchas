!TODO: finish documentation
!! HT_2D_IDAESOL_MODEL_TYPE
!!
!! This module defines an extension of the IDAESOL_MODEL abstract class that
!! implements the methods required by the ODE integrator. It bundles several
!! different computational pieces for the 2D heat transfer model and adapts
!! them to the required interface.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! May 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HT_2d_idaesol_model_type

  use kinds, only: r8
  use idaesol_type, only: idaesol_model
  use HT_2d_model_type
  use HT_2d_precon_type
  use HT_2d_norm_type
  implicit none
  private

  type, extends(idaesol_model), public :: HT_2d_idaesol_model
    type(HT_2d_model),  pointer :: model =>  null() ! reference only -- not owned
    type(HT_2d_precon), pointer :: precon => null() ! reference only -- not owned
    type(HT_2d_norm),   pointer :: norm   => null() ! reference only -- not owned
  contains
    procedure :: init
    !! Deferred procedures from IDAESOL_MODEL
    procedure :: size => model_size
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
  end type

contains

  subroutine init(this, model, precon, norm)
    class(HT_2d_idaesol_model), intent(out) :: this
    type(HT_2d_model),  intent(in), target :: model
    type(HT_2d_precon), intent(in), target :: precon
    type(HT_2d_norm),   intent(in), target :: norm
    this%model => model
    this%precon => precon
    this%norm => norm
    ASSERT(associated(this%model, precon%model))
  end subroutine init

  integer function model_size(this)
    class(HT_2d_idaesol_model), intent(in) :: this
    model_size = this%model%num_dof()
  end function model_size

  subroutine compute_f(this, t, u, udot, f)
    class(HT_2d_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), contiguous, intent(in)  :: u(:), udot(:)
    real(r8), contiguous, intent(out) :: f(:)
    call this%model%compute_f(t, u, udot, f)
  end subroutine compute_f

  !TODO: t and u are not used!
  subroutine apply_precon(this, t, u, f)
    class(HT_2d_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), contiguous, intent(in) :: u(:)
    real(r8), contiguous, intent(inout) :: f(:)
    call this%precon%apply(f)
  end subroutine apply_precon

  subroutine compute_precon(this, t, u, dt)
    class(HT_2d_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    real(r8), contiguous, intent(in) :: u(:)
    call this%precon%compute(t, u, dt)
  end subroutine compute_precon

  subroutine du_norm(this, t, u, du, error)
    class(HT_2d_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), contiguous, intent(in) :: u(:), du(:)
    real(r8), intent(out) :: error
    call this%norm%compute(u, du, error)
  end subroutine du_norm

end module HT_2d_idaesol_model_type
