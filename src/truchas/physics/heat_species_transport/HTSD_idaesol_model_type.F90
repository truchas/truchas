!!
!! HTSD_IDAESOL_MODEL_TYPE
!!
!! This module defines an extension of the IDAESOL_MODEL abstract class that
!! implements the methods required by the ODE integrator. It bundles several
!! different computational pieces for the heat transfer/species diffusion model
!! and adapts them to the required interface.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Updated to object oriented F2003, May 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module htsd_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use idaesol_type, only: idaesol_model
  use htsd_model_type
  use htsd_precon_type
  use htsd_norm_type
  implicit none
  private

  type, extends(idaesol_model), public :: htsd_idaesol_model
    type(htsd_model),  pointer :: model =>  null() ! reference only -- not owned
    type(htsd_precon), pointer :: precon => null() ! reference only -- not owned
    type(htsd_norm),   pointer :: norm   => null() ! reference only -- not owned
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
    class(htsd_idaesol_model), intent(out) :: this
    type(htsd_model),  intent(in), target :: model
    type(htsd_precon), intent(in), target :: precon
    type(htsd_norm),   intent(in), target :: norm
    this%model => model
    this%precon => precon
    this%norm => norm
    ASSERT(associated(this%model, precon%model))
  end subroutine init

  integer function model_size(this)
    class(htsd_idaesol_model), intent(in) :: this
    model_size = htsd_model_size(this%model)
  end function model_size

  subroutine compute_f(this, t, u, udot, f)
    class(htsd_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), contiguous, intent(in)  :: u(:), udot(:)
    real(r8), contiguous, intent(out) :: f(:)
    call htsd_model_compute_f(this%model, t, u, udot, f)
  end subroutine compute_f

  subroutine apply_precon(this, t, u, f)
    class(htsd_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), contiguous, intent(in) :: u(:)
    real(r8), contiguous, intent(inout) :: f(:)
    call htsd_precon_apply(this%precon, t, u, f)
  end subroutine apply_precon

  subroutine compute_precon(this, t, u, dt)
    class(htsd_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    real(r8), contiguous, intent(in) :: u(:)
    integer :: errc
    call htsd_precon_compute(this%precon, t, u, dt, errc)
    INSIST(errc == 0)
  end subroutine compute_precon

  subroutine du_norm(this, u, du, error)
    class(htsd_idaesol_model) :: this
    real(r8), contiguous, intent(in) :: u(:), du(:)
    real(r8), intent(out) :: error
    call htsd_norm_compute(this%norm, u, du, error)
  end subroutine du_norm

end module htsd_idaesol_model_type
