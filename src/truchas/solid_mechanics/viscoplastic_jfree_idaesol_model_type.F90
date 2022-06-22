!!
!! Jacobian-free viscoplastic idaesol model
!!
!! This version closely mirrors what the legacy solver did, but on top
!! of the object-oriented idaesol_type.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! May 2022
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module viscoplastic_jfree_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use idaesol_type
  use viscoplastic_model_type
  implicit none
  private

  type, extends(idaesol_model), public :: viscoplastic_jfree_idaesol_model
    private
    type(viscoplastic_model), pointer :: vp_model => null() ! unowned reference
    real(r8) :: atol, rtol
    real(r8) :: u0(6), h ! time step size
    logical :: updated_precon = .true.
  contains
    procedure :: init
    procedure :: size => model_size
    procedure :: du_norm
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
  end type viscoplastic_jfree_idaesol_model

contains

  subroutine init(this, vp_model, atol, rtol)
    class(viscoplastic_jfree_idaesol_model), intent(out) :: this
    type(viscoplastic_model), intent(in), target :: vp_model
    real(r8), intent(in) :: atol, rtol
    this%vp_model => vp_model
    this%atol = atol
    this%rtol = rtol
  end subroutine init

  subroutine du_norm(this, t, u, du, error)
    class(viscoplastic_jfree_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:), du(:)
    real(r8), intent(out) :: error
    error = maxval(abs(du)/(this%atol + this%rtol*abs(u)))
  end subroutine du_norm

  integer function model_size(this)
    class(viscoplastic_jfree_idaesol_model), intent(in) :: this
    model_size = 6
  end function model_size

  subroutine compute_f(this, t, u, udot, f)
    class(viscoplastic_jfree_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:), udot(:)
    real(r8), intent(out), contiguous :: f(:)
    call this%vp_model%compute_udot(t, u, f) ! dstrain_plastic / dt
    if (this%updated_precon) then
      this%u0 = u - this%h * udot
      this%updated_precon = .false.
    end if
  end subroutine compute_f

  subroutine compute_precon(this, t, u, dt)
    class(viscoplastic_jfree_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), contiguous :: u(:)
    this%updated_precon = .true.
    this%h = dt
  end subroutine compute_precon

  subroutine apply_precon(this, t, u, f)
    external dgesv ! LAPACK
    class(viscoplastic_jfree_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:)
    real(r8), intent(inout), contiguous :: f(:)
    f = u - this%u0 - this%h*f
  end subroutine apply_precon

end module viscoplastic_jfree_idaesol_model_type
