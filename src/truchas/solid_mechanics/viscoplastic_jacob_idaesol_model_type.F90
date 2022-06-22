!!
!! Zach Jibben <zjibben@lanl.gov>
!! October 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module viscoplastic_jacob_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use idaesol_type
  use viscoplastic_model_type
  implicit none
  private

  type, extends(idaesol_model), public :: viscoplastic_jacob_idaesol_model
    private
    type(viscoplastic_model), pointer :: vp_model => null() ! unowned reference
    real(r8) :: atol, rtol
    real(r8) :: precon(6,6)
    integer :: ipiv(6)
  contains
    procedure :: init
    procedure :: size => model_size
    procedure :: du_norm
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
  end type viscoplastic_jacob_idaesol_model

contains

  subroutine init(this, vp_model, atol, rtol)
    class(viscoplastic_jacob_idaesol_model), intent(out) :: this
    type(viscoplastic_model), intent(in), target :: vp_model
    real(r8), intent(in) :: atol, rtol
    this%vp_model => vp_model
    this%atol = atol
    this%rtol = rtol
  end subroutine init

  subroutine du_norm(this, t, u, du, error)
    class(viscoplastic_jacob_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:), du(:)
    real(r8), intent(out) :: error
    error = maxval(abs(du)/(this%atol + this%rtol*abs(u)))
  end subroutine du_norm

  integer function model_size(this)
    class(viscoplastic_jacob_idaesol_model), intent(in) :: this
    model_size = 6
  end function model_size


  subroutine compute_f(this, t, u, udot, f)
    class(viscoplastic_jacob_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:), udot(:)
    real(r8), intent(out), contiguous :: f(:)
    call this%vp_model%compute_udot(t, u, f) ! dstrain_plastic / dt
    f = f - udot ! residual
  end subroutine compute_f


  subroutine compute_precon(this, t, u, dt)
    external dgetrf ! LAPACK
    class(viscoplastic_jacob_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), contiguous :: u(:)
    integer :: j, stat
    call this%vp_model%compute_precon(t, u, this%precon)
    do j = 1, 6
      this%precon(j,j) = this%precon(j,j) - 1 / dt
    end do
    call dgetrf(6, 6, this%precon, 6, this%ipiv, stat)
    INSIST(stat == 0)
  end subroutine compute_precon


  subroutine apply_precon(this, t, u, f)
    external dgetrs ! LAPACK
    class(viscoplastic_jacob_idaesol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:)
    real(r8), intent(inout), contiguous :: f(:)
    integer :: stat
    call dgetrs('N', 6, 1, this%precon, 6, this%ipiv, f, 6, stat)
    INSIST(stat == 0)
  end subroutine apply_precon

end module viscoplastic_jacob_idaesol_model_type
