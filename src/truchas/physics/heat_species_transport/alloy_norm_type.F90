!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use alloy_vector_type
  use alloy_model_type
  use parallel_communication
  implicit none
  private

  type, public :: alloy_norm
    private
    type(alloy_model), pointer :: model => null()
    real(r8) :: temp_atol
    real(r8) :: temp_rtol
    real(r8) :: enth_atol
    real(r8) :: enth_rtol
  contains
    procedure :: init
    procedure :: compute
  end type alloy_norm

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(alloy_norm), intent(out) :: this
    type(alloy_model), intent(in), target :: model
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context

    this%model => model

    context = 'processing ' // params%path() // ': '
    call params%get('abs-t-tol', this%temp_atol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%temp_atol < 0) then
      stat = 1
      errmsg = context // '"abs-t-tol" must be >= 0.0'
      return
    end if

    call params%get('rel-t-tol', this%temp_rtol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%temp_rtol < 0) then
      stat = 1
      errmsg = context // '"rel-t-tol" must be >= 0.0'
      return
    end if

    if (this%temp_atol == 0 .and. this%temp_rtol == 0) then
      stat = 1
      errmsg = context // '"abs-t-tol" and "rel-t-tol" cannot both be 0.0'
      return
    end if

    call params%get('abs-h-tol', this%enth_atol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%enth_atol < 0) then
      stat = 1
      errmsg = context // '"abs-h-tol" must be >= 0.0'
      return
    end if

    call params%get('rel-h-tol', this%enth_rtol, stat, errmsg, default=this%temp_rtol)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%enth_rtol < 0) then
      stat = 1
      errmsg = context // '"rel-h-tol" must be >= 0.0'
      return
    end if

    if (this%enth_atol == 0 .and. this%enth_rtol == 0) then
      stat = 1
      errmsg = context // '"abs-h-tol" and "rel-h-tol" cannot both be 0.0'
      return
    end if

  end subroutine init


  subroutine compute(this, t, u, du, du_norm)

    class(alloy_norm), intent(in) :: this
    real(r8), intent(in) :: t
    type(alloy_vector), intent(in) :: u, du
    real(r8), intent(out) :: du_norm

    integer :: n
    real(r8) :: qerror
    real(r8), allocatable :: res(:)
    integer :: ncell, nface

    ncell = this%model%mesh%ncell_onP
    nface = this%model%mesh%nface_onP

    du_norm = 0.0_r8
    du_norm = max(du_norm, maxerr(u%hc(:ncell), du%hc(:ncell), this%enth_atol, this%enth_rtol))
    du_norm = max(du_norm, maxerr(u%tc(:ncell), du%tc(:ncell), this%temp_atol, this%temp_rtol))
    du_norm = max(du_norm, maxerr(u%tf(:nface), du%tf(:nface), this%temp_atol, this%temp_rtol))
    du_norm = global_maxval(du_norm)

  contains

    real(r8) function maxerr (u, du, atol, rtol)
      real(r8), intent(in) :: u(:), du(:), atol, rtol
      real(r8) :: array(size(du))
      array = abs(du) / (atol + rtol*abs(u))
      maxerr = maxval(array)
    end function maxerr

  end subroutine compute

end module alloy_norm_type
