!!
!! HT_2D_NORM_TYPE
!!
!! This module defines the correction norm used by the implicit integrator for
!! the 2D thermal transport time-step system. It combines temperature and
!! enthalpy error measures using the tolerances configured for the solver.
!!
!! David Neill-Asanza <davidhneill@gmail.com>, May 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ht_2d_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_2d_model_type
  use ht_2d_vector_type
  use parallel_communication, only: global_maxval
  implicit none
  private

  type, public :: ht_2d_norm
    private
    ! type(unstr_2d_mesh), pointer :: mesh => null()  ! reference only -- do not own
    type(ht_2d_model), pointer :: model => null()   ! reference only -- do not own
    real(r8) :: abs_T_tol   ! absolute temperature tolerance
    real(r8) :: rel_T_tol   ! relative temperature tolerance
    real(r8) :: abs_H_tol   ! absolute enthalpy tolerance
    real(r8) :: rel_H_tol   ! relative enthalpy tolerance
  contains
    procedure :: init
    procedure :: compute
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(ht_2d_norm), intent(out) :: this
    type(ht_2d_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context

    this%model => model

    context = 'processing ' // params%path() // ': '
    call params%get('abs-t-tol', this%abs_T_tol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%abs_T_tol < 0.0_r8) then
      stat = 1
      errmsg = context // '"abs-t-tol" must be >= 0.0'
      return
    end if

    call params%get('rel-t-tol', this%rel_T_tol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%rel_T_tol < 0.0_r8) then
      stat = 1
      errmsg = context // '"rel-t-tol" must be >= 0.0'
      return
    end if
    if (.not. valid_tol(this%abs_T_tol, this%rel_T_tol)) then
      stat = 1
      errmsg = context // '"abs-t-tol" and "rel-t-tol" cannot both be 0.0'
      return
    end if

    call params%get('abs-h-tol', this%abs_H_tol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%abs_H_tol < 0.0_r8) then
      stat = 1
      errmsg = context // '"abs-h-tol" must be >= 0.0'
      return
    end if

    call params%get('rel-h-tol', this%rel_H_tol, stat, errmsg, default=this%rel_T_tol)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%rel_H_tol < 0.0_r8) then
      stat = 1
      errmsg = context // '"rel-h-tol" must be >= 0.0'
      return
    end if
    if (.not. valid_tol(this%abs_H_tol, this%rel_H_tol)) then
      stat = 1
      errmsg = context // '"abs-h-tol" and "rel-h-tol" cannot both be 0.0'
      return
    end if
    stat = 0

  contains

    elemental logical function valid_tol(abs_tol, rel_tol)
      real(r8), intent(in) :: abs_tol, rel_tol
      valid_tol = (abs_tol >= 0.0_r8) .and. (rel_tol >= 0.0_r8) .and. &
                  (abs_tol > 0.0_r8 .or. rel_tol > 0.0_r8)
    end function

  end subroutine init


  subroutine compute(this, t, u, du, du_norm)

    class(ht_2d_norm), intent(in) :: this
    real(r8), intent(in) :: t
    type(ht_2d_vector), intent(in) :: u, du
    real(r8), intent(out) :: du_norm

    associate (ncell_onP => this%model%mesh%ncell_onP, nface_onP => this%model%mesh%nface_onP)
      du_norm = maxerr(u%hc(:ncell_onP), du%hc(:ncell_onP), this%abs_H_tol, this%rel_H_tol)
      du_norm = max(du_norm, &
          maxerr(u%tc(:ncell_onP), du%tc(:ncell_onP), this%abs_T_tol, this%rel_T_tol))
      du_norm = max(du_norm, &
          maxerr(u%tf(:nface_onP), du%tf(:nface_onP), this%abs_T_tol, this%rel_T_tol))
    end associate
    du_norm = global_maxval(du_norm)

  contains

    real(r8) function maxerr(u, du, atol, rtol)
      real(r8), intent(in) :: u(:), du(:), atol, rtol
      maxerr = maxval(abs(du)/(atol + rtol*abs(u)))
    end function

  end subroutine compute

end module ht_2d_norm_type
