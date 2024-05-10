!TODO: finish documentation
!! HT_2D_NORM_TYPE
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

module HT_2d_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use HT_2d_model_type
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  type, public:: HT_2d_norm
    private
    ! type(unstr_2d_mesh), pointer :: mesh => null()  ! reference only -- do not own
    type(HT_2d_model), pointer :: model => null()   ! reference only -- do not own
    real(r8) :: abs_T_tol   ! absolute temperature tolerance
    real(r8) :: rel_T_tol   ! relative temperature tolerance
    real(r8) :: abs_H_tol   ! absolute enthalpy tolerance
    real(r8) :: rel_H_tol   ! relative enthalpy tolerance
  contains
    procedure :: init
    procedure :: compute
  end type HT_2d_norm

contains

  subroutine init(this, model, params)

    class(HT_2d_norm), intent(out) :: this
    type(HT_2d_model), intent(in), target :: model
    type(parameter_list) :: params

    integer :: stat
    character(:), allocatable :: context, errmsg

    this%model => model

    context = 'processing ' // params%path() // ': '
    call params%get('temp-abs-tol', this%abs_T_tol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) call TLS_fatal(context//errmsg)
    call params%get('temp-rel-tol', this%rel_T_tol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) call TLS_fatal(context//errmsg)
    if (this%abs_T_tol < 0.0_r8) call TLS_fatal(context//'"temp-abs-tol" must be >= 0.0')
    if (this%rel_T_tol < 0.0_r8) call TLS_fatal(context//'"temp-rel-tol" must be >= 0.0')
    if (this%abs_T_tol == 0.0_r8 .and. this%rel_T_tol == 0.0_r8) &
        call TLS_fatal(context//'"temp-abs-tol" and "temp-rel-tol" cannot both be 0.0')

    call params%get('enth-abs-tol', this%abs_H_tol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) call TLS_fatal(context//errmsg)
    call params%get('enth-rel-tol', this%rel_H_tol, stat, errmsg, default=0.0_r8)
    if (stat /= 0) call TLS_fatal(context//errmsg)
    if (this%abs_H_tol < 0.0_r8) call TLS_fatal(context//'"enth-abs-tol" must be >= 0.0')
    if (this%rel_H_tol < 0.0_r8) call TLS_fatal(context//'"enth-rel-tol" must be >= 0.0')
    if (this%abs_H_tol == 0.0_r8 .and. this%rel_H_tol == 0.0_r8) &
        call TLS_fatal(context//'"enth-abs-tol" and "enth-rel-tol" cannot both be 0.0')

  end subroutine init


  subroutine compute(this, u, du, du_norm)

    use parallel_communication, only: global_maxval

    class(HT_2d_norm), intent(in) :: this
    real(r8), intent(in), target :: u(:), du(:)
    real(r8), intent(out) :: du_norm

    real(r8), pointer :: useg(:), duseg(:)

    ASSERT(size(u) == size(du))
    ASSERT(size(u) == this%model%num_dof())

    du_norm = 0.0_r8

    !! Cell temperature delta norm
    call this%model%get_cell_temp_view(u, useg)
    call this%model%get_cell_temp_view(du, duseg)
    du_norm = max(du_norm, maxerr(useg, duseg, this%abs_T_tol, this%rel_T_tol))

    !! Face temperature delta norm
    call this%model%get_face_temp_view(u, useg)
    call this%model%get_face_temp_view(du, duseg)
    du_norm = max(du_norm, maxerr(useg, duseg, this%abs_T_tol, this%rel_T_tol))

    !! Cell enthalpy delta norm
    call this%model%get_cell_heat_view(u, useg)
    call this%model%get_cell_heat_view(du, duseg)
    du_norm = max(du_norm, maxerr(useg, duseg, this%abs_H_tol, this%rel_H_tol))

    du_norm = global_maxval(du_norm)

  contains

    real(r8) function maxerr(u, du, atol, rtol)
      real(r8), intent(in) :: u(:), du(:), atol, rtol
      maxerr = maxval(abs(du)/(atol + rtol*abs(u)))
    end function

  end subroutine compute

end module HT_2d_norm_type
