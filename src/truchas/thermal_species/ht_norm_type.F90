!!
!! HT_NORM_TYPE
!!
!! This module defines the norm used by the implicit integrator for the thermal
!! transport time-step system. It combines thermal and optional view-factor
!! radiation error measures using the tolerances configured for the solver.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module ht_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_vector_type
  use ht_model_type
  use parallel_communication
  implicit none
  private

  type, public :: ht_norm
    private
    type(ht_model), pointer :: model => null() ! unowned reference
    real(r8) :: abs_T_tol
    real(r8) :: rel_T_tol
    real(r8) :: abs_H_tol
    real(r8) :: rel_H_tol
    real(r8) :: rad_tol
  contains
    procedure :: init
    procedure :: compute
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(ht_norm), intent(out) :: this
    type(ht_model), intent(in), target :: model
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
    if (.not.valid_tol(this%abs_T_tol, this%rel_T_tol)) then
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
    if (.not.valid_tol(this%abs_H_tol, this%rel_H_tol)) then
      stat = 1
      errmsg = context // '"abs-h-tol" and "rel-h-tol" cannot both be 0.0'
      return
    end if

    call params%get('rad-tol', this%rad_tol, stat, errmsg, default=1.0e-3_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%rad_tol <= 0.0_r8) then
      stat = 1
      errmsg = context // '"rad-tol" must be > 0.0'
      return
    end if
    stat = 0

  contains

    elemental logical function valid_tol(abs_tol, rel_tol)
      real(r8), intent(in) :: abs_tol, rel_tol
      valid_tol = (abs_tol >= 0.0_r8) .and. (rel_tol >= 0.0_r8) .and. &
                  ((abs_tol > 0.0_r8 .or. rel_tol > 0.0_r8))
    end function

  end subroutine init


  subroutine compute(this, t, u, du, du_norm)

    class(ht_norm), intent(in) :: this
    real(r8), intent(in) :: t
    type(ht_vector), intent(in) :: u, du
    real(r8), intent(out) :: du_norm

    real(r8) :: qerror

    du_norm = 0.0_r8
    associate (ncell_onP => this%model%mesh%ncell_onP, nface_onP => this%model%mesh%nface_onP)
      du_norm = max(du_norm, &
          maxerr(u%hc(:ncell_onP), du%hc(:ncell_onP), this%abs_H_tol, this%rel_H_tol, &
              this%model%void_cell))
      du_norm = max(du_norm, &
          maxerr(u%tc(:ncell_onP), du%tc(:ncell_onP), this%abs_T_tol, this%rel_T_tol, &
              this%model%void_cell))
      du_norm = max(du_norm, &
          maxerr(u%tf(:nface_onP), du%tf(:nface_onP), this%abs_T_tol, this%rel_T_tol, &
              this%model%void_face))
    end associate
    du_norm = global_maxval(du_norm)
    !! View-factor radiation algebraic residual norm. The BDF2 integrator
    !! asks for a correction norm, but the algebraic radiation block is tested
    !! with its global relative residual norm.
    if (this%model%vf_rad%is_active()) then
      call this%model%vf_rad%relative_residual_norm(t, u%tf, u%qrad, qerror)
      du_norm = max(du_norm, qerror/this%rad_tol)
    end if

  contains

    real(r8) function maxerr(u, du, atol, rtol, void)
      real(r8), intent(in) :: u(:), du(:), atol, rtol
      logical, intent(in), optional :: void(:)
      real(r8) :: array(size(du))
      if (present(void)) then
        where (void(:size(du)))
          array = 0.0_r8
        elsewhere
          array = abs(du) / (atol + rtol*abs(u))
        end where
      else
        array = abs(du) / (atol + rtol*abs(u))
      end if
      maxerr = maxval(array)
    end function

  end subroutine compute

end module ht_norm_type
