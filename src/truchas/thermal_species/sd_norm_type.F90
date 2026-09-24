!!
!! SD_NORM_TYPE
!!
!! This module defines the norm used by the implicit integrator for the species
!! transport time-step system. It combines the species component error measures
!! using the tolerances configured for the solver.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module sd_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sd_model_type
  use sd_vector_type
  use parallel_communication
  implicit none
  private

  type, public :: sd_norm
    private
    type(sd_model), pointer :: model ! unowned reference
    real(r8), allocatable :: abs_C_tol(:)
    real(r8), allocatable :: rel_C_tol(:)
  contains
    procedure :: init
    procedure :: compute
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type
    class(sd_norm), intent(out) :: this
    type(sd_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    character(:), allocatable :: context

    this%model => model

    context = 'processing ' // params%path() // ': '
    call params%get('abs-c-tol', this%abs_c_tol, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (size(this%abs_C_tol) /= model%num_comp) then
      stat = 1
      errmsg = context // '"abs-c-tol" must have one value for each species component'
      return
    end if

    call params%get('rel-c-tol', this%rel_c_tol, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (size(this%rel_C_tol) /= model%num_comp) then
      stat = 1
      errmsg = context // '"rel-c-tol" must have one value for each species component'
      return
    end if
    if (any(this%abs_C_tol < 0.0_r8)) then
      stat = 1
      errmsg = context // '"abs-c-tol" values must be >= 0.0'
      return
    else if (any(this%rel_C_tol < 0.0_r8)) then
      stat = 1
      errmsg = context // '"rel-c-tol" values must be >= 0.0'
      return
    else if (.not.all(valid_tol(this%abs_C_tol, this%rel_C_tol))) then
      stat = 1
      errmsg = context // '"abs-c-tol" and "rel-c-tol" cannot both be 0.0 for any species component'
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

    class(sd_norm), intent(in) :: this
    real(r8), intent(in) :: t
    type(sd_vector), intent(in) :: u, du
    real(r8), intent(out) :: du_norm

    integer :: n
    real(r8) :: sd_du_norm

    sd_du_norm = 0.0_r8
    associate (ncell_onP => this%model%mesh%ncell_onP, nface_onP => this%model%mesh%nface_onP)
      do n = 1, this%model%num_comp
        sd_du_norm = max(sd_du_norm, maxerr(u%cc(:ncell_onP,n), du%cc(:ncell_onP,n), &
                                            this%abs_C_tol(n), this%rel_C_tol(n), &
                                            this%model%void_cell))
        sd_du_norm = max(sd_du_norm, maxerr(u%cf(:nface_onP,n), du%cf(:nface_onP,n), &
                                            this%abs_C_tol(n), this%rel_C_tol(n), &
                                            this%model%void_face))
        sd_du_norm = global_maxval(sd_du_norm)
      end do
    end associate
    du_norm = sd_du_norm

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

end module sd_norm_type
