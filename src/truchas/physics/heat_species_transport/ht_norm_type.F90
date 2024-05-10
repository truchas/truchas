!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ht_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_vector_type
  use ht_model_type
  use rad_problem_type
  use parallel_communication
  implicit none
  private

  type, public :: ht_norm
    private
    type(ht_model), pointer :: model => null()
    real(r8) :: temp_atol
    real(r8) :: temp_rtol
    real(r8) :: enth_atol
    real(r8) :: enth_rtol
  contains
    procedure :: init
    procedure :: compute
  end type ht_norm

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(ht_norm), intent(out) :: this
    type(ht_model), intent(in), target :: model
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

    class(ht_norm), intent(in) :: this
    real(r8), intent(in) :: t
    type(ht_vector), intent(in) :: u, du
    real(r8), intent(out) :: du_norm

    integer :: n
    real(r8) :: qerror
    real(r8), allocatable :: res(:)
    integer :: ncell, nface

    ncell = this%model%mesh%ncell_onP
    nface = this%model%mesh%nface_onP

    du_norm = 0.0_r8
    du_norm = max(du_norm, maxerr(u%hc(:ncell), du%hc(:ncell), this%enth_atol, this%enth_rtol, this%model%void_cell))
    du_norm = max(du_norm, maxerr(u%tc(:ncell), du%tc(:ncell), this%temp_atol, this%temp_rtol, this%model%void_cell))
    du_norm = max(du_norm, maxerr(u%tf(:nface), du%tf(:nface), this%temp_atol, this%temp_rtol, this%model%void_face))
    du_norm = global_maxval(du_norm)
    !! Enclosure radiation system error norms
    !TODO! This is a quick hack that needs to be fixed.  The tolerance
    !TODO! is hardwired, and  we are (re)computing the residual.  The BDF2
    !TODO! integrator thinks it's computing the norm of the correction du but
    !TODO! in this case it is the actual residual norm.  We need more general
    !TODO! norms in the integrator.
    if (associated(this%model%vf_rad_prob)) then
      !if (is_IOP) write(*,'(a)',advance='no')'ER error: ||res||/||rhs||='
      do n = 1, size(this%model%vf_rad_prob)
        associate (q => u%encl(n)%qrad, faces => this%model%vf_rad_prob(n)%faces)
          allocate(res(size(faces)))
          call this%model%vf_rad_prob(n)%residual(t, q, u%tf(faces), res)
          qerror = sqrt(global_sum(norm2(res)**2))
          call this%model%vf_rad_prob(n)%rhs(t, u%tf(faces), res)
          qerror = qerror / sqrt(global_sum(norm2(res)**2))
          !if (is_IOP) write(*,'(e10.3)',advance='no') qerror
          du_norm = max(du_norm, qerror/1.0d-3)
          deallocate(res)
        end associate
      end do
      !if (is_IOP) write(*,*)
    end if

  contains

    real(r8) function maxerr (u, du, atol, rtol, void)
      real(r8), intent(in) :: u(:), du(:), atol, rtol
      logical, pointer :: void(:)
      real(r8) :: array(size(du))
      if (associated(void)) then
        where (void(:size(du)))
          array = 0.0_r8
        elsewhere
          array = abs(du) / (atol + rtol*abs(u))
        end where
      else
        array = abs(du) / (atol + rtol*abs(u))
      end if
      maxerr = maxval(array)
    end function maxerr

  end subroutine compute

end module ht_norm_type
