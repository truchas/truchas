!!
!! FHT_NORM_TYPE
!!
!! This module defines the residual norm used by the nonadaptive integrator
!! for the flow-coupled thermal transport solver. It combines the thermal
!! residual norm on non-void cells/faces with the optional view-factor
!! radiation residual norm.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module fht_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fht_model_type
  use fht_vector_type
  use parallel_communication
  implicit none
  private

  type, public :: fht_norm
    private
    type(fht_model), pointer :: model => null() ! unowned reference
    real(r8) :: abs_tol
    real(r8) :: rel_tol
    real(r8) :: rad_tol
    real(r8) :: err0
    logical  :: verbose
    integer  :: unit
  contains
    procedure :: init
    procedure :: fnorm
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(fht_norm), intent(out) :: this
    type(fht_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context

    this%model => model

    context = 'processing ' // params%path() // ': '
    call params%get('abs-tol', this%abs_tol, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%abs_tol < 0.0_r8) then
      stat = 1
      errmsg = context // '"abs-tol" must be >= 0.0'
      return
    end if

    call params%get('rel-tol', this%rel_tol, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    else if (this%rel_tol < 0.0_r8) then
      stat = 1
      errmsg = context // '"rel-tol" must be >= 0.0'
      return
    end if
    if (.not.valid_tol(this%abs_tol, this%rel_tol)) then
      stat = 1
      errmsg = context // '"abs-tol" and "rel-tol" cannot both be 0.0'
      return
    end if

    call params%get('verbose', this%verbose, stat, errmsg, default=.false.)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    this%verbose = (is_IOP .and. this%verbose)

    call params%get('unit', this%unit, stat, errmsg, default=-1)
    if (stat /= 0) then
      errmsg = context // errmsg
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

  subroutine fnorm(this, t, u, hdot, f, error)

    class(fht_norm), intent(inout) :: this
    real(r8), intent(in) :: t, hdot(:)
    type(fht_vector), intent(in) :: u, f
    real(r8), intent(out), optional :: error

    integer :: num1, num2
    real(r8) :: sum1, sum2, err

    ASSERT(associated(this%model%void_cell))
    ASSERT(associated(this%model%void_face))

    sum1 = global_sum(f%tc(:this%model%mesh%ncell_onP)**2, &
        mask=.not.this%model%void_cell(:this%model%mesh%ncell_onP))
    num1 = global_count(.not.this%model%void_cell(:this%model%mesh%ncell_onP))
    sum2 = global_sum(f%tf(:this%model%mesh%nface_onP)**2, &
        mask=.not.this%model%void_face(:this%model%mesh%nface_onP))
    num2 = global_count(.not.this%model%void_face(:this%model%mesh%nface_onP))

    err = sqrt((sum1+sum2)/(num1+num2))

    if (present(error)) then

      error = err / (this%abs_tol + this%rel_tol * this%err0)
      if (this%verbose) then
        write(this%unit,'(2(a,es10.3))') '  HC error: ||F||_2 =', err, ', scaled =', error
      end if

      !! Enclosure radiation system error norms
      if (this%model%vf_rad%is_active()) then
        if (this%verbose) then
          write(this%unit,'(a)',advance='no')'  ER error: scaled ||res||_2/||rhs||_2 ='
        end if
        call this%model%vf_rad%relative_residual_norm(t, u%tf, u%qrad, err)
        err = err / this%rad_tol
        if (this%verbose) write(this%unit,'(es10.3)',advance='no') err
        error = max(error, err)
        if (this%verbose) write(this%unit,*)
      end if

    else

      this%err0 = err
      if (this%verbose) then
        write(this%unit,'(a)') 'HEAT TRANSFER FUNCTION NORMS'
        write(this%unit,'(a,es10.3,a)') '  HC error: ||F||_2 =', err, ' (initial)'
      end if

    end if

  end subroutine fnorm

end module fht_norm_type
