!!
!! MATL_PROP_CLASS
!!
!! This module defines the abstract base class MATL_PROP that defines an
!! interface for evaluating a general material property and its derivatives
!! with respect to one or more state variables.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!TODO: For efficiency we may want to add versions that act on a vector of states
!      and combine into a generic interface.

#include "f90_assert.fpp"

module matl_prop_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: matl_prop
  contains
    procedure(compute_value), deferred :: compute_value
    procedure :: compute_deriv
    procedure :: write_plotfile
  end type

  abstract interface
    subroutine compute_value(this, state, value)
      import matl_prop, r8
      class(matl_prop), intent(in) :: this
      real(r8), intent(in)  :: state(:)
      real(r8), intent(out) :: value
    end subroutine
  end interface

contains

  subroutine compute_deriv(this, state, n, deriv)

    class(matl_prop), intent(in) :: this
    real(r8), intent(in) :: state(:)
    integer,  intent(in) :: n
    real(r8), intent(out) :: deriv

    real(r8) :: fdinc, v1, v2, pstate(size(state))

    ASSERT(n >= 1 .and. n <= size(state))

    !TODO: make the computation of fdinc more efficient
    fdinc = max(1.0_r8, abs(state(n))) * sqrt(epsilon(1.0_r8))
    fdinc = scale(1.0_r8,exponent(fdinc))
    pstate = state
    pstate(n) = state(n) + fdinc
    call this%compute_value(pstate, v2)
    pstate(n) = state(n) - fdinc
    call this%compute_value(pstate, v1)
    deriv = (v2 - v1) / (2*fdinc)

  end subroutine compute_deriv

  subroutine write_plotfile(this, filename, digits, state_min, state_max, npoints, iostat)

    class(matl_prop), intent(in) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: digits, npoints(:)
    real(r8), intent(in) :: state_min(:), state_max(:)
    integer, intent(out) :: iostat

    integer :: narg, lun, j
    character(12) :: fmt
    real(r8), allocatable :: ds(:), s(:)
    real(r8) :: p, dpds

    narg = size(state_min)
    ASSERT(size(state_max) == narg)
    ASSERT(size(npoints) == narg)
    ASSERT(all(npoints > 0))
    ASSERT(digits > 1)

    open(newunit=lun,file=filename,status='replace',action='write',iostat=iostat)
    if (iostat /= 0) return
    write(fmt,'("(*(es",i0,".",i0,"))")') digits+7, digits-1

    ds = (state_max - state_min) / npoints
    allocate(s, mold=ds)

    select case (narg)
    case (1)

      do j = 0, npoints(1)
        s(1) = state_min(1) + j*ds(1)
        call this%compute_value(s, p)
        call this%compute_deriv(s, 1, dpds)
        write(lun,fmt) s, p, dpds
      end do

    case default
      INSIST(.false.) ! not yet implemented
    end select

  end subroutine write_plotfile

end module matl_prop_class

