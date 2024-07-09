!!
!! PCSR_PRECON_SSOR_TYPE
!!
!! A concrete implementation of the abstract base class PCSR_PRECON that
!! uses SSOR preconditioning.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, February 2015.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived type PCSR_PRECON_SSOR which is an
!!  extension of the abstract base class PCSR_PRECON that implements SSOR
!!  preconditioning.  See the base class comments for a description of the
!!  common type bound procedures.
!!
!!  The INIT procedure expects to find the following parameters in the
!!  TYPE(PARAMETER_LIST) argument PARAMS.  Parameters with a default value
!!  are optional; the others are required.
!!
!!    'num-sweeps' - The number of SSOR sweeps; > 0.
!!    'omega'      - The (over) relaxation factor (default 1.0); > 0.0
!!
!! PRECONDITIONER ALGORITHM
!!
!!  This preconditioner applies block-Jacobi SSOR sweeps to a vector.  The
!!  Jacobi blocks correspond to the row partioning of the matrix.  Each sweep
!!  consists of independent forward GS on each block, communication of new
!!  iterate, independent backward GS on each block, and communication of new
!!  iterate.  An over-relaxation is possible if it makes sense.  In serial
!!  this is precisely SSOR.  Note that the GS steps are performed only for
!!  the on-process rows, which are assumed to be the complete rows of the
!!  global matrix.  The off-process rows, which are ignored, exists simply
!!  to facilitate the FE-like matrix assembly.
!!

#include "f90_assert.fpp"

module pcsr_precon_ssor_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use pcsr_matrix_type
  use pcsr_precon_class
  use parameter_list_type
  use truchas_timers
  implicit none
  private

  type, extends(pcsr_precon), public :: pcsr_precon_ssor
    private
    integer  :: num_iter
    real(r8) :: omega
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type

contains

  subroutine init(this, A, params, stat, errmsg)

    class(pcsr_precon_ssor), intent(out) :: this
    type(pcsr_matrix), intent(in), target :: A
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context

    this%A => A

    context = 'processing ' // params%path() // ': '

    call params%get('num-cycles', this%num_iter, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%num_iter <= 0) then
      stat = 1
      errmsg = context // '"num-cycles" must be > 0'
      return
    end if

    call params%get('omega', this%omega, stat, errmsg, default=1.0_r8)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%omega <= 0.0_r8) then
      stat = 1
      errmsg = context // '"omega" must be > 0.0'
      return
    end if

  end subroutine init


  subroutine compute(this)
    class(pcsr_precon_ssor), intent(inout) :: this
    call start_timer('ssor-setup')
    call this%A%kdiag_init
    call stop_timer('ssor-setup')
  end subroutine compute


  subroutine apply(this, x)

    class(pcsr_precon_ssor), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: n, j, k
    real(r8) :: s, u(this%A%nrow)

    ASSERT(size(x) >= this%A%nrow_onP)

    call start_timer('ssor-solve')

    u = 0.0_r8
    do n = 1, this%num_iter
      !! Forward sweep.
      do j = 1, this%A%nrow_onP
        s = x(j)
        do k = this%A%graph%xadj(j), this%A%graph%xadj(j+1)-1
          s = s - this%A%values(k) * u(this%A%graph%adjncy(k))
        end do
        u(j) = u(j) + this%omega * (s / this%A%values(this%A%kdiag(j)))
      end do
      call this%A%graph%row_imap%gather_offp(u)
      !! Backward sweep.
      do j = this%A%nrow_onP, 1, -1
        s = x(j)
        do k = this%A%graph%xadj(j), this%A%graph%xadj(j+1)-1
          s = s - this%A%values(k) * u(this%A%graph%adjncy(k))
        end do
        u(j) = u(j) + this%omega * (s / this%A%values(this%A%kdiag(j)))
      end do
      call this%A%graph%row_imap%gather_offp(u)
    end do

    !! Copy solution to return array.
    x(:this%A%nrow_onP) = u(:this%A%nrow_onP)

    call stop_timer('ssor-solve')

  end subroutine apply

end module pcsr_precon_ssor_type
