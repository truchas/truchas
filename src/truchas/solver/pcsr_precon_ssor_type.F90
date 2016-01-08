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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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
!! IMPLEMENTATION NOTES
!!
!!  The INIT procedure includes some commented-out code from the Pececillo
!!  mini-app for error checking the parameter list.  This code, especially
!!  its error messages, make little sense in the current context, but are
!!  included for future reference.  The provided parameter list is created
!!  by the client code using namelist data read from the input file, and it
!!  is expected that the parameter list is fully vetted at that point.  In
!!  the future I expect the parameter list to be read directly from the input
!!  file, and when that happens this included code may become apropos.
!!

#include "f90_assert.fpp"

module pcsr_precon_ssor_type

  use kinds, only: r8
  use pcsr_matrix_type
  use pcsr_precon_class
  use parameter_list_type
  use timing_tree
  implicit none
  private

  type, extends(pcsr_precon), public :: pcsr_precon_ssor
    private
    real(r8), allocatable :: diag(:)
    integer  :: num_iter
    real(r8) :: omega
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type pcsr_precon_ssor

contains

  subroutine init (this, A, params)

    class(pcsr_precon_ssor), intent(out) :: this
    type(pcsr_matrix), intent(in), target :: A
    type(parameter_list) :: params

    this%A => A
    allocate(this%diag(A%nrow_onP))

    call params%get ('num-sweeps', this%num_iter)
    INSIST(this%num_iter > 0)
    call params%get ('omega', this%omega, default=1.0_r8)
    INSIST(this%omega > 0.0_r8)

!NNC    !! Process the parameters.
!NNC    use truchas_logging_services
!NNC    integer :: stat
!NNC    character(:), allocatable :: context, errmsg
!NNC    context = 'processing ' // params%name() // ': '
!NNC    call params%get ('num-sweeps', this%num_iter, stat=stat, errmsg=errmsg)
!NNC    if (stat /= 0) call TLS_fatal (context//errmsg)
!NNC    if (this%num_iter <= 0) call TLS_fatal (context//'"num-sweeps" must be > 0')
!NNC    call params%get ('omega', this%omega, default=1.0_r8, stat=stat, errmsg=errmsg)
!NNC    if (stat /= 0) call TLS_fatal (context//errmsg)
!NNC    if (this%omega <= 0.0_r8) call TLS_fatal (context//'"omega" must be > 0.0')

  end subroutine init


  subroutine compute (this)
    class(pcsr_precon_ssor), intent(inout) :: this
    call start_timer ('ssor-setup')
    call this%A%get_diag_copy (this%diag)
    call stop_timer ('ssor-setup')
  end subroutine compute


  subroutine apply (this, x)

    use index_partitioning, only: gather_boundary

    class(pcsr_precon_ssor), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: n, j, k
    real(r8) :: s, u(this%A%nrow)

    ASSERT(size(x) >= this%A%nrow_onP)

    call start_timer ('ssor-solve')

    u = 0.0_r8
    do n = 1, this%num_iter
      !! Forward sweep.
      do j = 1, this%A%nrow_onP
        s = x(j)
        do k = this%A%graph%xadj(j), this%A%graph%xadj(j+1)-1
          s = s - this%A%values(k) * u(this%A%graph%adjncy(k))
        end do
        u(j) = u(j) + this%omega * (s / this%diag(j))
      end do
      call gather_boundary (this%A%graph%row_ip, u)
      !! Backward sweep.
      do j = this%A%nrow_onP, 1, -1
        s = x(j)
        do k = this%A%graph%xadj(j), this%A%graph%xadj(j+1)-1
          s = s - this%A%values(k) * u(this%A%graph%adjncy(k))
        end do
        u(j) = u(j) + this%omega * (s / this%diag(j))
      end do
      call gather_boundary (this%A%graph%row_ip, u)
    end do

    !! Copy solution to return array.
    x(:this%A%nrow_onP) = u(:this%A%nrow_onP)

    call stop_timer ('ssor-solve')

  end subroutine apply

end module pcsr_precon_ssor_type
