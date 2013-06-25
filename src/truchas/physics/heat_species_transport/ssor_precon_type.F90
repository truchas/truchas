!!
!! SSOR_PRECON_TYPE
!!
!! This module defines a derived type and associated procedures that describe
!! an SSOR preconditioner for a parallel CSR matrix.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived data type SSOR_PRECON (private components)
!!  and the following procedures that operate on instances of the type passed
!!  as the THIS argument.  An instance describes an SSOR preconditioner for a
!!  parallel CSR matrix of type PCSR_MATRIX.
!!
!!  CALL SSOR_PRECON_INIT (THIS, A, PARAMS) initializes the object to be a
!!    preconditioner for the parallel CSR matrix A of type PCSR_MATRIX.  The
!!    object holds a reference to the matrix A, and so the matrix must never
!!    go out of scope during the lifetime of the object.  Morover the actual
!!    argument must be a pointer or have the target attribute.  Only the
!!    structure of the matrix A needs to be defined at this point; the matrix
!!    values are not referenced.  Note that the matrix will not be modified
!!    in any way by this, or the other procedures.
!!
!!    PARAMS is an intent-in argument of type SSOR_PRECON_PARAMS which has
!!    the following components:
!!      NUM_ITER    The number of SSOR sweeps to apply; default 1
!!      OMEGA       The over/under relaxation factor; default 1.0
!!
!!  CALL SSOR_PRECON_COMPUTE (THIS) performs the final setup and configuration
!!    of the preconditioner.  It is at this point the values of the matrix A
!!    are referenced.  It must be called before using the APPLY procedure and
!!    after the matrix values are defined, and must be called again whenever
!!    the matrix values are modified.
!!
!!  CALL SSOR_PRECON_APPLY (THIS, X) applies the preconditioner to the vector X.
!!    The size of X must be at least A%NROW_ONP; only the the initial A%NROW_ONP
!!    elements will be referenced or modified.
!!
!!  SSOR_PRECON_MATRIX(THIS) returns a pointer to the parallel CSR matrix A
!!    with which the preconditioner was initialized.
!!
!!  CALL SSOR_PRECON_DELETE (THIS) frees resources allocated by the object.
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
!!  The params derived type is intended as a temporary stand-in for the use of
!!  a generic parameter list capability similiar to the Teuchos::ParameterList
!!  class from Trilinos.  A partial F2003 implementation exists.
!!
!!  Transitioning to F2003. The design is intended make the change to an OO
!!  implementation straightforward:
!!  * All procedures become type-bound; delete is a final procedure.
!!  * Drop the SSOR_PRECON_ prefix from the method names and don't export
!!    them as loose procedures.
!!  * The diag array component should be made allocatable.
!!  * The interface is identical in form to that for BOOMER_AMG_PRECON.
!!    Both should extend the same abstract base type, which implements
!!    the MATRIX function and defines the deferred subroutines COMPUTE
!!    and APPLY.  It could also define a deferred subroutine INIT if it
!!    were modified to accept a generic 'parameter list' argument as
!!    described above.  Otherwise the INIT method would need to be
!!    specific to a particular implementation of the base type.
!!

#include "f90_assert.fpp"

module ssor_precon_type

  use kinds, only: r8
  use parallel_csr_matrix
  implicit none
  private

  public :: ssor_precon_init
  public :: ssor_precon_delete
  public :: ssor_precon_compute
  public :: ssor_precon_apply
  public :: ssor_precon_matrix

  type, public :: ssor_precon
    private
    type(pcsr_matrix), pointer :: A => null()
    real(r8), pointer :: diag(:) => null()
    integer  :: num_iter
    real(r8) :: omega
  end type ssor_precon

  type, public :: ssor_precon_params
    integer  :: num_iter = 1
    real(r8) :: omega = 1.0_r8
  end type ssor_precon_params

contains

  subroutine ssor_precon_init (this, A, params)
    type(ssor_precon), intent(out) :: this
    type(pcsr_matrix), target, intent(in) :: A
    type(ssor_precon_params), intent(in) :: params
    this%A => A
    allocate(this%diag(A%nrow_onP))
    !! Process the parameters.
    INSIST(params%num_iter > 0)
    this%num_iter = params%num_iter
    INSIST(params%omega > 0.0_r8)
    this%omega = params%omega
  end subroutine ssor_precon_init

  function ssor_precon_matrix (this) result (matrix)
    type(ssor_precon), intent(in) :: this
    type(pcsr_matrix), pointer :: matrix
    matrix => this%A
  end function ssor_precon_matrix

  subroutine ssor_precon_delete (this)
    type(ssor_precon), intent(inout) :: this
    if (associated(this%diag)) deallocate(this%diag)
  end subroutine ssor_precon_delete

  subroutine ssor_precon_compute (this)
    type(ssor_precon), intent(inout) :: this
    integer :: i, k
    !! Ensure the diagonal array is allocated with the correct size.
    if (associated(this%diag)) then
      if (size(this%diag) /= this%A%nrow_onP) deallocate(this%diag)
    end if
    if (.not.associated(this%diag)) allocate(this%diag(this%A%nrow_onP))
    !! Extract the diagonal elements from the sparse matrix.
    do i = 1, this%A%nrow_onP
      do k = this%A%graph%xadj(i), this%A%graph%xadj(i+1)-1
        if (this%A%graph%adjncy(k) == i) then
          this%diag(i) = this%A%data(k)
          exit
        end if
      end do
    end do
  end subroutine ssor_precon_compute

  subroutine ssor_precon_apply (this, b)
    use index_partitioning, only: gather_boundary
    type(ssor_precon), intent(in) :: this
    real(r8), intent(inout) :: b(:)
    integer :: n, j, k
    real(r8) :: s, u(this%A%nrow)
    ASSERT(size(b) >= this%A%nrow_onP)
    u = 0.0_r8
    do n = 1, this%num_iter
      !! Forward sweep.
      do j = 1, this%A%nrow_onP
        s = b(j)
        do k = this%A%graph%xadj(j), this%A%graph%xadj(j+1)-1
          s = s - this%A%data(k) * u(this%A%graph%adjncy(k))
        end do
        u(j) = u(j) + this%omega * (s / this%diag(j))
      end do
      call gather_boundary (this%A%graph%row_ip, u)
      !! Backward sweep.
      do j = this%A%nrow_onP, 1, -1
        s = b(j)
        do k = this%A%graph%xadj(j), this%A%graph%xadj(j+1)-1
          s = s - this%A%data(k) * u(this%A%graph%adjncy(k))
        end do
        u(j) = u(j) + this%omega * (s / this%diag(j))
      end do
      call gather_boundary (this%A%graph%row_ip, u)
    end do
    !! Copy solution to return array.
    b(:this%A%nrow_onP) = u(:this%A%nrow_onP)
  end subroutine ssor_precon_apply

end module ssor_precon_type
