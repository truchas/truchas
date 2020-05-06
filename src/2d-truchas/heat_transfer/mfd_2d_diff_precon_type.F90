!TODO: finish documentation
!!
!! MFD_2D_DIFF_PRECON_TYPE
!!
!! This module defines a derived type that implements the preconditioner
!! for a 2D local mimetic finite difference diffusion matrix.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! April 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The preconditioner is derived from a frozen-coefficient diffusion matrix
!! that is double-sized, involving both the primary cell-based unknowns and
!! the face-based Lagrange multiplier unknowns.  The procedure implemented
!! here first eliminates the cell unknowns leaving a face-based Schur
!! complement system.  A preconditioner (SSOR, Hypre BoomerAMG, etc.) is
!! applied to that reduced system to obtain the result for the face unknowns.
!! The result for the cell unknowns is then obtained directly by back
!! substitution.
!!
!! The derived type MFD_2D_DIFF_PRECON has the following type bound procedures.
!!
!!  INIT(DM, PARAMS) configures the object to be a preconditioner for the
!!    specified diffusion matrix DM.  DM is an allocatable variable of type
!!    MFD_2D_DIFF_MATRIX.  The object takes ownership of DM by moving its
!!    allocation to an internal component, and DM is returned unallocated.
!!    Use the function MATRIX to get a pointer to this internal component.
!!    DM must have been initialized to define its structure, but its values
!!    need not have been defined at this point. Note that the matrix will not
!!    be modified in any way by this, or the other procedures.
!!
!!    The PARAMETER_LIST type argument PARAMS gives the parameters associated
!!    with the preconditioner.  It expects two parameters:
!!
!!      'method'      - The choice of preconditioner for the face Schur
!!                      complement system.  Valid values are either
!!                      'SSOR', or 'BoomerAMG'
!!      'parameters'  - This is a sublist of the parameters for the
!!                      specified preconditioner.  See the comments for
!!                      CSR_PRECON_SSOR_TYPE and CSR_PRECON_BOOMER_TYPE
!!                      for a description.  This sublist is passed to
!!                      the INIT subroutine for those types.
!!
!!  COMPUTE() performs the final setup and configuration of the preconditioner.
!!    It must be called before calling APPLY and after the diffusion matrix
!!    values have been set, and must be called again whenever the matrix is modified.
!!
!!  APPLY(V1, V2) applies the diffusion preconditioner to the vector (V1, V2),
!!    where V1 and V2 are the cell and face-based parts of the vector.
!!
!!  MATRIX() returns a pointer to the MFD_2D_DIFF_MATRIX object DM with which the
!!    the preconditioner was initialized.  This can be used to access the matrix
!!    and redefine its values.  If the values are redefined, the COMPUTE method
!!    must be called again before applying the updated preconditioner.
!!
!! IMPLEMENTATION NOTES
!!
!!  The Sff_precon component holds an internal reference to the Sff component.
!!  This requires Sff to have the target attribute, and the easiest way to do
!!  this is to make it a pointer -- that is the only reason it is a pointer.
!!  The alternative is to give the passed object argument to INIT the target
!!  attribute, and so forth -- a cascade up the calling chain.
!!

#include "f90_assert.fpp"

module mfd_2d_diff_precon_type

  use kinds, only: r8
  use mfd_2d_diff_matrix_type
  use pcsr_matrix_type
  use pcsr_precon_class
  use index_partitioning
  use parameter_list_type
  implicit none
  private

  type, public :: mfd_2d_diff_precon
    private
    type(mfd_2d_diff_matrix), allocatable :: dm
    type(pcsr_matrix), pointer :: Sff => null()  ! needs the target attribute
    class(pcsr_precon), allocatable :: Sff_precon
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
    procedure :: matrix
  end type mfd_2d_diff_precon

contains

  !TODO: switch to using target for Sff?
  !! Final subroutine for MFD_2D_DIFF_PRECON objects
  subroutine mfd_2d_diff_precon_delete(this)
    type(mfd_2d_diff_precon) :: this
    if (associated(this%Sff)) deallocate(this%Sff)
  end subroutine mfd_2d_diff_precon_delete

  subroutine init(this, dm, params)

    use pcsr_precon_ssor_type
    use pcsr_precon_boomer_type
    use truchas_logging_services

    class(mfd_2d_diff_precon), intent(out) :: this
    type(mfd_2d_diff_matrix), allocatable, intent(inout) :: dm
    type(parameter_list) :: params

    integer :: stat
    character(:), allocatable :: context, errmsg, method
    type(parameter_list), pointer :: plist

    !! Take the diffusion matrix.
    call move_alloc(dm, this%dm)

    !! Initialize the Schur complement matrix; same structure as A22.
    allocate(this%Sff)
    call this%Sff%init(mold=this%dm%a22)

    !! Instantiate the requested preconditioner for the Schur complement matrix.
    context = 'processing ' // params%name() // ': '
    call params%get('method', method, stat=stat, errmsg=errmsg)
    if (stat /= 0) call TLS_fatal(context//errmsg)
    select case (method)
    case ('SSOR')
      allocate(pcsr_precon_ssor :: this%Sff_precon)
    case ('BoomerAMG')
      allocate(pcsr_precon_boomer :: this%Sff_precon)
    case default
      call TLS_fatal(context//'unknown "method": '//method)
    end select
    if (params%is_sublist('parameters')) then
      plist => params%sublist('parameters')
      call this%Sff_precon%init(this%Sff, plist)
    else
      call TLS_fatal(context//'missing "parameters" sublist parameter')
    end if

  end subroutine init

  function matrix(this)
    class(mfd_2d_diff_precon), intent(in), target :: this
    type(mfd_2d_diff_matrix), pointer :: matrix
    matrix => this%dm
  end function matrix

  subroutine compute(this)
    class(mfd_2d_diff_precon), intent(inout) :: this
    call this%dm%compute_face_schur_matrix(this%Sff)
    call this%Sff_precon%compute
  end subroutine compute

  subroutine apply(this, f1, f2)
    class(mfd_2d_diff_precon), intent(in) :: this
    real(r8), intent(inout) :: f1(:), f2(:)
    ASSERT(size(f1) == this%dm%mesh%ncell)
    ASSERT(size(f2) == this%dm%mesh%nface)
    !! Eliminate the cell unknowns.
    call this%dm%forward_elimination(f1, f2)
    !! Approximately solve the Schur complement system for the face unknowns.
    call this%Sff_precon%apply(f2)
    call gather_boundary(this%dm%mesh%face_ip, f2)
    !! Solve for the cell unknowns by back substitution.
    call this%dm%backward_substitution(f1, f2)
  end subroutine apply

end module mfd_2d_diff_precon_type
