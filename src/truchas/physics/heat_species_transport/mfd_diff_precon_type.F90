!!
!! MFD_DIFF_PRECON_TYPE
!!
!! This module provides a derived type and associated procedures that implement
!! a preconditioner for the local mimetic finite difference diffusion matrix.
!!
!! Neil Carlson <nnc@lanl.gov>
!! Updated for F2008, January 2021
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
!! here first eliminates the cell unknowns leaving a face-based system.  A
!! preconditioner (SSOR, Hypre BoomerAMG, etc.) is applied to this reduced
!! Schur complement system to obtain the result for the face unknowns.  The
!! result for the cell unknowns is then obtained directly by back substitution.
!!
!! This module defines the derived type DIFF_PRECON (private components) with
!! the following type bound procedures.
!!
!! INIT(DM, PARAMS, STAT, ERRMSG) configures the object to be a preconditioner
!!   for the diffusion matrix DM of type mfd_diff_matrix. DM is allocatable
!!   and its allocation is moved into the object. DM must have been initialized
!!   to define its structure, but its values need not have been defined at this
!!   point. Note that the matrix will not be modified in any way by this, or
!!   the other procedures. The integer STAT returns a nonzero value if an error
!!   occurs, and in such a case the allocatable character ERRMSG returns an
!!   explanatory error message. The parameter list PARAMS is expected to define
!!   the "method" parameter, which can take the values "ssor" or "boomeramg",
!!   and the "params" sublist, which defines the parameters for the selected
!!   method. See pcsr_precon_ssor_type.F90 and pcsr_precon_boomer_type.F90 for
!!   those expected parameters.
!!
!! COMPUTE() performs the final setup and configuration of the diffusion
!!   preconditioner. It must be called before applying the preconditioner and
!!   after the diffusion matrix values have been set, and must be called again
!!   whenever the diffusion matrix is modified.
!!
!! APPLY(V1, V2) applies the diffusion preconditioner to the vector (V1, V2),
!!   where V1 is the cell-based part of the vector and V2 is the face-based
!!   part.  Both V1 and V2 are extended vectors that include off-process
!!   components with consistent values.
!!
!! MATRIX_REF() returns a pointer to the diffusion matrix DM with which the
!!   preconditioner was initialized. This can be used to access the matrix and
!!   redefine its values. If the values are redefined, the COMPUTE method must
!!   be called again before applying the updated preconditioner.
!!
!! IMPLEMENTATION NOTES
!!
!!  (1) The preconditioner of the Schur complement matrix holds a reference to
!!  that matrix, which is a component of the object. To avoid having to give
!!  the object the target attribute in the INIT procedure and the consequent
!!  risk of dangling pointers up the call chain (the rules on targets are
!!  complex), the Schur complement matrix is made a pointer.
!!

#include "f90_assert.fpp"

module mfd_diff_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfd_diff_matrix_type
  use pcsr_matrix_type
  use pcsr_precon_class
  use index_partitioning
  implicit none
  private

  type, public :: mfd_diff_precon
    private
    type(mfd_diff_matrix), allocatable :: dm
    type(pcsr_matrix), pointer :: Sff => null()
    class(pcsr_precon), allocatable :: Sff_precon
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
    procedure :: matrix_ref
    final :: mfd_diff_precon_delete
  end type

contains

  !! Final subroutine for DIFF_PRECON objects
  subroutine mfd_diff_precon_delete(this)
    type(mfd_diff_precon), intent(inout) :: this
    if (associated(this%Sff)) deallocate(this%Sff)
  end subroutine

  subroutine init(this, dm, params, stat, errmsg)

    use pcsr_precon_factory
    use parameter_list_type

    class(mfd_diff_precon), intent(out), target :: this
    type(mfd_diff_matrix), allocatable, intent(inout) :: dm
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    allocate(this%Sff)
    call this%Sff%init(mold=dm%a22)
    call alloc_pcsr_precon(this%Sff_precon, this%Sff, params, stat, errmsg)
    if (stat /= 0) return
    call move_alloc(dm, this%dm)

  end subroutine init

  function matrix_ref(this) result(matrix)
    class(mfd_diff_precon), intent(in), target :: this
    type(mfd_diff_matrix), pointer :: matrix
    matrix => this%dm
  end function

  subroutine compute(this)
    class(mfd_diff_precon), intent(inout) :: this
    call this%dm%compute_face_schur_matrix(this%Sff)
    call this%Sff_precon%compute
  end subroutine

  subroutine apply(this, f1x, f2x)

    class(mfd_diff_precon), intent(in) :: this
    real(r8), intent(inout) :: f1x(:), f2x(:)

    ASSERT(size(f1x) == this%dm%mesh%ncell)
    ASSERT(size(f2x) == this%dm%mesh%nface)

    !! Eliminate the cell unknowns.
    call forward_elimination(this%dm, f1x, f2x)

    !! Approximately solve the Schur complement system for the face unknowns.
    call this%Sff_precon%apply(f2x)
    call gather_boundary(this%dm%mesh%face_ip, f2x)

    !! Solve for the cell unknowns by back substitution.
    call backward_substitution(this%dm, f1x, f2x)

  end subroutine

  subroutine forward_elimination(this, b1x, b2x)

    type(mfd_diff_matrix), intent(in) :: this
    real(r8), intent(in) :: b1x(:)
    real(r8), intent(inout) :: b2x(:)

    integer :: j, k, n
    real(r8) :: s
    real(r8), allocatable :: b2x_dir(:)

    ASSERT(size(b1x) == this%mesh%ncell)
    ASSERT(size(b2x) == this%mesh%nface)

    if (allocated(this%dir_faces)) then
      allocate(b2x_dir(size(this%dir_faces)))
      b2x_dir = b2x(this%dir_faces)
    end if

    do j = 1, this%mesh%ncell
      s = b1x(j) / this%a11(j)
      associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1), &
                 a12 => this%a12_val(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        b2x(cface) = b2x(cface) - a12 * s
      end associate
    end do

    if (allocated(this%dir_faces)) then
      b2x(this%dir_faces) = b2x_dir
      deallocate(b2x_dir)
    end if

    call gather_boundary(this%mesh%face_ip, b2x)

  end subroutine forward_elimination

  subroutine backward_substitution(this, b1x, u2x)

    type(mfd_diff_matrix), intent(in) :: this
    real(r8), intent(inout) :: b1x(:), u2x(:)

    integer :: j, k
    real(r8) :: s
    real(r8), allocatable :: u2x_dir(:)

    if (allocated(this%dir_faces)) then
      allocate(u2x_dir(size(this%dir_faces)))
      u2x_dir = u2x(this%dir_faces)
      u2x(this%dir_faces) = 0.0_r8
    end if

    do j = 1, this%mesh%ncell_onP
      associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1), &
                 a12 => this%a12_val(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        b1x(j) = (b1x(j) - dot_product(a12, u2x(cface))) / this%a11(j)
      end associate
    end do
    
    if (allocated(this%dir_faces)) then
      u2x(this%dir_faces) = u2x_dir
      deallocate(u2x_dir)
    end if

    call gather_boundary(this%mesh%cell_ip, b1x)

  end subroutine backward_substitution

end module mfd_diff_precon_type
