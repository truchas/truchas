!!
!! DIFF_PRECON_TYPE
!!
!! This module provides a derived type and associated procedures that implement
!! a preconditioner for the local mimetic finite difference diffusion matrix.
!!
!! Neil Carlson <nnc@lanl.gov>
!!
!! The preconditioner is derived from a frozen-coefficient diffusion matrix
!! that is double-sized, involving both the primary cell-based unknowns and
!! the face-based Lagrange multiplier unknowns.  The procedure implemented
!! here first eliminates the cell unknowns leaving a face-based system.  A
!! preconditioner (SSOR, Hypre BoomerAMG, etc.) is applied to this reduced
!! Schur complement system to obtain the result for the face unknowns.  The
!! result for the cell unknowns is then obtained directly by back substitution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived type DIFF_PRECON (private components) and
!!  the following procedures that operate on instances of the type passed as
!!  the initial THIS argument.
!!
!!  CALL DIFF_PRECON_INIT (THIS, DM, PARAMS) configures the object to be a
!!    preconditioner for the  diffusion matrix DM of type DIST_DIFF_MAT.
!!    The object holds a reference to DM, and so the actual argument must
!!    either be a pointer or have the target attribute and must persist for
!!    the lifetime of the object.  DM must have been initialized, to define
!!    its structure, but its values need not have been defined at this point.
!!    Note that the matrix will not be modified in any way by this, or the
!!    other procedures.
!!
!!    PARAMS is an intent-in argument of type DIFF_PRECON_PARAMS which has
!!    the following components:
!!      SOLVER        The choice of solver to precondition the face Schur
!!                    complement system: either 'SSOR' or 'BoomerAMG'.
!!      SSOR_PARAMS   The parameters for the SSOR preconditioner.  See the
!!                    description of the PARAMS argument for SSOR_PRECON_INIT
!!                    for the components of this derived type value.
!!      BAMG_PARAMS   The parameters for the BoomerAMG preconditioner. See the
!!                    description of the PARAMS argument for BAMG_PRECON_INIT
!!                    for the components of this derived type value.
!!
!!  CALL DIFF_PRECON_COMPUTE (THIS) performs the final setup and configuration
!!    of the diffusion preconditioner.  It must be called before applying
!!    the preconditioner and after the diffusion matrix values have been set,
!!    and must be called again whenever the diffusion matrix is modified.
!!
!!  CALL DIFF_PRECON_APPLY (THIS, V1, V2) applies the diffusion preconditioner
!!    to the vector (V1, V2), where V1 is the cell-based part of the vector and
!!    V2 is the face-based part.  Both V1 and V2 are extended vectors that
!!    include off-process components with consistent values.
!!
!!  DIFF_PRECON_MATRIX(THIS) returns a pointer to the diffusion matrix DM with
!!    which the preconditioner was initialized.  This can be used to access the
!!    matrix and redefine its values.  If the values are redefined, the COMPUTE
!!    procedure must be called again before applying the updated preconditioner.
!!
!!  CALL DIFF_PRECON_DELETE (THIS) frees resources allocated by the object.
!!
!! IMPLEMENTATION NOTES
!!
!!  The params derived type is intended as a temporary stand-in for the use of
!!  a generic parameter list capability similiar to the Teuchos::ParameterList
!!  class from Trilinos.  A partial F2003 implementation exists.
!!
!!  Transitioning to F2003. The design is intended make the change to an OO
!!  implementation straightforward:
!!  * All public procedures become type-bound; delete is a final procedure.
!!  * Drop the DIFF_PRECON_ prefix from the method names and don't export
!!    them as unbound procedures.
!!  * The interface to the SSOR and BoomerAMG preconditioners are identical
!!    (except perhaps for the INIT procedure.)  That interface should be
!!    abstracted out as an abstract base type and the specific solver
!!    components replaced by a polymorphic variable of that type, eliminating
!!    the need for the all the select case constructs.
!!

#include "f90_assert.fpp"

module diff_precon_type

  use kinds
  use diffusion_matrix
  use pcsr_matrix_type
  use pcsr_precon_class
  use index_partitioning
  use unstr_mesh_type
  use parameter_list_type
  implicit none
  private

  type, public :: diff_precon
    private
    type(dist_diff_matrix), pointer :: dm => null()
    type(pcsr_matrix), pointer :: Sff => null()
    class(pcsr_precon), allocatable :: Sff_precon
  end type diff_precon

  public :: diff_precon_init, diff_precon_compute, diff_precon_delete
  public :: diff_precon_apply, diff_precon_matrix

contains

  subroutine diff_precon_init (this, dm, params)

    use pcsr_precon_ssor_type
    use pcsr_precon_boomer_type
    use parameter_list_type

    type(diff_precon), intent(out) :: this
    type(dist_diff_matrix), pointer :: dm
    type(parameter_list) :: params

    type(parameter_list), pointer :: precon_params
    character(:), allocatable :: method

    this%dm => dm ! take ownership
    allocate(this%Sff)
    call this%Sff%init (mold=dm%a22)
    !! Process the parameters.
    call params%get('method', method)
    select case (method)
    case ('SSOR')
      allocate(pcsr_precon_ssor :: this%Sff_precon)
    case ('BoomerAMG')
      allocate(pcsr_precon_boomer :: this%Sff_precon)
    case default
      INSIST(.false.)
    end select
    precon_params => params%sublist('params')
    call this%Sff_precon%init (this%Sff, precon_params)
    nullify(dm) ! object takes ownership of the matrix

  end subroutine diff_precon_init

  subroutine diff_precon_delete (this)
    type(diff_precon), intent(inout) :: this
    if (associated(this%Sff)) deallocate(this%Sff)
    if (associated(this%dm)) deallocate(this%dm)
  end subroutine diff_precon_delete

  function diff_precon_matrix (this) result (matrix)
    type(diff_precon), intent(in) :: this
    type(dist_diff_matrix), pointer :: matrix
    matrix => this%dm
  end function diff_precon_matrix

  subroutine diff_precon_compute (this)
    type(diff_precon), intent(inout) :: this
    call this%dm%compute_face_schur_matrix (this%Sff)
    call this%Sff_precon%compute
  end subroutine diff_precon_compute

  subroutine diff_precon_apply (this, f1x, f2x)

    type(diff_precon), intent(in) :: this
    real(r8), intent(inout)  :: f1x(:), f2x(:)

    ASSERT(size(f1x) == this%dm%mesh%ncell)
    ASSERT(size(f2x) == this%dm%mesh%nface)

    !! Eliminate the cell unknowns.
    call forward_elimination (this%dm, f1x, f2x)

    !! Approximately solve the Schur complement system for the face unknowns.
    call this%Sff_precon%apply (f2x)
    call gather_boundary (this%dm%mesh%face_ip, f2x)

    !! Solve for the cell unknowns by back substitution.
    call backward_substitution (this%dm, f1x, f2x)

  end subroutine diff_precon_apply

  subroutine forward_elimination (this, b1x, b2x)

    type(dist_diff_matrix), intent(in) :: this
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

    call gather_boundary (this%mesh%face_ip, b2x)

  end subroutine forward_elimination

  subroutine backward_substitution (this, b1x, u2x)

    type(dist_diff_matrix), intent(in) :: this
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

    call gather_boundary (this%mesh%cell_ip, b1x)

  end subroutine backward_substitution

end module diff_precon_type
