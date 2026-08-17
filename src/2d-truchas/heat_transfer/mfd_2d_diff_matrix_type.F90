!!
!! MFD_2D_DIFF_MATRIX_TYPE
!!
!! This module defines a derived type that constructs and stores the local
!! frozen-coefficient 2D mimetic finite difference diffusion operator. It is
!! a double-sized block matrix involving primary cell-based unknowns and
!! face-based Lagrange multipliers. The type supplies matrix assembly,
!! boundary-condition modification, and Schur-complement construction for the
!! associated diffusion preconditioner.
!!
!! David Neill <davidhneill@gmail.com>, May 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! NOTES
!!
!! INIT(DISC) initializes the matrix structure from an MFD discretization.
!! INIT(MOLD) initializes an independent matrix with the same discretization
!! as MOLD, but sharing MOLD's face-matrix graph. MOLD must persist for the
!! lifetime of the matrix.
!!
!! For each preconditioner assembly, call COMPUTE(COEF) first. It defines the
!! raw diffusion matrix and discards all modifications made by the preceding
!! assembly. Next, call INCR_CELL_DIAG as needed to add time derivative terms,
!! then call SET_DIR_FACES for each distinct Dirichlet face set. Finally, call
!! COMPUTE_FACE_SCHUR_MATRIX to form the face Schur complement after all
!! matrix modifications are complete.
!!
!! SET_DIR_FACES both projects the face block and records the specified faces.
!! The record is used when forming the Schur complement and by the diffusion
!! preconditioner during elimination and back substitution. It remains valid
!! until the next call to COMPUTE.
!!

#include "f90_assert.fpp"

module mfd_2d_diff_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use mfd_2d_disc_type
  use pcsr_matrix_type
  use parallel_communication
  use index_map_type
  implicit none
  private

  type, public :: mfd_2d_diff_matrix
    type(mfd_2d_disc),   pointer :: disc => null() ! unowned reference
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    real(r8), allocatable :: a11(:)     ! the cell-cell submatrix
    real(r8), allocatable :: a12_val(:) ! the cell-face submatrix
    type(pcsr_matrix)     :: a22        ! the face-face submatrix
    integer, allocatable  :: dir_faces(:)
  contains
    procedure, private :: init_disc
    procedure, private :: init_mold
    generic   :: init => init_disc, init_mold
    procedure :: compute
    procedure :: set_dir_faces
    procedure :: incr_cell_diag
    procedure :: compute_face_schur_matrix
  end type mfd_2d_diff_matrix

contains

  subroutine init_disc(this, disc)

    class(mfd_2d_diff_matrix), intent(out) :: this
    type(mfd_2d_disc), intent(in), target :: disc

    integer :: j
    type(pcsr_graph), pointer :: g
    type(index_map), pointer :: row_imap

    this%disc => disc
    this%mesh => disc%mesh

    allocate(this%a11(this%mesh%ncell))
    allocate(this%a12_val(size(this%mesh%cface)))

    !! Create a CSR matrix graph for the A22 submatrix.
    allocate(g)
    row_imap => this%mesh%face_imap
    call g%init(row_imap)
      do j = 1, this%mesh%ncell
        associate (cface => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1))
          call g%add_clique(cface)
        end associate
      end do
    do j = 1, this%mesh%nlink
      call g%add_clique(this%mesh%lface(:,j))
    end do
    call g%add_complete

    !! Create the A22 submatrix.
    call this%a22%init(g, take_graph=.true.)

  end subroutine init_disc


  subroutine init_mold(this, mold)
    class(mfd_2d_diff_matrix), intent(out) :: this
    class(mfd_2d_diff_matrix), intent(in)  :: mold
    this%disc => mold%disc
    this%mesh => mold%mesh
    allocate(this%a11(size(mold%a11)))
    allocate(this%a12_val(size(mold%a12_val)))
    call this%a22%init(mold%a22)
  end subroutine


  subroutine compute(this, coef)

    use upper_packed_matrix_procs, only: upm_col_sum

    class(mfd_2d_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: coef(:)

    integer :: j, l, ir, ic, n, nface_max
    real(r8), allocatable :: w(:)

    ASSERT(size(coef) == this%mesh%ncell)

    if (this%mesh%ncell == 0) then
      nface_max = 0
    else
      nface_max = maxval(this%mesh%cstart(2:) - this%mesh%cstart(:this%mesh%ncell))
    end if
    allocate(w(nface_max))

    call this%a22%set_all(0.0_r8)

    do j = 1, this%mesh%ncell
      associate (a12 => this%a12_val(this%mesh%cstart(j):this%mesh%cstart(j+1)-1), &
                 minv => this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1))
        if (coef(j) == 0.0_r8) then
          this%a11(j) = 0.0_r8
          a12 = 0.0_r8
          cycle
        end if
        !! Fill the A11 and A12 submatrices
        n = size(a12)
        call upm_col_sum(minv, w(:n))
        w(:n) = coef(j) * w(:n)
        this%a11(j) = sum(w(:n))
        a12 = -w(:n)
      end associate
      !! Assemble the A22 CSR submatrix.
      associate (index => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1), &
                 minv => this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1))
        l = 1
        do ic = 1, size(index)
          do ir = 1, ic-1
            call this%a22%add_to(index(ir), index(ic), coef(j)*minv(l))
            call this%a22%add_to(index(ic), index(ir), coef(j)*minv(l))
            l = l + 1
          end do
          call this%a22%add_to(index(ic), index(ic), coef(j)*minv(l))
          l = l + 1
        end do
      end associate
    end do

    if (allocated(this%dir_faces)) deallocate(this%dir_faces)

  end subroutine compute

 !! Set the specified faces to be Dirichlet faces. These are added to the
 !! existing set of Dirichlet faces, if any. Unknowns associated with these
 !! faces are not actually unknowns and should not be included in the diffusion
 !! matrix. This is handled by projecting out the corresponding face rows and
 !! columns and setting unit diagonal values for those faces, effectively
 !! replacing equations for those unknowns with dummy equations that are
 !! decoupled from all other unknowns. The face-face submatrix is directly
 !! modified. The modified cell-face submatrix is kept in factored form as
 !! the product of the original cell-face submatrix and the projection matrix
 !! described by the list of Dirichlet faces; this is more efficient than
 !! directly modifying the cell-face submatrix.

  subroutine set_dir_faces(this, dir_faces)

    class(mfd_2d_diff_matrix), intent(inout) :: this
    integer, intent(in) :: dir_faces(:)

    integer :: j, n
    integer, allocatable :: tmp(:)

    ASSERT(minval(dir_faces) >= 1)
    ASSERT(maxval(dir_faces) <= this%mesh%nface)

    if (allocated(this%dir_faces)) then
      if (size(dir_faces) > 0) then
        n = size(this%dir_faces) + size(dir_faces)
        allocate(tmp(n))
        n = size(this%dir_faces)
        tmp(:n) = this%dir_faces
        tmp(n+1:) = dir_faces
        deallocate(this%dir_faces)
        call move_alloc(tmp, this%dir_faces)
      end if
    else
      allocate(this%dir_faces(size(dir_faces)))
      this%dir_faces = dir_faces
    end if

    do j = 1, size(dir_faces)
      n = dir_faces(j)
      call this%a22%project_out(n)
      call this%a22%set(n, n, 1.0_r8)
    end do

  end subroutine set_dir_faces

  !! Increment the (entire) diagonal cell-cell diffusion submatrix with
  !! the specified values.  The intended use is the incorporation of time
  !! derivative terms into the base diffusion matrix.

  subroutine incr_cell_diag(this, values)
    class(mfd_2d_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    ASSERT(size(values) == size(this%a11))
    this%a11 = this%a11 + values
  end subroutine

  !! Compute the face Schur complement matrix SFF. SFF is intent(inout). It
  !! must be defined on input (structure only), and its values are overwritten
  !! with the values of the Schur complement. Its structure must match the
  !! structure of the face-face diffusion submatrix (same PCSR_GRAPH component).

  subroutine compute_face_schur_matrix(this, Sff)

    class(mfd_2d_diff_matrix), intent(in) :: this
    type(pcsr_matrix), intent(inout) :: Sff

    integer :: j, n, ir, ic
    real(r8) :: tmp

    ASSERT(associated(this%a22%graph, Sff%graph))

    Sff%values = this%a22%values
    do j = 1, this%mesh%ncell
      associate (indices => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1), &
                   a12 => this%a12_val(this%mesh%cstart(j):this%mesh%cstart(j+1)-1))
        do ir = 1, size(indices)
          do ic = 1, size(indices)
            tmp = -a12(ir)*a12(ic)/this%a11(j)
            call Sff%add_to(indices(ir), indices(ic), tmp)
          end do
        end do
      end associate
    end do

    !! Apply the Dirichlet projections.
    if (allocated(this%dir_faces)) then
      do j = 1, size(this%dir_faces)
        n = this%dir_faces(j)
        call Sff%project_out(n)
        call Sff%set(n, n, 1.0_r8)
      end do
    end if

  end subroutine compute_face_schur_matrix

end module mfd_2d_diff_matrix_type
