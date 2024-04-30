!TODO: finish documentation
!! MFD_2D_DIFF_MATRIX_TYPE
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! May 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    type(mfd_2d_disc),   pointer :: disc => null()  ! reference only -- do not own
    type(unstr_2d_mesh), pointer :: mesh => null()  ! reference only -- do not own
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
    procedure :: forward_elimination
    procedure :: backward_substitution
  end type mfd_2d_diff_matrix

contains

  subroutine init_disc (this, disc)

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
          call g%add_clique (cface)
        end associate
      end do
    do j = 1, this%mesh%nlink
      call g%add_clique (this%mesh%lface(:,j))
    end do
    call g%add_complete

    !! Create the A22 submatrix.
    call this%a22%init (g, take_graph=.true.)

  end subroutine init_disc


  subroutine init_mold (this, mold)

    class(mfd_2d_diff_matrix), intent(out) :: this
    class(mfd_2d_diff_matrix), intent(in)  :: mold

    type(pcsr_graph), pointer :: g

    this%disc => mold%disc
    this%mesh => mold%mesh
    allocate(this%a11(size(mold%a11)))
    allocate(this%a12_val(size(mold%a12_val)))
    g => mold%a22%graph_ptr()
    call this%a22%init (g, take_graph=.false.)

  end subroutine init_mold


  !! Compute the diffusion matrix given the cell-based array of diffusion
  !! coefficients.  Any existing dirichlet faces are dropped
  subroutine compute (this, coef)

    use upper_packed_matrix_procs, only: upm_col_sum

    class(mfd_2d_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: coef(:)

    integer :: j, l, ir, ic
    real(r8), allocatable :: w(:), minv(:)

    ASSERT(size(coef) == this%mesh%ncell)

    call this%a22%set_all (0.0_r8)

    do j = 1, this%mesh%ncell
      associate(a12 => this%a12_val(this%mesh%cstart(j):this%mesh%cstart(j+1)-1))
        if (coef(j) == 0.0_r8) then
          this%a11(j) = 0.0_r8
          a12 = 0.0_r8
          cycle
        end if
        minv = coef(j) * this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1)
        !! Fill the A11 and A12 submatrices
        allocate(w(size(a12)))
        call upm_col_sum (minv, w)
        this%a11(j) = sum(w)
        a12 = -w
        deallocate(w)
      end associate
      !! Assemble the A22 CSR submatrix.
      minv = coef(j) * this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1)
      associate(index => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1))
        l = 1
        do ic = 1, size(index)
          do ir = 1, ic-1
            call this%a22%add_to (index(ir), index(ic), minv(l))
            call this%a22%add_to (index(ic), index(ir), minv(l))
            l = l + 1
          end do
          call this%a22%add_to (index(ic), index(ic), minv(l))
          l = l + 1
        end do
      end associate
    end do

    if (allocated(this%dir_faces)) deallocate(this%dir_faces)

  end subroutine compute


 !! Set the specified faces to be Dirichlet faces.  These are added to the
 !! existing set of Dirichlet faces, if any.  Unknowns associated with these
 !! faces are not actually unknowns and should no be included in the diffusion
 !! matrix.  This is handled by projecting out the corresponding face rows and
 !! columns and setting unit diagonal values for those faces, effectively
 !! replacing equations for those unknowns with dummy equations that are
 !! decoupled from all other unknowns.  The face-face submatrix is directly
 !! modified.  The modified cell-face submatrix is kept in factored form as
 !! the product of the original cell-face submatrix and the projection matrix
 !! described by the list of Dirichlet faces; this is more efficient than
 !! directly modifying the cell-face submatrix.

  subroutine set_dir_faces (this, dir_faces)

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
        call move_alloc (tmp, this%dir_faces)
      end if
    else
      allocate(this%dir_faces(size(dir_faces)))
      this%dir_faces = dir_faces
    end if

    do j = 1, size(dir_faces)
      n = dir_faces(j)
      call this%a22%project_out (n)
      call this%a22%set (n, n, 1.0_r8)
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
  end subroutine incr_cell_diag


  !! Compute the face Schur complement matrix SFF.  SFF is intent inout.  It
  !! must be defined on input (structure only), and its values are overwritten
  !! with the values of the Schur complement.  Its structure must match the
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
            call Sff%add_to (indices(ir), indices(ic), tmp)
          end do
        end do
      end associate
    end do

    !! Apply the Dirichlet projections.
    if (allocated(this%dir_faces)) then
      do j = 1, size(this%dir_faces)
        n = this%dir_faces(j)
        call Sff%project_out (n)
        call Sff%set (n, n, 1.0_r8)
      end do
    end if

  end subroutine compute_face_schur_matrix


  subroutine forward_elimination(this, b1, b2)

    class(mfd_2d_diff_matrix), intent(in) :: this
    real(r8), intent(in) :: b1(:)
    real(r8), intent(inout) :: b2(:)

    integer :: j
    real(r8) :: s
    real(r8), allocatable :: b2_dir(:)

    ASSERT(size(b1) == this%mesh%ncell)
    ASSERT(size(b2) == this%mesh%nface)

    if (allocated(this%dir_faces)) b2_dir = b2(this%dir_faces)

    do j = 1, this%mesh%ncell
      associate (a12 => this%a12_val(this%mesh%cstart(j):this%mesh%cstart(j+1)-1), &
                 cface => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1))
        s = b1(j) / this%a11(j)
        b2(cface) = b2(cface) - a12 * s
      end associate
    end do

    if (allocated(this%dir_faces)) b2(this%dir_faces) = b2_dir

  end subroutine forward_elimination


  subroutine backward_substitution(this, b1, u2)

    class(mfd_2d_diff_matrix), intent(in) :: this
    real(r8), intent(inout) :: b1(:), u2(:)

    integer :: j
    real(r8) :: s
    real(r8), allocatable :: u2_dir(:)

    if (allocated(this%dir_faces)) then
      u2_dir = u2(this%dir_faces)
      u2(this%dir_faces) = 0.0_r8
    end if

    do j = 1, this%mesh%ncell
      associate (a12 => this%a12_val(this%mesh%cstart(j):this%mesh%cstart(j+1)-1), &
                 cface => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1))
        s = b1(j) - sum(a12 * u2(cface))
        b1(j) = s / this%a11(j)
      end associate
    end do

    if (allocated(this%dir_faces)) u2(this%dir_faces) = u2_dir

  end subroutine backward_substitution

end module mfd_2d_diff_matrix_type
