!!
!! MFD_DIFF_MATRIX_TYPE
!!
!! This module provides a derived type that constructs and stores the local
!! mimetic finite difference diffusion matrix. It is a double-sized matrix
!! involving both the primary cell-based unknowns and the face-based Lagrange
!! multiplier unknowns.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! INIT(DISC) or INIT(MOLD) initializes the object. In the first form DISC is
!! a MFD_DISC object that defines the MFD discretization over an unstructured
!! mesh. In the second form MOLD is another MFD_DIFF_MATRIX object. In this
!! case the object is initialized using the same underlying MFD_DISC object
!! but it uses a reference to the face adjacency graph of the MOLD object
!! rather constructing its own independent graph.
!!
!! COMPUTE(D) computes the values of the matrix given the cell-based array D
!! diffusion coefficients. This is the raw matrix assuming natural boundary
!! conditions. Following methods are used to modify the values to account for
!! boundary conditions.
!!
!! SET_DIR_FACES(DIR_FACES) imposes Dirichlet boundary conditions on the faces
!! specified in the rank-1 integer array DIR_FACES. This can be called multiple
!! times (with distinct lists of faces). Dirichlet faces must be respecified
!! after a call to COMPUTE.
!!
!! INCR_CELL_DIAG(VALUES) increments the diagonal elements of the diagonal cell
!! block with the specified values. This is for incorporating a time derivative
!! term into the matrix, for example. VALUES is a complete cell-based vector.
!!
!! INCR_FACE_DIAG(INDICES, VALUES) increments the specified diagonal elements
!! of the diagonal face block with the corresponding values. This is used for
!! the discretization of certain flux-type boundary conditions.
!!
!! INCR_INTERFACE_FLUX(FACES, VALUES)
!! INCR_INTERFACE_FLUX2(FACES, VALUES)
!! INCR_INTERFACE_FLUX3(FACES, VALUES)
!! These increment specified components of the diagonal face block to implement
!! the discretization of several types of internal interface conditions.
!!
!! COMPUTE_FACE_SCHUR_MATRIX(SFF) computes the elements of the face-based Schur
!! complement matrix obtained by eliminating the cell unknowns of the diffusion
!! matrix. SFF is a PCSR_MATRIX object whose structure must have previously been
!! defined (it is the same as the face diagonal block of the diffusion matrix).
!!

#include "f90_assert.fpp"

module mfd_diff_matrix_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use mfd_disc_type
  use pcsr_matrix_type
  use parallel_communication
  use index_partitioning
  implicit none
  private

  type, public :: mfd_diff_matrix
    type(mfd_disc),   pointer :: disc => null()  ! reference only -- do not own
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    real(r8), allocatable :: a11(:)     ! the cell-cell submatrix
    real(r8), allocatable :: a12(:,:)   ! the cell-face submatrix
    type(pcsr_matrix)     :: a22        ! the face-face submatrix
    integer, allocatable  :: dir_faces(:)
    !! replacement for a12 when using unstr_mesh
    !! a12(:,j) --> a12_val(mesh%xcface(j):mesh%xcface(j+1)-1)
    real(r8), allocatable :: a12_val(:)
  contains
    procedure, private :: init_disc
    procedure, private :: init_mold
    generic   :: init => init_disc, init_mold
    procedure :: compute
    procedure :: set_dir_faces
    procedure :: incr_cell_diag
    procedure :: incr_face_diag
    procedure :: incr_interface_flux
    procedure :: incr_interface_flux2
    procedure :: incr_interface_flux3
    procedure :: compute_face_schur_matrix
  end type mfd_diff_matrix

contains

  subroutine init_disc(this, disc)

    class(mfd_diff_matrix), intent(out) :: this
    type(mfd_disc), intent(in), target :: disc

    integer :: j
    type(pcsr_graph), pointer :: g
    type(ip_desc), pointer :: row_ip

    this%disc => disc
    this%mesh => disc%mesh

    allocate(this%a11(this%mesh%ncell))
    allocate(this%a12_val(size(this%mesh%cface)))

    !! Create a CSR matrix graph for the A22 submatrix.
    allocate(g)
    row_ip => this%mesh%face_ip
    call g%init(row_ip)
      do j = 1, this%mesh%ncell
        associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
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

    class(mfd_diff_matrix), intent(out) :: this
    class(mfd_diff_matrix), intent(in)  :: mold

    type(pcsr_graph), pointer :: g

    this%disc => mold%disc
    this%mesh => mold%mesh
    allocate(this%a11(size(mold%a11)))
    allocate(this%a12_val(size(mold%a12_val)))
    g => mold%a22%graph_ptr()
    call this%a22%init(g, take_graph=.false.)

  end subroutine init_mold

  subroutine compute(this, d)

    use upper_packed_matrix, only: upm_col_sum

    class(mfd_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: d(:)

    integer :: j, l, ir, ic
    real(r8), allocatable :: w(:), minv(:)

    ASSERT(size(d) == this%mesh%ncell)

    call this%a22%set_all(0.0_r8)

    do j = 1, this%mesh%ncell
      associate(a12 => this%a12_val(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        if (d(j) == 0.0_r8) then
          this%a11(j) = 0.0_r8
          a12 = 0.0_r8
          cycle
        end if
        minv = d(j) * this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1)
        !! Fill the A11 and A12 submatrices
        allocate(w(size(a12)))
        call upm_col_sum(minv, w)
        this%a11(j) = sum(w)
        a12 = -w
        deallocate(w)
      end associate
      !! Assemble the A22 CSR submatrix.
      associate(index => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        l = 1
        do ic = 1, size(index)
          do ir = 1, ic-1
            call this%a22%add_to(index(ir), index(ic), minv(l))
            call this%a22%add_to(index(ic), index(ir), minv(l))
            l = l + 1
          end do
          call this%a22%add_to(index(ic), index(ic), minv(l))
          l = l + 1
        end do
      end associate
    end do

    if (allocated(this%dir_faces)) deallocate(this%dir_faces)

  end subroutine compute

  !! This fix-up subroutine sets selected faces to be Dirichlet faces.  The
  !! implication is that the unknowns associated with those faces are not
  !! actually unknowns and should not be included in the diffusion matrix.
  !! This is handled by projecting out the corresponding face rows and columns
  !! and setting unit diagonal values for those faces (effectively replacing
  !! equations for those unknows with dummy equations that are decoupled from
  !! all other unknowns).  The face-face submatrix is directly modified.  The
  !! modified cell-face submatrix is kept in factored form as the product of
  !! the original cell-face submatrix and the projection matrix described by
  !! the list of Dirichlet faces.  The array DIR_FACES is a list of local
  !! indices of the Dirichlet faces.

  subroutine set_dir_faces(this, dir_faces)

    class(mfd_diff_matrix), intent(inout) :: this
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

  !! This subroutine increments the (entire) diagonal cell-cell diffusion
  !! submatrix with the specified values.  The intended use is to modify
  !! the base diffusion matrix to account for the discretization of a time
  !! derivative term.

  subroutine incr_cell_diag(this, values)
    class(mfd_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    ASSERT(size(values) == size(this%a11))
    this%a11 = this%a11 + values
  end subroutine

  !! This subroutine increments selected diagonal elements of the face-face
  !! diffusion submatrix.  The array INDICES is a list of local face indices
  !! and VALUES the corresponding array of values to add to the diagonal.
  !! The intended use is to modify the base diffusion matrix to account for
  !! the discretization of certain types of flux boundary conditions.

  subroutine incr_face_diag(this, indices, values)

    class(mfd_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: indices(:)
    real(r8), intent(in) :: values(:)

    integer :: j

    ASSERT(size(indices) == size(values))
    ASSERT(minval(indices) >= 1)
    ASSERT(maxval(indices) <= this%mesh%nface)

    do j = 1, size(indices)
      call this%a22%add_to(indices(j), indices(j), values(j))
    end do

  end subroutine incr_face_diag

  !! This subroutine increments certain elements of the face-face diffusion
  !! submatrix in a specific manner to account for the contribution from the
  !! internal interface flux conditions.  For a matching pair of faces, I and J,
  !! spanning the iterface, the linearization of the condition adds to the
  !! {I,J} submatrix the 2x2 matrix [A, -A; -A, A], for some scalar A.

  subroutine incr_interface_flux(this, faces, values)

    class(mfd_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: faces(:,:)
    real(r8), intent(in) :: values(:)

    integer :: j, n1, n2

    ASSERT(size(faces,1) == 2)
    ASSERT(size(faces,2) == size(values))
    ASSERT(minval(faces) >= 1)
    ASSERT(maxval(faces) <= this%mesh%nface)

    do j = 1, size(faces,2)
      n1 = faces(1,j)
      n2 = faces(2,j)
      call this%a22%add_to(n1, n1,  values(j))
      call this%a22%add_to(n1, n2, -values(j))
      call this%a22%add_to(n2, n1, -values(j))
      call this%a22%add_to(n2, n2,  values(j))
    end do

  end subroutine incr_interface_flux

  !! This subroutine increments certain elements of the face-face diffusion
  !! submatrix in a specific manner to account for the contribution from the
  !! internal gap radiation conditions.  For a matching pair of faces, I and J,
  !! spanning the iterface, the linearization of the condition adds to the
  !! {I,J} submatrix the 2x2 matrix [A, -B; -A, B], for some scalars A and B.

  subroutine incr_interface_flux2(this, faces, values)

    class(mfd_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: faces(:,:)
    real(r8), intent(in) :: values(:,:)

    integer :: j, n1, n2

    ASSERT(size(faces,1) == 2)
    ASSERT(size(values,1) == 2)
    ASSERT(size(faces,2) == size(values,2))
    ASSERT(minval(faces) >= 1)
    ASSERT(maxval(faces) <= this%mesh%nface)

    do j = 1, size(faces,2)
      n1 = faces(1,j)
      n2 = faces(2,j)
      call this%a22%add_to(n1, n1,  values(1,j))
      call this%a22%add_to(n1, n2, -values(2,j))
      call this%a22%add_to(n2, n1, -values(1,j))
      call this%a22%add_to(n2, n2,  values(2,j))
    end do

  end subroutine incr_interface_flux2

  subroutine incr_interface_flux3(this, faces, values)

    class(mfd_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: faces(:,:)
    real(r8), intent(in) :: values(:,:)

    integer :: j, n1, n2

    ASSERT(size(faces,1) == 2)
    ASSERT(size(values,1) == 2)
    ASSERT(size(faces,2) == size(values,2))
    ASSERT(minval(faces) >= 1)
    ASSERT(maxval(faces) <= this%mesh%nface)

    do j = 1, size(faces,2)
      n1 = faces(1,j)
      n2 = faces(2,j)
      call this%a22%add_to(n1, n1,  values(1,j))
      call this%a22%add_to(n1, n2,  values(2,j))
      call this%a22%add_to(n2, n1, -values(1,j))
      call this%a22%add_to(n2, n2, -values(2,j))
    end do

  end subroutine incr_interface_flux3

  subroutine compute_face_schur_matrix (this, Sff)

    class(mfd_diff_matrix), intent(in) :: this
    type(pcsr_matrix), intent(inout) :: Sff

    integer :: j, n, ir, ic
    real(r8) :: val

    ASSERT(associated(this%a22%graph, Sff%graph))

    Sff%values = this%a22%values
    do j = 1, this%mesh%ncell
      associate (indices => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1), &
                   a12 => this%a12_val(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        do ir = 1, size(indices)
          do ic = 1, size(indices)
            val = -a12(ir)*a12(ic)/this%a11(j)
            call Sff%add_to(indices(ir), indices(ic), val)
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

end module mfd_diff_matrix_type
