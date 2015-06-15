!!
!! MFD_DISC_TYPE
!!
!! This module defines a derived type and associated procedures that provide
!! some core functionality for implementing mimetic finite difference (MFD)
!! discretizations of differential operators such as the diffusion operator.
!!
!! Neil Carlson <nnc@lanl.gov>
!! Revised April 2015
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived type MFD_DISC which composes a reference
!!  to the computational mesh with additional persistent mesh-dependent data
!!  that is specfic to MFD.  This type is intended to be used by trusted code
!!  and so the data components are public.  However, the components should
!!  be regarded as read-only, just as the mesh, because it is also expected
!!  that an instance of this type will be referenced by multiple higher-level
!!  objects.  The available components are:
!!
!!    MESH -- a pointer to the underlying computational mesh
!!    MINV -- the inverse of the local flux mass matrices.  MINV(:,j) is the
!!      inverse of the (symmetric) local flux mass matrix on cell j in upper
!!      packed matrix format.
!!
!!  The MFD_DISC type has the following type bound procedures.
!!
!!  INIT(MESH, MINV) initializes the object with MESH as the underlying
!!    computational mesh.  The object holds a reference to MESH, and so the
!!    actual argument must either be a pointer or have the target attribute,
!!    and must persist for the lifetime of the object.  The MINV component
!!    of the object is defined by this call unless the optional logical
!!    argument MINV is specified with the value false.
!!
!!  APPLY_DIFF(COEF, UCELL, UFACE, RCELL, RFACE) applies the local MFD
!!    diffusion operator to a vector.  The vector is partitioned into cell
!!    and face-based parts in the UCELL and UFACE vectors, and the like-
!!    partitioned result is returned in the RCELL and RFACE vectors.  COEF
!!    is the cell-based vector of diffusion coefficients.  Note that this
!!    procedure does not invoke any communication, so that the off-process
!!    face elements of the result are not correct.  The off-process cell
!!    elements are correct, however.
!!
!!  COMPUTE_CELL_GRAD(UFACE, GRAD [,MASK]) computes an approximation to the
!!    average gradient of a field u on each cell given average values of u
!!    on the faces of the mesh.  The face values are passed in the array UFACE
!!    and the computed gradient in the array GRAD; GRAD(:,j) is the gradient
!!    on cell j.  Note that this ignores the average cell values of u which
!!    are also available (typically).  The results are really only intended
!!    to be used for output, and may not be at all suitable for use in a
!!    discretization scheme.
!!
!! NOTES
!!
!!  NNC, Aug 2014.  This module has grown unruly and needs to be re-thought.
!!  It is essentially a collection of (type bound) procedures requiring no
!!  state data other than the mesh pointer, plus other procedures that do
!!  expect state (MINV) that is costly (space and time).  In addition, the
!!  data components of the type are public, so the object cannot automatically
!!  generate that state itself on demand.  Thus the addition of the optional
!!  argument to the initialization call.  Users of this type that know they
!!  will not need that state data can opt to not create it; the default is
!!  to create it.
!!

#include "f90_assert.fpp"

module mfd_disc_type

  use kinds
  use base_mesh_class
  use dist_mesh_type
  implicit none
  private

  type, public :: mfd_disc
    class(base_mesh), pointer :: mesh => null()  ! reference only - do not own
    real(r8), allocatable :: minv(:,:)
    ! new stuff for mixed cell type meshes
  contains
    procedure :: init => mfd_disc_init
    procedure :: apply_diff => mfd_disc_apply_diff
    procedure :: compute_flux_matrix => mfd_disc_compute_flux_matrix
    procedure, private :: mfd_disc_compute_cell_grad1
    procedure, private :: mfd_disc_compute_cell_grad2
    generic :: compute_cell_grad => mfd_disc_compute_cell_grad1, mfd_disc_compute_cell_grad2
  end type mfd_disc

  !! Private types for internal use.
  type :: mfd_hex
    real(r8) :: volume
    real(r8) :: corner_volumes(8)
    real(r8) :: face_normals(3,6)
  contains
    procedure :: init => mfd_hex_init
    procedure :: compute_flux_matrix => mfd_hex_compute_flux_matrix
  end type mfd_hex

  type :: mfd_tet
    real(r8) :: volume
    real(r8) :: face_normals(3,4)
  contains
    procedure :: init => mfd_tet_init
    procedure :: compute_flux_matrix => mfd_tet_compute_flux_matrix
  end type mfd_tet

contains

  subroutine mfd_disc_init (this, mesh, minv)
    class(mfd_disc), intent(out) :: this
    class(base_mesh), intent(in), target :: mesh
    logical, intent(in), optional :: minv
    logical :: define_minv
    this%mesh => mesh
    define_minv = .true.
    if (present(minv)) define_minv = minv
    if (define_minv) call init_minv (this)
  end subroutine mfd_disc_init

  subroutine mfd_disc_apply_diff (this, coef, ucell, uface, rcell, rface)

    use upper_packed_matrix, only: sym_matmul

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in)  :: coef(:)
    real(r8), intent(in)  :: ucell(:), uface(:)
    real(r8), intent(out) :: rcell(:), rface(:)

    integer :: j
    real(r8), allocatable :: flux(:)

    ASSERT(size(coef) == this%mesh%ncell)
    ASSERT(size(ucell) == size(coef))
    ASSERT(size(rcell) == size(ucell))
    ASSERT(size(uface) == this%mesh%nface)
    ASSERT(size(rface) == size(uface))

    select type (mesh => this%mesh)
    type is (dist_mesh)
      allocate(flux(size(mesh%cface,dim=1)))
      rface = 0.0_r8
      do j = 1, mesh%ncell
        flux = coef(j) * sym_matmul(this%minv(:,j), ucell(j) - uface(mesh%cface(:,j)))
        rface(mesh%cface(:,j)) = rface(mesh%cface(:,j)) - flux
        rcell(j) = sum(flux)
      end do
    class default
      INSIST(.false.)
    end select

  end subroutine mfd_disc_apply_diff

  !! This auxillary procedure allocates and initializes the MINV component.

  subroutine init_minv (this)

    type(mfd_disc), intent(inout) :: this

    integer :: j
    type(mfd_tet), allocatable :: tet
    type(mfd_hex), allocatable :: hex

    select type (mesh => this%mesh)
    type is (dist_mesh)
      select case (mesh%cell_type)
      case ('TET')
        allocate(tet, this%minv(10,mesh%ncell))
        do j = 1, mesh%ncell
          call tet%init (mesh%x(:,mesh%cnode(:,j)))
          call tet%compute_flux_matrix (1.0_r8, this%minv(:,j), invert=.true.)
        end do
      case ('HEX')
        allocate(hex, this%minv(21,mesh%ncell))
        do j = 1, mesh%ncell
          call hex%init (mesh%x(:,mesh%cnode(:,j)))
          call hex%compute_flux_matrix (1.0_r8, this%minv(:,j), invert=.true.)
        end do
      case default
        INSIST(.false.)
      end select
    class default
      INSIST(.false.)
    end select

  end subroutine init_minv

 !! This procedure computes a cell-wise average gradient of a scalar field
 !! given its values at the mesh faces.  The average gradient for a cell C
 !! is computed using the following formula:
 !!
 !!    1/V \int_C \nabla u dx = 1/V \int_{\partial C} u n\dot dS
 !!                           \approx 1/V \sum_{f\in C} u_f n_f
 !!
 !! where V is the volume of C, u_f is the (average) value of u on face f,
 !! and n_f is the (average) outward normal to face f with magnitude equal
 !! to the area of the face.

  subroutine mfd_disc_compute_cell_grad1 (this, uface, grad)

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    real(r8), intent(out) :: grad(:,:)

    integer :: j
    type(mfd_tet), allocatable :: tet
    type(mfd_hex), allocatable :: hex

    INSIST(size(grad,1) == 3)
    INSIST(size(uface) == this%mesh%nface)
    INSIST(size(grad,2) <= this%mesh%ncell)

    select type (mesh => this%mesh)
    type is (dist_mesh)
      select case (mesh%cell_type)
      case ('TET')
        allocate(tet)
        do j = 1, size(grad,2)
          call tet%init (mesh%x(:,mesh%cnode(:,j)))
          grad(:,j) = matmul(tet%face_normals, uface(mesh%cface(:,j))) / tet%volume
        end do
      case ('HEX')
        allocate(hex)
        do j = 1, size(grad,2)
          call hex%init (mesh%x(:,mesh%cnode(:,j)))
          grad(:,j) = matmul(hex%face_normals, uface(mesh%cface(:,j))) / hex%volume
        end do
      case default
        INSIST(.false.)
      end select
    class default
      INSIST(.false.)
    end select

  end subroutine mfd_disc_compute_cell_grad1

  subroutine mfd_disc_compute_cell_grad2 (this, uface, mask, grad)

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    logical,  intent(in) :: mask(:)
    real(r8), intent(out) :: grad(:,:)

    integer :: j
    type(mfd_tet), allocatable :: tet
    type(mfd_hex), allocatable :: hex

    INSIST(size(grad,1) == 3)
    INSIST(size(uface) == this%mesh%nface)
    INSIST(size(grad,2) <= this%mesh%ncell)
    INSIST(size(mask) == size(grad,2))

    select type (mesh => this%mesh)
    type is (dist_mesh)
      select case (mesh%cell_type)
      case ('TET')
        allocate(tet)
        do j = 1, size(grad,2)
          if (mask(j)) then
            call tet%init (mesh%x(:,mesh%cnode(:,j)))
            grad(:,j) = matmul(tet%face_normals, uface(mesh%cface(:,j))) / tet%volume
          else
            grad(:,j) = 0.0_r8
          end if
        end do
      case ('HEX')
        allocate(hex)
        do j = 1, size(grad,2)
          if (mask(j)) then
            call hex%init (mesh%x(:,mesh%cnode(:,j)))
            grad(:,j) = matmul(hex%face_normals, uface(mesh%cface(:,j))) / hex%volume
          else
            grad(:,j) = 0.0_r8
          end if
        end do
      case default
        INSIST(.false.)
      end select
    class default
      INSIST(.false.)
    end select

  end subroutine mfd_disc_compute_cell_grad2

  !!!! MFD_TET type bound procedures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mfd_tet_init (this, vertices)
    use cell_geometry, only: tet_volume, tet_face_normals
    class(mfd_tet), intent(out) :: this
    real(r8), intent(in) :: vertices(:,:)
    ASSERT(size(vertices,1) == 3 .and. size(vertices,2) == 4)
    this%volume = tet_volume(vertices)
    this%face_normals = tet_face_normals(vertices)
  end subroutine mfd_tet_init

  subroutine mfd_tet_compute_flux_matrix (this, coef, matrix, invert)

    use cell_topology, only: TETRA4_VERT_FACE
    use upper_packed_matrix, only: invert_upm

    class(mfd_tet), intent(in) :: this
    real(r8), intent(in) :: coef
    real(r8), intent(out) :: matrix(:)
    logical, intent(in), optional :: invert

    integer :: c, i, j, ii, jj, loc
    real(r8) :: s, Nc(3,3), Mc(3,3)

    ASSERT(size(matrix) == 10)

    matrix = 0.0_r8
    do c = 1, 4
      Nc = this%face_normals(:,TETRA4_VERT_FACE(:,c))
      call invert_sym_3x3 (matmul(transpose(Nc),Nc), Mc)
      !! Scatter the corner matrix into the full cell flux matrix.
      !! It is essential that TETRA4_VERT_FACE(:,c) is an increasing sequence of indices.
      s = (0.25/coef)*this%volume
      do j = 1, 3
        jj = TETRA4_VERT_FACE(j,c)
        do i = 1, j
          ii = TETRA4_VERT_FACE(i,c)
          loc = ii + jj*(jj - 1)/2
          matrix(loc) = matrix(loc) + s*Mc(i,j)
        end do
      end do
    end do

    if (present(invert)) then
      if (invert) call invert_upm (matrix)
    end if

  end subroutine mfd_tet_compute_flux_matrix

  !!!! MFD_HEX type bound procedures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mfd_hex_init (this, vertices)
    use cell_geometry, only: eval_hex_volumes, hex_face_normals
    class(mfd_hex), intent(out) :: this
    real(r8), intent(in) :: vertices(:,:)
    ASSERT(size(vertices,1) == 3 .and. size(vertices,2) == 8)
    call eval_hex_volumes (vertices, this%volume, this%corner_volumes)
    this%face_normals = hex_face_normals(vertices)
  end subroutine mfd_hex_init

  subroutine mfd_hex_compute_flux_matrix (this, coef, matrix, invert)

    use cell_topology, only: HEX8_VERT_FACE
    use upper_packed_matrix, only: invert_upm

    class(mfd_hex), intent(in) :: this
    real(r8), intent(in) :: coef
    real(r8), intent(out) :: matrix(:)
    logical, intent(in), optional :: invert

    integer :: c, i, j, ii, jj, loc
    real(r8) :: s, cwt(8), Nc(3,3), Mc(3,3)

    ASSERT(size(matrix) == 21)

    matrix = 0.0_r8
    cwt = this%corner_volumes / sum(this%corner_volumes)
    do c = 1, 8
      Nc = this%face_normals(:,HEX8_VERT_FACE(:,c))
      call invert_sym_3x3(matmul(transpose(Nc),Nc), Mc)
      !! Scatter the corner matrix into the full cell flux matrix.
      !! It is essential that HEX8_VERT_FACE(:,c) is an increasing sequence of indices.
      s = this%volume * cwt(c) / coef
      do j = 1, 3
        jj = HEX8_VERT_FACE(j,c)
        do i = 1, j
          ii = HEX8_VERT_FACE(i,c)
          loc = ii + jj*(jj - 1)/2
          matrix(loc) = matrix(loc) + s*Mc(i,j)
        end do
      end do
    end do

    if (present(invert)) then
      if (invert) call invert_upm (matrix)
    end if

  end subroutine mfd_hex_compute_flux_matrix

 !! Direct inversion of a 3x3 symmetrix matrix using the formula that the
 !! inverse equals the transponse of the matrix of cofactors divided by
 !! the determinant.

  subroutine invert_sym_3x3 (A, Ainv)
    real(r8), intent(in)  :: A(3,3)
    real(r8), intent(out) :: Ainv(3,3)
    !! Transpose of the matrix of cofactors ...
    Ainv(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
    Ainv(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
    Ainv(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
    Ainv(1,2) = Ainv(2,1)
    Ainv(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    Ainv(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
    Ainv(1,3) = Ainv(3,1)
    Ainv(2,3) = Ainv(3,2)
    Ainv(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    !! and scale by the determinant to get the inverse.
    Ainv = (1.0_r8/(A(1,1)*Ainv(1,1) + A(2,1)*Ainv(2,1) + A(3,1)*Ainv(3,1))) * Ainv
  end subroutine invert_sym_3x3

  !! NNC, August 2014.  Need external access to a MFD_HEX_COMPUTE_FLUX_MATRIX.
  !! This needs to be revisited/generalized when this module is refactored.
  subroutine mfd_disc_compute_flux_matrix (this, n, matrix)
    class(mfd_disc), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(out) :: matrix(:)
    type(mfd_tet), allocatable :: tet
    type(mfd_hex), allocatable :: hex
    ASSERT(n >=1 .and. n <= this%mesh%ncell)
    select type (mesh => this%mesh)
    type is (dist_mesh)
      select case (mesh%cell_type)
      case ('TET')
        allocate(tet)
        call tet%init (mesh%x(:,mesh%cnode(:,n)))
        call tet%compute_flux_matrix (1.0_r8, matrix, invert=.false.)
      case ('HEX')
        allocate(hex)
        call hex%init (mesh%x(:,mesh%cnode(:,n)))
        call hex%compute_flux_matrix (1.0_r8, matrix, invert=.false.)
      case default
        INSIST(.false.)
      end select
    class default
      INSIST(.false.)
    end select
  end subroutine mfd_disc_compute_flux_matrix

end module mfd_disc_type
