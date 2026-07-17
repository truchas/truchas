!!
!! MFD_DISC_TYPE
!!
!! This module defines the MFD_DISC type used by thermal/species diffusion
!! models to represent a mimetic finite difference (MFD) discretization on an
!! unstructured mesh. Its primary state is the cell-local inverse flux mass
!! matrix MINV, stored in upper packed format. The standard construction
!! follows the generic polyhedral-cell MFD method described in:
!!
!!   Konstantin Lipnikov, Gianmarco Manzini, and Mikhail Shashkov. Mimetic
!!   finite difference method. Journal of Computational Physics, 257:1163-1227,
!!   2014.
!!
!!   V. Gyrya, K. Lipnikov, G. Manzini, and D. Svyatskiy. M-adaptation in the
!!   mimetic finite difference method. Mathematical Models & Methods in Applied
!!   Sciences, 24(8):1621-1663, 2014.
!!
!! The type is used by trusted code; its public data components should be
!! treated as read-only after initialization. It treats its mesh as a local
!! serial mesh, with no ownership or on/off-process semantics. Callers are
!! responsible for ensuring input values are valid on the full local mesh
!! extent they pass in, and for interpreting which entries of any result are
!! meaningful for their distributed algorithm.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module mfd_disc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  implicit none
  private

  integer, parameter, public :: MFD_STANDARD = 1
  integer, parameter, public :: MFD_LEGACY = 2
  integer, parameter :: MFD_CELL_NFACE_MAX = 6

  type, public :: mfd_disc
    type(unstr_mesh), pointer :: mesh => null()  ! reference only - do not own
    integer, allocatable :: xminv(:)
    real(r8), allocatable :: minv(:)
    integer :: method = MFD_STANDARD
  contains
    procedure :: init
    procedure :: apply_diff
    procedure :: add_diff
    generic :: compute_cell_grad => compute_cell_grad1, compute_cell_grad2
    procedure, private :: compute_cell_grad1
    procedure, private :: compute_cell_grad2
  end type

  !! Private types for internal use.
  type :: mfd_hex
    real(r8) :: volume
    real(r8) :: corner_volumes(8)
    real(r8) :: face_normals(3,6)
  contains
    procedure :: init => mfd_hex_init
    procedure :: compute_flux_matrix => mfd_hex_compute_flux_matrix
  end type

  type :: mfd_wedge
    real(r8) :: volume
    real(r8) :: corner_volumes(6)
    real(r8) :: face_normals(3,5)
  contains
    procedure :: init => mfd_wedge_init
    procedure :: compute_flux_matrix => mfd_wedge_compute_flux_matrix
  end type

  type :: mfd_tet
    real(r8) :: volume
    real(r8) :: face_normals(3,4)
  contains
    procedure :: init => mfd_tet_init
    procedure :: compute_flux_matrix => mfd_tet_compute_flux_matrix
  end type
  
  type, public :: mfd_cell
    integer :: nfaces = 0
    real(r8) :: face_centers(3,MFD_CELL_NFACE_MAX)
    real(r8) :: cell_center(3)
    real(r8) :: face_area(MFD_CELL_NFACE_MAX)
    real(r8) :: face_normals(3,MFD_CELL_NFACE_MAX)
    real(r8) :: volume
  contains
    procedure :: init => mfd_cell_init
    procedure :: compute_flux_matrix_inv => mfd_cell_compute_flux_matrix_inv
    procedure :: dump => mfd_cell_dump
  end type

contains

  !! Initialize the discretization for a mesh. The object holds an unowned
  !! reference to MESH, which must persist for the object's lifetime. When MINV
  !! is absent or true, the persistent inverse flux mass matrices are allocated
  !! and initialized. METHOD defaults to the standard direct construction.
  !! MFD_LEGACY is retained only for one-off developer comparisons.

  subroutine init(this, mesh, method, minv)
    class(mfd_disc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in), optional :: method
    logical, intent(in), optional :: minv
    logical :: define_minv
    this%mesh => mesh
    if (present(method)) this%method = method
    ASSERT(this%method == MFD_STANDARD .or. this%method == MFD_LEGACY)
    define_minv = .true.
    if (present(minv)) define_minv = minv
    if (define_minv) call init_minv (this)
  end subroutine

  !! Apply the local MFD diffusion operator as a serial operator over the full
  !! mesh held by this object. No communication is performed; COEF, UCELL, and
  !! UFACE must be valid over their full local extents. RCELL and RFACE are
  !! filled over those same extents, and callers decide which entries to use.

  subroutine apply_diff(this, coef, ucell, uface, rcell, rface)

    use upper_packed_matrix_procs, only: upm_sym_matvec

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in)  :: coef(:)
    real(r8), intent(in)  :: ucell(:), uface(:)
    real(r8), intent(out) :: rcell(:), rface(:)

    integer :: j, n, nface_max
    real(r8), allocatable :: du(:), flux(:)

    ASSERT(size(coef) == this%mesh%ncell)
    ASSERT(size(ucell) == size(coef))
    ASSERT(size(rcell) == size(ucell))
    ASSERT(size(uface) == this%mesh%nface)
    ASSERT(size(rface) == size(uface))

    nface_max = max_cell_faces(this%mesh)
    allocate(du(nface_max), flux(nface_max))
    rface = 0.0_r8
    do j = 1, this%mesh%ncell
      associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1), &
                 minv  => this%minv(this%xminv(j):this%xminv(j+1)-1))
        n = size(cface)
        du(:n) = ucell(j) - uface(cface)
        call upm_sym_matvec(minv, du(:n), flux(:n))
        flux(:n) = coef(j) * flux(:n)
        rface(cface) = rface(cface) - flux(:n)
        rcell(j) = sum(flux(:n))
      end associate
    end do

  end subroutine apply_diff

  !! Add the local MFD diffusion operator result into RCELL and RFACE. This
  !! has the same input contract as APPLY_DIFF, but accumulates into the output
  !! arrays instead of overwriting them.

  subroutine add_diff(this, coef, ucell, uface, rcell, rface)

    use upper_packed_matrix_procs, only: upm_sym_matvec

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in)    :: coef(:)
    real(r8), intent(in)    :: ucell(:), uface(:)
    real(r8), intent(inout) :: rcell(:), rface(:)

    integer :: j, n, nface_max
    real(r8), allocatable :: du(:), flux(:)

    ASSERT(size(coef) == this%mesh%ncell)
    ASSERT(size(ucell) == size(coef))
    ASSERT(size(rcell) == size(ucell))
    ASSERT(size(uface) == this%mesh%nface)
    ASSERT(size(rface) == size(uface))

    nface_max = max_cell_faces(this%mesh)
    allocate(du(nface_max), flux(nface_max))
    do j = 1, this%mesh%ncell
      associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1), &
                 minv  => this%minv(this%xminv(j):this%xminv(j+1)-1))
        n = size(cface)
        du(:n) = ucell(j) - uface(cface)
        call upm_sym_matvec(minv, du(:n), flux(:n))
        flux(:n) = coef(j) * flux(:n)
        rface(cface) = rface(cface) - flux(:n)
        rcell(j) = rcell(j) + sum(flux(:n))
      end associate
    end do

  end subroutine add_diff

  !! This auxiliary procedure allocates and initializes the MINV component.

  subroutine init_minv(this)

    use cell_topology, only: num_cell_faces
    type(mfd_disc), intent(inout) :: this

    integer :: j, n

    allocate(this%xminv(this%mesh%ncell+1))
    this%xminv(1) = 1
    do j = 1, this%mesh%ncell
      associate (cnode => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        n = num_cell_faces(cnode)
        this%xminv(j+1) = this%xminv(j) + (n*(n+1))/2
      end associate
    end do
    allocate(this%minv(this%xminv(this%mesh%ncell+1)-1))

    select case (this%method)
    case (MFD_STANDARD)
      call init_standard_minv(this)
    case (MFD_LEGACY)
      call init_legacy_minv(this)
    end select

  end subroutine init_minv

  subroutine init_standard_minv(this)

    type(mfd_disc), intent(inout) :: this

    integer :: j
    type(mfd_cell) :: cell

    do j = 1, this%mesh%ncell
      associate (cnode => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1), &
                 minv => this%minv(this%xminv(j):this%xminv(j+1)-1))
        call cell%init (this%mesh%x(:,cnode))
        call cell%compute_flux_matrix_inv (1.0_r8, minv)
      end associate
    end do

  end subroutine init_standard_minv

  integer function max_cell_faces(mesh)
    type(unstr_mesh), intent(in) :: mesh
    if (mesh%ncell == 0) then
      max_cell_faces = 0
    else
      max_cell_faces = maxval(mesh%xcface(2:) - mesh%xcface(:mesh%ncell))
    end if
  end function

  !! Compute a cell-wise average gradient of a scalar field from face values,
  !! primarily for output and diagnostics.  This ignores any available cell
  !! averages and is not intended as a discretization primitive.  Gradients are
  !! computed for the leading cell range selected by size(GRAD,2), which may be
  !! any value up to the full local cell count.

  subroutine compute_cell_grad1(this, uface, grad)
    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    real(r8), intent(out) :: grad(:,:)
    select case (this%method)
    case (MFD_STANDARD)
      call compute_cell_grad_standard(this, uface, grad)
    case (MFD_LEGACY)
      call compute_cell_grad_legacy(this, uface, grad)
    end select
  end subroutine

  !! Compute the average cell gradient from face values using the Gauss-Green
  !! formula over the caller-selected leading cell range.

  subroutine compute_cell_grad_standard(this, uface, grad)

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    real(r8), intent(out) :: grad(:,:)

    integer :: j
    type(mfd_cell) :: cell

    INSIST(size(grad,1) == 3)
    INSIST(size(uface) == this%mesh%nface)
    INSIST(size(grad,2) <= this%mesh%ncell)

    do j = 1, size(grad,2)
      associate (cnode => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1), &
                 cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        call cell%init (this%mesh%x(:,cnode))
        grad(:,j) = matmul(cell%face_normals(:,:cell%nfaces), uface(cface)) / cell%volume
      end associate
    end do

  end subroutine compute_cell_grad_standard

  !! Masked overload of COMPUTE_CELL_GRAD.  Over the caller-selected leading
  !! cell range, gradients are computed only where MASK is true; inactive
  !! entries of GRAD are set to zero.

  subroutine compute_cell_grad2(this, uface, mask, grad)
    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    logical,  intent(in) :: mask(:)
    real(r8), intent(out) :: grad(:,:)
    select case (this%method)
    case (MFD_STANDARD)
      call compute_cell_grad_masked_standard(this, uface, mask, grad)
    case (MFD_LEGACY)
      call compute_cell_grad_masked_legacy(this, uface, mask, grad)
    end select
  end subroutine

  subroutine compute_cell_grad_masked_standard(this, uface, mask, grad)

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    logical,  intent(in) :: mask(:)
    real(r8), intent(out) :: grad(:,:)

    integer :: j
    type(mfd_cell) :: cell

    INSIST(size(grad,1) == 3)
    INSIST(size(uface) == this%mesh%nface)
    INSIST(size(grad,2) <= this%mesh%ncell)
    INSIST(size(mask) == size(grad,2))

    do j = 1, size(grad,2)
      if (mask(j)) then
        associate (cnode => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1), &
                   cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
          call cell%init (this%mesh%x(:,cnode))
          grad(:,j) = matmul(cell%face_normals(:,:cell%nfaces), uface(cface)) / cell%volume
        end associate
      else
        grad(:,j) = 0.0_r8
      end if
    end do

  end subroutine compute_cell_grad_masked_standard

  !!!! MFD_CELL type bound procedures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mfd_cell_init(this, vertices)
    use cell_geometry, only: cell_volume, get_cell_face_centers, get_cell_face_normals, &
      vector_length
    class(mfd_cell), intent(inout) :: this
    real(r8), intent(in) :: vertices(:,:)
    integer :: j, nface
    ASSERT(size(vertices,dim=1) == 3)
    call get_cell_face_centers(vertices, this%face_centers, this%nfaces)
    call get_cell_face_normals(vertices, this%face_normals, nface)
    INSIST(this%nfaces > 0)
    ASSERT(nface == this%nfaces)
    this%cell_center  = sum(vertices,dim=2) / size(vertices,dim=2)
    do j = 1, this%nfaces
      this%face_area(j) = vector_length(this%face_normals(:,j))
    end do
    this%volume = cell_volume(vertices)
  end subroutine

  subroutine mfd_cell_dump(this)
    class(mfd_cell), intent(in) :: this
    integer :: j
    write(*,'(a,es10.2,a,3es10.2)') 'VOLUME=', this%volume, ', CENTER=', this%cell_center
    write(*,'((a,i0,a,es10.2,a,3es10.2))') &
        ('FACE(', j, '): AREA=', this%face_area(j), ', CENTER=', this%face_centers(:,j), j=1,this%nfaces)
    write(*,'((a,i0,a,3es10.2))') &
        ('FACE(', j, '): NORMAL=', this%face_normals(:,j), j=1,this%nfaces)
  end subroutine

  !! Compute the inverse MFD mass matrix for a polyhedral cell using the
  !! standard parameter choice; see references listed above.  With N the matrix
  !! of face normals and Q an orthonormal basis for range(R), the inverse is
  !! W = (coef/volume) * (N*N' + I - Q*Q').

  subroutine mfd_cell_compute_flux_matrix_inv(this, coef, matrix)

    class(mfd_cell), intent(in) :: this
    real(r8), intent(in) :: coef
    real(r8), intent(out) :: matrix(:)

    integer :: i, j, k, n, loc
    real(r8) :: stab_val, rv(3)
    real(r8) :: Q(3,MFD_CELL_NFACE_MAX), R(3,MFD_CELL_NFACE_MAX)

    matrix = 0.0_r8

    n = this%nfaces
    ASSERT(n <= MFD_CELL_NFACE_MAX)

    do i = 1, n
      R(:, i) = (this%face_centers(:,i) - this%cell_center)*this%face_area(i)
    end do

    !! Modified Gram-Schmidt orthogonalization.

    do k = 1, 3
      rv(k) = 0
      do i = 1, n
        rv(k) = rv(k) + R(k,i)*R(k,i)
      end do
      rv(k) = sqrt(rv(k))
      Q(k,:n) = R(k,:n)/rv(k)
      do j = k+1, 3
        rv(j) = 0
        do i = 1, n
          rv(j) = rv(j) + Q(k,i)*R(j,i)
        end do
        R(j,:n) = R(j,:n) - Q(k,:n)*rv(j)
      end do
    end do

    do i = 1, n
      do j = i, n
        loc = i + j*(j-1)/2
        matrix(loc) = 0
        do k = 1, 3
          matrix(loc) = matrix(loc) + this%face_normals(k,i)*this%face_normals(k,j)
        end do
        stab_val = 0
        do k = 1, 3
          stab_val = stab_val - Q(k,i)*Q(k,j)
        end do
        if (i == j) stab_val = stab_val + 1
        stab_val = stab_val * this%face_area(i) * this%face_area(j)
        matrix(loc) = (matrix(loc) + stab_val)*(coef/this%volume)
      end do
    end do

  end subroutine mfd_cell_compute_flux_matrix_inv

  !! Legacy MFD procedures
  !!
  !! These procedures implement the historical shape-specific MFD construction
  !! used only when mfd_disc%init is called with METHOD=MFD_LEGACY.

  subroutine init_legacy_minv(this)

    type(mfd_disc), intent(inout) :: this

    integer :: j
    type(mfd_tet) :: tet
    type(mfd_hex) :: hex
    type(mfd_wedge) :: wedge
    type(mfd_cell) :: cell

    do j = 1, this%mesh%ncell
      associate (cnode => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1), &
                 minv => this%minv(this%xminv(j):this%xminv(j+1)-1))
        select case (size(cnode))
        case (4)  ! tet
          call tet%init (this%mesh%x(:,cnode))
          call tet%compute_flux_matrix (1.0_r8, minv, invert=.true.)
        case (5)  ! pyramid
          call cell%init (this%mesh%x(:,cnode))
          call cell%compute_flux_matrix_inv (1.0_r8, minv)
        case (6)  ! wedge
          call wedge%init (this%mesh%x(:,cnode))
          call wedge%compute_flux_matrix (1.0_r8, minv, invert=.true.)
        case (8)  ! hex
          call hex%init (this%mesh%x(:,cnode))
          call hex%compute_flux_matrix (1.0_r8, minv, invert=.true.)
        case default
          INSIST(.false.)
        end select
      end associate
    end do

  end subroutine init_legacy_minv

  !! Compute the average cell gradient from face values using the Gauss-Green
  !! formula with legacy shape-specific geometry helpers over the
  !! caller-selected leading cell range.

  subroutine compute_cell_grad_legacy(this, uface, grad)

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    real(r8), intent(out) :: grad(:,:)

    integer :: j
    type(mfd_tet) :: tet
    type(mfd_wedge) :: wedge
    type(mfd_hex) :: hex
    type(mfd_cell) :: cell

    INSIST(size(grad,1) == 3)
    INSIST(size(uface) == this%mesh%nface)
    INSIST(size(grad,2) <= this%mesh%ncell)

    do j = 1, size(grad,2)
      associate (cnode => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1), &
                 cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        select case (size(cnode))
        case (4)
          call tet%init (this%mesh%x(:,cnode))
          grad(:,j) = matmul(tet%face_normals, uface(cface)) / tet%volume
        case (5)  ! pyramid
          call cell%init (this%mesh%x(:,cnode))
          grad(:,j) = matmul(cell%face_normals(:,:cell%nfaces), uface(cface)) / cell%volume
        case (6)
          call wedge%init (this%mesh%x(:,cnode))
          grad(:,j) = matmul(wedge%face_normals, uface(cface)) / wedge%volume
        case (8)
          call hex%init (this%mesh%x(:,cnode))
          grad(:,j) = matmul(hex%face_normals, uface(cface)) / hex%volume
        case default
          INSIST(.false.)
        end select
      end associate
    end do

  end subroutine compute_cell_grad_legacy

  subroutine compute_cell_grad_masked_legacy(this, uface, mask, grad)

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    logical,  intent(in) :: mask(:)
    real(r8), intent(out) :: grad(:,:)

    integer :: j
    type(mfd_tet) :: tet
    type(mfd_wedge) :: wedge
    type(mfd_hex) :: hex
    type(mfd_cell) :: cell

    INSIST(size(grad,1) == 3)
    INSIST(size(uface) == this%mesh%nface)
    INSIST(size(grad,2) <= this%mesh%ncell)
    INSIST(size(mask) == size(grad,2))

    do j = 1, size(grad,2)
      if (mask(j)) then
        associate (cnode => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1), &
                   cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
          select case (size(cnode))
          case (4)
            call tet%init (this%mesh%x(:,cnode))
            grad(:,j) = matmul(tet%face_normals, uface(cface)) / tet%volume
          case (5)  ! pyramid
            call cell%init (this%mesh%x(:,cnode))
            grad(:,j) = matmul(cell%face_normals(:,:cell%nfaces), uface(cface)) / cell%volume
          case (6)
            call wedge%init (this%mesh%x(:,cnode))
            grad(:,j) = matmul(wedge%face_normals, uface(cface)) / wedge%volume
          case (8)
            call hex%init (this%mesh%x(:,cnode))
            grad(:,j) = matmul(hex%face_normals, uface(cface)) / hex%volume
          case default
            INSIST(.false.)
          end select
        end associate
      else
        grad(:,j) = 0.0_r8
      end if
    end do

  end subroutine compute_cell_grad_masked_legacy

  subroutine mfd_tet_init (this, vertices)
    use cell_geometry, only: tet_volume, tet_face_normals
    class(mfd_tet), intent(out) :: this
    real(r8), intent(in) :: vertices(:,:)
    ASSERT(size(vertices,1) == 3 .and. size(vertices,2) == 4)
    this%volume = tet_volume(vertices)
    this%face_normals = tet_face_normals(vertices)
  end subroutine

  subroutine mfd_tet_compute_flux_matrix(this, coef, matrix, invert)

    use cell_topology, only: TET4_VERT_FACE
    use upper_packed_matrix_procs, only: upm_invert

    class(mfd_tet), intent(in) :: this
    real(r8), intent(in) :: coef
    real(r8), intent(out) :: matrix(:)
    logical, intent(in), optional :: invert

    integer :: c, i, j, ii, jj, loc
    real(r8) :: s, Nc(3,3), Mc(3,3)

    ASSERT(size(matrix) == 10)

    matrix = 0.0_r8
    do c = 1, 4
      Nc = this%face_normals(:,TET4_VERT_FACE(:,c))
      call invert_sym_3x3 (matmul(transpose(Nc),Nc), Mc)
      !! Scatter the corner matrix into the full cell flux matrix.
      !! It is essential that TET4_VERT_FACE(:,c) is an increasing sequence of indices.
      s = (0.25/coef)*this%volume
      do j = 1, 3
        jj = TET4_VERT_FACE(j,c)
        do i = 1, j
          ii = TET4_VERT_FACE(i,c)
          loc = ii + jj*(jj - 1)/2
          matrix(loc) = matrix(loc) + s*Mc(i,j)
        end do
      end do
    end do

    if (present(invert)) then
      if (invert) call upm_invert (matrix)
    end if

  end subroutine mfd_tet_compute_flux_matrix

  subroutine mfd_wedge_init(this, vertices)
    use cell_geometry, only: cell_volume, wedge_face_normals, tet_volume
    class(mfd_wedge), intent(out) :: this
    real(r8), intent(in) :: vertices(:,:)
    ASSERT(size(vertices,1) == 3 .and. size(vertices,2) == 6)
    this%volume = cell_volume(vertices)
    this%face_normals = wedge_face_normals(vertices)
    this%corner_volumes(1) = tet_volume(vertices(:,[1,2,3,4]))
    this%corner_volumes(2) = tet_volume(vertices(:,[2,3,1,5]))
    this%corner_volumes(3) = tet_volume(vertices(:,[3,1,2,6]))
    this%corner_volumes(4) = tet_volume(vertices(:,[4,6,5,1]))
    this%corner_volumes(5) = tet_volume(vertices(:,[5,4,6,2]))
    this%corner_volumes(6) = tet_volume(vertices(:,[6,5,4,3]))
  end subroutine

  subroutine mfd_wedge_compute_flux_matrix(this, coef, matrix, invert)

    use cell_topology, only: WED6_VERT_FACE
    use upper_packed_matrix_procs, only: upm_invert

    class(mfd_wedge), intent(in) :: this
    real(r8), intent(in) :: coef
    real(r8), intent(out) :: matrix(:)
    logical, intent(in), optional :: invert

    integer :: c, i, j, ii, jj, loc
    real(r8) :: s, cwt(6), Nc(3,3), Mc(3,3)

    ASSERT(size(matrix) == 15)

    matrix = 0.0_r8
    cwt = this%corner_volumes / sum(this%corner_volumes)
    do c = 1, 6
      Nc = this%face_normals(:,WED6_VERT_FACE(:,c))
      call invert_sym_3x3(matmul(transpose(Nc),Nc), Mc)
      !! Scatter the corner matrix into the full cell flux matrix.
      !! It is essential that WED6_VERT_FACE(:,c) is an increasing sequence of indices.
      s = this%volume * cwt(c) / coef
      do j = 1, 3
        jj = WED6_VERT_FACE(j,c)
        do i = 1, j
          ii = WED6_VERT_FACE(i,c)
          loc = ii + jj*(jj - 1)/2
          matrix(loc) = matrix(loc) + s*Mc(i,j)
        end do
      end do
    end do

    if (present(invert)) then
      if (invert) call upm_invert (matrix)
    end if

  end subroutine mfd_wedge_compute_flux_matrix

  subroutine mfd_hex_init (this, vertices)
    use cell_geometry, only: eval_hex_volumes, hex_face_normals
    class(mfd_hex), intent(out) :: this
    real(r8), intent(in) :: vertices(:,:)
    ASSERT(size(vertices,1) == 3 .and. size(vertices,2) == 8)
    call eval_hex_volumes (vertices, this%volume, this%corner_volumes)
    this%face_normals = hex_face_normals(vertices)
  end subroutine

  subroutine mfd_hex_compute_flux_matrix(this, coef, matrix, invert)

    use cell_topology, only: HEX8_VERT_FACE
    use upper_packed_matrix_procs, only: upm_invert

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
      if (invert) call upm_invert (matrix)
    end if

  end subroutine mfd_hex_compute_flux_matrix

 !! Direct inversion of a 3x3 symmetrix matrix using the formula that the
 !! inverse equals the transponse of the matrix of cofactors divided by
 !! the determinant.

  subroutine invert_sym_3x3(A, Ainv)
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
  end subroutine

end module mfd_disc_type
