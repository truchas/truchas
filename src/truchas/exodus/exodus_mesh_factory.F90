!!
!! EXODUS_MESH_FACTORY
!!
!! Some procedures for generating simple 2D and 3D rectilinear Exodus meshes
!! of rectangle and brick domains.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The generic subroutine CREATE_RECT_EXODUS_MESH provides several interfaces
!! for initializing the contents of an EXODUS_MESH class object MESH. In all
!! cases the mesh consists of a single element block with ID 1, and defines a
!! side set for each side of the rectangle or brick: IDs 1, 2, 3, 4 for the
!! x=x_min, x=x_max, y=y_min, y=y_max sides, and in 3D also IDs 5 and 6 for
!! the z=z_min and z=z_max sides.
!!
!!  CALL CREATE_RECT_EXODUS_MESH(MESH, X, Y, Z [,EPS])
!!    REAL(INT64), INTENT(IN) :: X(:), Y(:), Z(:), EPS
!!
!!  This creates a 3D hexahedral mesh whose node coordinates are the tensor
!!  product of the discretization arrays X, Y, and Z for each of the coordinate
!!  axes. The first and last elements of each array define the extent of the
!!  domain in that coordinate. If EPS is specified, the node coordinates will
!!  be randomly perturbed by an amount whose magnitude is at most EPS times
!!  the local cell size. Nodes on the boundary will not be perturbed in
!!  directions normal to the sides, so that the domain itself is preserved.
!!
!!  CALL CREATE_RECT_EXODUS_MESH(MESH, X, Y, [,EPS])
!!
!!  Same as the preceding case except that it creates a 2D quadrilateral mesh.
!!
!!  CALL CREATE_RECT_EXODUS_MESH(MESH, XMIN, XMAX, NCELL [,RATIO] [,EPS])
!!    REAL(INT64), INTENT(IN) :: XMIN(:), XMAX(:), RATIO(:), EPS
!!    INTEGER, INTENT(IN) :: NCELL(:)
!!
!!  This creates a rectilinear mesh of the rectangle or brick domain defined
!!  by the coordinates XMIN and XMAX of its diametrically opposed extreme
!!  corners. The array arguments must all have the same size, either 2 or 3,
!!  corresponding to the dimension of the mesh. The NCELL array specifies
!!  the number of intervals into which the domain should be subdivided in each
!!  of the coordinate directions. By default the intervals are equally sized,
!!  but if RATIO is present, it defines a size biasing to use for each
!!  coordinate: the ratio of the lengths of successive intervals will equal
!!  the value of RATIO for that coordinate. The optional EPS argument is the
!!  same as previously described.
!!
!!  CALL CREATE_RECT_EXODUS_MESH(MESH, PARAMS, STAT, ERRMSG)
!!    TYPE(PARAMETER_LIST), INTENT(INOUT) :: PARAMS
!!    INTEGER, INTENT(OUT) :: STAT
!!    CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: ERRMSG
!!
!!  This creates a rectilinear mesh described by the parameter list PARAMS.
!!  If an error is encountered, STAT returns a non-zero value and ERRMSG an
!!  explanatory message. The parameter list has the following format:
!!
!!    {
!!      "dimen": 2 or 3 (default 3)
!!      "x-axis": axis-sublist
!!      "y-axis": axis-sublist
!!      "z-axis": axis-sublist (if dimen == 3)
!!      "noise-factor": real (positive, default 0)
!!    }
!!
!!    where axis-sublist describes the grid for a coordinate axis and has the
!!    following format:
!!
!!    {
!!      "coarse-grid": real-vector (strictly increasing, size n+1)
!!      "intervals":   integer-vector (positive, size n)
!!      "ratio":       real-vector (positive, size n, default 1.0)
!!    }
!!
!!  Each interval of the coarse coordinate grid given by the coarse-grid
!!  parameter is subdivided into the number of subintervals prescribed by
!!  the corresponding element of the intervals parameter. By default the
!!  subintervals are equal size, but if the ratio parameter is specified
!!  its value gives the ratio of the lengths of consecutive subintervals.
!!  If the noise-factor parameter is specified with a positive value, the
!!  node coordinates of the tensor product grid will be perturbed by a
!!  random amount whose magnitude will not exceed this value times the
!!  local cell size at each node.
!!

#include "f90_assert.fpp"

module exodus_mesh_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use exodus_mesh_type
  implicit none
  private

  public :: exodus_mesh, create_rect_exodus_mesh

  interface create_rect_exodus_mesh
    procedure create_2d_rect_exodus_mesh
    procedure create_3d_rect_exodus_mesh
    procedure create_rect_exodus_mesh_block
    procedure create_rect_exodus_mesh_params
  end interface

  type :: coord_grid
    private
    real(r8), allocatable :: grid(:)
    integer,  allocatable :: nint(:)
    real(r8), allocatable :: ratio(:)
  contains
    procedure :: init => coord_grid_init
    procedure :: get_grid
  end type

contains

  !! Create a 2D or 3D rectilinear mesh according to parameter list input.

  subroutine create_rect_exodus_mesh_params(mesh, params, stat, errmsg)

    use parameter_list_type

    class(exodus_mesh), intent(out) :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: dim
    real(r8) :: eps
    real(r8), allocatable :: x(:), y(:), z(:)

    call params%get('dimen', dim, default=3, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (dim < 2 .or. dim > 3) then
      stat = 1
      errmsg = 'dimen must be either 2 or 3'
      return
    end if

    call params%get('noise-factor', eps, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (eps < 0.0_r8) then
      stat = 1
      errmsg = 'noise-factor must be >= 0'
      return
    else if (eps > 0.3_r8) then
      stat = 1
      errmsg = 'noise-factor must be <= 0.3'
      return
    end if

    call get_axis_grid(params, 'x-axis', x, stat, errmsg)
    if (stat /= 0) return

    call get_axis_grid(params, 'y-axis', y, stat, errmsg)
    if (stat /= 0) return

    select case (dim)
    case (2)
      call create_2d_rect_exodus_mesh(mesh, x, y, eps)
    case (3)
      call get_axis_grid(params, 'z-axis', z, stat, errmsg)
      if (stat /= 0) return
      call create_3d_rect_exodus_mesh(mesh, x, y, z, eps)
    case default
      ASSERT(.false.)
    end select

  contains

    subroutine get_axis_grid(params, key, x, stat, errmsg)
      type(parameter_list), intent(inout) :: params
      character(*), intent(in) :: key
      real(r8), allocatable, intent(out) :: x(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      type(coord_grid) :: grid
      type(parameter_list), pointer :: plist
      if (params%is_sublist(key)) then
        plist => params%sublist(key)
        call grid%init(plist, 'coarse-grid', 'intervals', 'ratio', stat, errmsg)
        if (stat /= 0) then
          errmsg = key // ': ' // errmsg
          return
        end if
        call grid%get_grid(x)
      else
        stat = 1
        errmsg = 'missing ' // key // ' sublist'
      end if
    end subroutine

  end subroutine create_rect_exodus_mesh_params

  !! Generate a simple, 2D or 3D, single-block rectilinear mesh, given the
  !! coordinates of the extreme corners of the domain, the number of intervals
  !! in each direction, and an optional node perturbation factor.

  subroutine create_rect_exodus_mesh_block(mesh, xmin, xmax, nx, ratio, eps)

    class(exodus_mesh), intent(out) :: mesh
    real(r8), intent(in) :: xmin(:), xmax(:)
    integer,  intent(in) :: nx(:)
    real(r8), intent(in), optional :: ratio(:)
    real(r8), intent(in), optional :: eps

    integer :: dim
    real(r8), allocatable :: x(:), y(:), z(:), r(:)

    dim = size(xmin)

    INSIST(dim >=2 .and. dim <=3)
    INSIST(size(xmax) == dim)
    INSIST(all(xmax > xmin))
    INSIST(size(nx) == dim)
    INSIST(all(nx > 0))

    allocate(r(dim))
    r = 1.0_r8
    if (present(ratio)) then
      INSIST(size(ratio) == 2)
      INSIST(all(ratio > 0))
      r = ratio
    end if

    allocate(x(0:nx(1)))
    x(0) = xmin(1); x(nx(1)) = xmax(1)
    call discretize_interval(r(1), x)

    allocate(y(0:nx(2)))
    y(0) = xmin(2); y(nx(2)) = xmax(2)
    call discretize_interval(r(2), y)

    select case (dim)
    case (2)
      call create_2d_rect_exodus_mesh(mesh, x, y, eps)
    case (3)
      allocate(z(0:nx(3)))
      z(0) = xmin(3); z(nx(3)) = xmax(3)
      call discretize_interval(r(3), z)
      call create_3d_rect_exodus_mesh(mesh, x, y, z, eps)
    end select

  end subroutine create_rect_exodus_mesh_block

  !! Generate a 2D rectilinear quad mesh given grids in each of the coordinate
  !! directions, with option to randomly perturb nodes.

  subroutine create_2d_rect_exodus_mesh(mesh, x, y, eps)

    class(exodus_mesh), intent(out) :: mesh
    real(r8), intent(in) :: x(0:), y(0:)
    real(r8), intent(in), optional :: eps

    integer :: i, j, n, nzone(2)

    nzone(1) = size(x) - 1
    nzone(2) = size(y) - 1

    !! Mesh is a tensor mesh with these node coordinates
    mesh%num_dim = 2
    mesh%num_elem = product(nzone)
    mesh%num_node = product(nzone+1)
    allocate(mesh%coord(mesh%num_dim,mesh%num_node))
    do concurrent (i = 0:nzone(1), j = 0:nzone(2))
      n = node_index(i,j)
      mesh%coord(1,n) = x(i)
      mesh%coord(2,n) = y(j)
    end do

    !! Perturb the node coordinates if requested
    if (present(eps)) call randomize_coord(eps)

    !! All cells go into a single element block (ID=1)
    mesh%num_eblk = 1
    allocate(mesh%eblk(mesh%num_eblk))
    mesh%eblk(1)%id = 1
    mesh%eblk(1)%num_elem = mesh%num_elem
    mesh%eblk(1)%num_nodes_per_elem = 4
    mesh%eblk(1)%elem_type = 'QUAD'

    !! Generate the regular quad connectivity
    allocate(mesh%eblk(1)%connect(4,mesh%eblk(1)%num_elem))
    do concurrent (i=1:nzone(1), j = 1:nzone(2))
      n = cell_index(i,j)
      mesh%eblk(1)%connect(1,n) = node_index(i-1,j-1)
      mesh%eblk(1)%connect(2,n) = node_index(i,j-1)
      mesh%eblk(1)%connect(3,n) = node_index(i,j)
      mesh%eblk(1)%connect(4,n) = node_index(i-1,j)
    end do
    if (present(eps)) call randomize_coord(eps)

    !! Generate a side set for each of the sides of the rectangle
    mesh%num_sset = 4
    allocate(mesh%sset(mesh%num_sset))

    !! Side sets 1 and 2 (x = xmin, xmax)
    mesh%sset(1)%id = 1
    mesh%sset(2)%id = 2
    n = nzone(2)
    mesh%sset(1)%num_side = n
    mesh%sset(2)%num_side = n
    allocate(mesh%sset(1)%elem(n), mesh%sset(1)%face(n))
    allocate(mesh%sset(2)%elem(n), mesh%sset(2)%face(n))
    do concurrent (j = 1:nzone(2))
      mesh%sset(1)%elem(j) = cell_index(1,j)
      mesh%sset(2)%elem(j) = cell_index(nzone(1),j)
    end do
    mesh%sset(1)%face = 4
    mesh%sset(2)%face = 2

    !! Side sets 3 and 4 (y = ymin, ymax)
    mesh%sset(3)%id = 3
    mesh%sset(4)%id = 4
    n = nzone(1)
    mesh%sset(3)%num_side = n
    mesh%sset(4)%num_side = n
    allocate(mesh%sset(3)%elem(n), mesh%sset(3)%face(n))
    allocate(mesh%sset(4)%elem(n), mesh%sset(4)%face(n))
    do concurrent (i = 1:nzone(1))
      mesh%sset(3)%elem(i) = cell_index(i,1)
      mesh%sset(4)%elem(i) = cell_index(i,nzone(2))
    end do
    mesh%sset(3)%face = 1
    mesh%sset(4)%face = 3

    !! No node sets
    mesh%num_nset = 0
    allocate(mesh%nset(0))

  contains

    !! Usual Fortran column-major ordering map
    pure integer function cell_index(i, j)
      integer, intent(in) :: i, j
      cell_index = i + (j-1)*nzone(1)
    end function

    !! Usual Fortran column-major ordering map
    pure integer function node_index(i, j)
      integer, intent(in) :: i, j
      node_index = 1 + i + j*(nzone(1)+1)
    end function

    !! Randomize the coordinates of the tensor product mesh
    subroutine randomize_coord(eps)

      real(r8), intent(in) :: eps

      integer :: i, j, n
      real(r8), allocatable :: dx(:), dy(:)
      logical :: mask(2)
      real(r8) :: d(2)

      if (eps == 0) return

      !! Get the local interval size in each coordinate direction. The local
      !! cell size at a node will be the minimum of the local interval sizes
      !! in each direction. Coordinates will be randomly perturbed by a
      !! random amount up to the fraction EPS of this size. Nodes on the
      !! boundary are not perturbed in a direction normal to the boundary.

      allocate(dx, mold=x)
      allocate(dy, mold=y)
      call dmin(x, dx)
      call dmin(y, dy)

      do j = 0, nzone(2)
        mask(2) = (j > 0 .and. j < nzone(2))
        do i = 0, nzone(1)
          mask(1) = (i > 0 .and. i < nzone(1))
          call random_number(d) ! in [0,1)
          d = eps*min(dx(i),dy(j))*(2*d - 1)
          n = node_index(i,j)
          mesh%coord(:,n) = mesh%coord(:,n) + merge(d, 0.0_r8, mask)
        end do
      end do

    end subroutine randomize_coord

    !! Compute a local interval size for each point in a discretization X(1:N).
    !! This is the min of the lengths of the intervals on either side of point.

    subroutine dmin(x, dx)
      real(r8), intent(in)  :: x(:)
      real(r8), intent(out) :: dx(:)
      integer :: j, n
      real(r8) :: d1, d2
      n = ubound(x,1)
      d2 = x(2) - x(1)
      dx(1) = d2
      do j = 2, n-1
        d1 = d2
        d2 = x(j+1) - x(j)
        dx(j) = min(d1, d2)
      end do
      dx(n) = d2
    end subroutine dmin

  end subroutine

  !! Generate a 3D rectilinear quad mesh given grids in each of the coordinate
  !! directions, with option to randomly perturb nodes.

  subroutine create_3d_rect_exodus_mesh(mesh, x, y, z, eps)

    class(exodus_mesh), intent(out) :: mesh
    real(r8), intent(in) :: x(0:), y(0:), z(0:)
    real(r8), intent(in), optional :: eps

    integer :: i, j, k, n, nzone(3)

    nzone(1) = size(x) - 1
    nzone(2) = size(y) - 1
    nzone(3) = size(z) - 1

    !! Mesh is a tensor mesh with these node coordinates
    mesh%num_dim = 3
    mesh%num_elem = product(nzone)
    mesh%num_node = product(nzone+1)
    allocate(mesh%coord(mesh%num_dim,mesh%num_node))
    do concurrent (i = 0:nzone(1), j = 0:nzone(2), k = 0:nzone(3))
      n = node_index(i,j,k)
      mesh%coord(1,n) = x(i)
      mesh%coord(2,n) = y(j)
      mesh%coord(3,n) = z(k)
    end do

    !! Perturb the node coordinates if requested
    if (present(eps)) call randomize_coord(eps)

    !! All cells go into a single element block (ID=1)
    mesh%num_eblk = 1
    allocate(mesh%eblk(mesh%num_eblk))
    mesh%eblk(1)%id = 1
    mesh%eblk(1)%num_elem = mesh%num_elem
    mesh%eblk(1)%num_nodes_per_elem = 8
    mesh%eblk(1)%elem_type = 'HEX8'

    !! Generate the regular hex connectivity
    allocate(mesh%eblk(1)%connect(8,mesh%eblk(1)%num_elem))
    do concurrent (i=1:nzone(1), j = 1:nzone(2), k = 1:nzone(3))
      n = cell_index(i,j,k)
      mesh%eblk(1)%connect(1,n) = node_index(i-1,j-1,k-1)
      mesh%eblk(1)%connect(2,n) = node_index(i,j-1,k-1)
      mesh%eblk(1)%connect(3,n) = node_index(i,j,k-1)
      mesh%eblk(1)%connect(4,n) = node_index(i-1,j,k-1)
      mesh%eblk(1)%connect(5,n) = node_index(i-1,j-1,k)
      mesh%eblk(1)%connect(6,n) = node_index(i,j-1,k)
      mesh%eblk(1)%connect(7,n) = node_index(i,j,k)
      mesh%eblk(1)%connect(8,n) = node_index(i-1,j,k)
    end do
    if (present(eps)) call randomize_coord(eps)

    !! Generate a side set for each of the sides of the brick
    mesh%num_sset = 6
    allocate(mesh%sset(mesh%num_sset))

    !! Side sets 1 and 2 (x = xmin, xmax)
    mesh%sset(1)%id = 1
    mesh%sset(2)%id = 2
    n = nzone(2)*nzone(3)
    mesh%sset(1)%num_side = n
    mesh%sset(2)%num_side = n
    allocate(mesh%sset(1)%elem(n), mesh%sset(1)%face(n))
    allocate(mesh%sset(2)%elem(n), mesh%sset(2)%face(n))
    do concurrent (j = 1:nzone(2), k = 1:nzone(3))
      n = j + (k-1)*nzone(2)
      mesh%sset(1)%elem(n) = cell_index(1,j,k)
      mesh%sset(2)%elem(n) = cell_index(nzone(1),j,k)
    end do
    mesh%sset(1)%face = 4
    mesh%sset(2)%face = 2

    !! Side sets 3 and 4 (y = ymin, ymax)
    mesh%sset(3)%id = 3
    mesh%sset(4)%id = 4
    n = nzone(1)*nzone(3)
    mesh%sset(3)%num_side = n
    mesh%sset(4)%num_side = n
    allocate(mesh%sset(3)%elem(n), mesh%sset(3)%face(n))
    allocate(mesh%sset(4)%elem(n), mesh%sset(4)%face(n))
    do concurrent (i = 1:nzone(1), k = 1:nzone(3))
      n = i + (k-1)*nzone(1)
      mesh%sset(3)%elem(n) = cell_index(i,1,k)
      mesh%sset(4)%elem(n) = cell_index(i,nzone(2),k)
    end do
    mesh%sset(3)%face = 1
    mesh%sset(4)%face = 3

    !! Side sets 5 and 6 (z = zmin, zmax)
    mesh%sset(5)%id = 5
    mesh%sset(6)%id = 6
    n = nzone(1)*nzone(2)
    mesh%sset(5)%num_side = n
    mesh%sset(6)%num_side = n
    allocate(mesh%sset(5)%elem(n), mesh%sset(5)%face(n))
    allocate(mesh%sset(6)%elem(n), mesh%sset(6)%face(n))
    do concurrent (i = 1:nzone(1), j = 1:nzone(2))
      n = i + (j-1)*nzone(1)
      mesh%sset(5)%elem(n) = cell_index(i,j,1)
      mesh%sset(6)%elem(n) = cell_index(i,j,nzone(3))
    end do
    mesh%sset(5)%face = 5
    mesh%sset(6)%face = 6

    !! No node sets
    mesh%num_nset = 0
    allocate(mesh%nset(0))

  contains

    !! Usual Fortran column-major ordering map
    pure integer function cell_index(i, j, k)
      integer, intent(in) :: i, j, k
      cell_index = i + ((j-1) + (k-1)*nzone(2))*nzone(1)
    end function

    !! Usual Fortran column-major ordering map
    pure integer function node_index(i, j, k)
      integer, intent(in) :: i, j, k
      node_index = 1 + i + (j + k*(nzone(2)+1))*(nzone(1)+1)
    end function

    !! Randomize the coordinates of the tensor product mesh
    subroutine randomize_coord(eps)

      real(r8), intent(in) :: eps

      integer :: i, j, k, n
      real(r8), allocatable :: dx(:), dy(:), dz(:)
      logical :: mask(3)
      real(r8) :: d(3)

      if (eps == 0) return

      !! Get the local interval size in each coordinate direction. The local
      !! cell size at a node will be the minimum of the local interval sizes
      !! in each direction. Coordinates will be randomly perturbed by a
      !! random amount up to the fraction EPS of this size. Nodes on the
      !! boundary are not perturbed in a direction normal to the boundary.

      allocate(dx, mold=x)
      allocate(dy, mold=y)
      allocate(dz, mold=z)
      call dmin(x, dx)
      call dmin(y, dy)
      call dmin(z, dz)

      do k = 0, nzone(3)
        mask(3) = (k > 0 .and. k < nzone(3))
        do j = 0, nzone(2)
          mask(2) = (j > 0 .and. j < nzone(2))
          do i = 0, nzone(1)
            mask(1) = (i > 0 .and. i < nzone(1))
            call random_number(d) ! in [0,1)
            d = eps*min(dx(i),dy(j),dz(k))*(2*d - 1)
            n = node_index(i,j,k)
            mesh%coord(:,n) = mesh%coord(:,n) + merge(d, 0.0_r8, mask)
          end do
        end do
      end do

    end subroutine randomize_coord

    !! Compute a local interval size for each point in a discretization X(1:N).
    !! This is the min of the lengths of the intervals on either side of point.

    subroutine dmin(x, dx)
      real(r8), intent(in)  :: x(:)
      real(r8), intent(out) :: dx(:)
      integer :: j, n
      real(r8) :: d1, d2
      n = ubound(x,1)
      d2 = x(2) - x(1)
      dx(1) = d2
      do j = 2, n-1
        d1 = d2
        d2 = x(j+1) - x(j)
        dx(j) = min(d1, d2)
      end do
      dx(n) = d2
    end subroutine dmin

  end subroutine

!!!! INTERVAL_PARTITION TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coord_grid_init(this, params, grid_key, nint_key, ratio_key, stat, errmsg)

    use parameter_list_type
    use string_utilities, only: i_to_c

    class(coord_grid), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    character(*), intent(in) :: grid_key, nint_key, ratio_key
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: rarray(:)
    integer,  allocatable :: iarray(:)
    integer :: n

    !! The coarse grid points (required)
    call params%get(grid_key, rarray, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    n = size(rarray)
    if (n < 2) then
      stat = 1
      errmsg = grid_key // ' requires at least two values'
      return
    else if (any(rarray(2:n) <= rarray(1:n-1))) then
      stat = 1
      errmsg = grid_key // ' values must be strictly increasing'
      return
    end if
    call move_alloc(rarray, this%grid)

    n = n - 1 ! expected size of the following arrays

    !! The number of subintervals for each coarse grid segment (required)
    call params%get(nint_key, iarray, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (size(iarray) /= n) then
      stat = 1
      errmsg = i_to_c(n) // ' values required for ' // nint_key
      return
    else if (any(iarray < 1)) then
      stat = 1
      errmsg = nint_key // ' values must be > 0'
      return
    end if
    call move_alloc(iarray, this%nint)

    !! The ratio for biased subdivision for each segment (optional)
    if (params%is_parameter(ratio_key)) then
      call params%get(ratio_key, rarray, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      if (size(rarray) /= n) then
        stat = 1
        errmsg = i_to_c(n) // ' values required for ' // ratio_key
        return
      else if (any(rarray <= 0)) then
        stat = 1
        errmsg = ratio_key // ' values must be > 0'
        return
      end if
      call move_alloc(rarray, this%ratio)
    else
      allocate(this%ratio(n))
      this%ratio = 1.0_r8
    end if

  end subroutine coord_grid_init

  !! Given the arrays GRID(1:N+1), NINT(1:N), and RATIO(1:N) this subroutine
  !! discretizes the interval [GRID(1), GRID(N+1)], returning the results in
  !! the array X(0:M), where M = SUM(NINT). The segment [GRID(j), GRID(j+1)]
  !! is divided into NINT(j) subintervals, with the ratio of the lengths of
  !! successive subintervals equal to RATIO(j).

  subroutine get_grid(this, x)

    class(coord_grid), intent(in) :: this
    real(r8), allocatable, intent(out) :: x(:)

    integer :: j, n1, n2

    allocate(x(0:sum(this%nint)))

    n2 = 0
    x(0) = this%grid(1)
    do j = 1, size(this%nint)
      n1 = n2
      n2 = n1 + this%nint(j)
      x(n2) = this%grid(j+1)
      call discretize_interval(this%ratio(j), x(n1:n2))
    end do

  end subroutine get_grid

  !! Given an array X(0:N) this auxiliary subroutine discretizes the interval
  !! [X(0), X(N)], defining the intermediate values X(j) for j = 1 to N-1. The
  !! ratio of the lengths of successive subintervals is RATIO.

  subroutine discretize_interval(ratio, x)
    real(r8), intent(in) :: ratio
    real(r8), intent(inout) :: x(0:)
    integer :: j, n
    real(r8) :: a, b
    n = ubound(x,1)
    a = 1
    do j = 1, n-1
      a = ratio*a + 1
    end do
    a = (x(n)-x(0))/a
    b = 1
    do j = 1, n-1
      x(j) = x(0) + a*b
      b = ratio*b + 1
    end do
  end subroutine discretize_interval

end module exodus_mesh_factory
