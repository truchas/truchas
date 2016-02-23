!!
!! EN_GATHER_IMPL
!!
!! Implements the legacy element-node gather/scatter procedures.
!!

#include "f90_assert.fpp"

module en_gather_impl

  use kinds, only: r8
  use common_impl, only: ncells, nnodes, new_mesh
  implicit none
  private

  public :: init_en_gather_impl
  public :: en_gather, gather_vertex_coord, en_min_gather, en_or_gather
  public :: en_sum_scatter, en_or_scatter, en_min_scatter, en_max_scatter

  interface gather_vertex_coord
    procedure gather_vertex_coord_dim
    procedure gather_vertex_coord_all
  end interface

  interface en_gather
    module procedure en_gather_real64
  end interface

  interface en_min_gather
    procedure en_min_gather_real64
  end interface

  interface en_or_gather
    procedure en_or_gather_log
  end interface

  !! NB: The original en_sum_scatter fails to ignore degenerate nodes for
  !! non-hex cells, resulting in multiple cell values being summed to the same
  !! node.  This cannot be correct (see issue #245).  Nevertheless the current
  !! implementation preserves this behavior.  Be aware that the same scatter
  !! operation done on the new mesh, which does not use degenerate hexes, will
  !! produce different results (likely the correct results).
  !!
  !! NNC, 15 Feb 2016.  Commenting out the validating versions.  Several SM
  !! tests exhibit occasional validation failures with Intel 16, but not
  !! Intel 15 or NAG 6.0.  The differences are always small.  This is very
  !! disturbing, but no time to try and track down the problem now.
  interface en_sum_scatter
    procedure en_sum_scatter_real64
    procedure cn_sum_scatter_real64
  end interface

  interface en_or_scatter
    procedure en_or_scatter_log
    procedure cn_or_scatter_log
  end interface

  interface en_min_scatter
    procedure en_min_scatter_real64
  end interface

  interface en_max_scatter
    procedure cn_max_scatter_real64
  end interface

  real(r8), allocatable :: coord_boundary(:,:)

contains

  subroutine init_en_gather_impl
    use vertex_impl, only: vertex
    use index_partitioning, only: gather_boundary
    integer :: n
    allocate(coord_boundary(3,new_mesh%node_ip%offP_size()))
    do n = 1, 3
      call gather_boundary (new_mesh%node_ip, vertex%coord(n), coord_boundary(n,:))
    end do
  end subroutine init_en_gather_impl

  !! This are replacements for EN_GATHER applied to node coordinates.
  !! These are the calls that used the optional BOUNDARY argument.

  subroutine gather_vertex_coord_dim (coord, dim)

    use mesh_impl, only: mesh
    use vertex_impl, only: vertex

    real(r8), intent(out) :: coord(:,:)
    integer, intent(in) :: dim

    integer :: j, k, n

    ASSERT(size(coord,1) == 8)
    ASSERT(size(coord,2) == ncells)

    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          coord(k,j) = vertex(n)%coord(dim)
        case (:-1)
          coord(k,j) = coord_boundary(dim,-n)
        end select
      end do
    end do

  end subroutine gather_vertex_coord_dim


  subroutine gather_vertex_coord_all (coord)

    use mesh_impl, only: mesh
    use vertex_impl, only: vertex

    real(r8), intent(out) :: coord(:,:,:)

    integer :: j, k, n

    ASSERT(size(coord,1) == 3)
    ASSERT(size(coord,2) == 8)
    ASSERT(size(coord,3) == ncells)

    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          coord(:,k,j) = vertex(n)%coord
        case (:-1)
          coord(:,k,j) = coord_boundary(:,-n)
        end select
      end do
    end do

  end subroutine gather_vertex_coord_all


  subroutine en_gather_real64 (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: gather_boundary

    real(r8), intent(out) :: dest(:,:)
    real(r8), intent(in)  :: src(:)

    integer :: j, k, n
    real(r8) :: boundary(new_mesh%node_ip%offP_size())

    call gather_boundary (new_mesh%node_ip, src, boundary)

    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          dest(k,j) = src(n)
        case (:-1)
          dest(k,j) = boundary(-n)
        case default
          dest(k,j) = 0.0_r8
        end select
      end do
    end do

  end subroutine en_gather_real64


  subroutine en_min_gather_real64 (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: gather_boundary

    real(r8), intent(out) :: dest(:)
    real(r8), intent(in)  :: src(:)

    integer :: j, k, n
    real(r8) :: tmp, boundary(new_mesh%node_ip%offP_size())

    call gather_boundary (new_mesh%node_ip, src, boundary)

    do j = 1, ncells
      tmp = huge(tmp)
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          tmp = min(tmp, src(n))
        case (:-1)
          tmp = min(tmp, boundary(-n))
        end select
      end do
      dest(j) = tmp
    end do

  end subroutine en_min_gather_real64


  subroutine en_or_gather_log (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: gather_boundary

    logical, intent(out) :: dest(:)
    logical, intent(in)  :: src(:)

    integer :: j, k, n
    logical :: tmp, boundary(new_mesh%node_ip%offP_size())

    call gather_boundary (new_mesh%node_ip, src, boundary)

    do j = 1, ncells
      tmp = .false.
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          tmp = (tmp .or. src(n))
        case (:-1)
          tmp = (tmp .or. boundary(-n))
        end select
      end do
      dest(j) = tmp
    end do

  end subroutine en_or_gather_log


  subroutine en_sum_scatter_real64 (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: scatter_boundary_sum

    real(r8), intent(out) :: dest(:)
    real(r8), intent(in)  :: src(:)

    integer  :: j, k, n
    real(r8) :: boundary(new_mesh%node_ip%offP_size())

    ASSERT(size(src) == ncells)
    ASSERT(size(dest) == nnodes)

    dest = 0.0_r8
    boundary = 0.0_r8
    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          dest(n) = dest(n) + src(j)
        case (:-1)
          boundary(-n) = boundary(-n) + src(j)
        end select
      end do
    end do

    call scatter_boundary_sum (new_mesh%node_ip, dest, boundary)

  end subroutine en_sum_scatter_real64


  subroutine cn_sum_scatter_real64 (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: scatter_boundary_sum

    real(r8), intent(out) :: dest(:)
    real(r8), intent(in)  :: src(:,:)

    integer  :: j, k, n
    real(r8) :: boundary(new_mesh%node_ip%offP_size())

    ASSERT(size(src,1) == 8)
    ASSERT(size(src,2) == ncells)
    ASSERT(size(dest) == nnodes)

    dest = 0.0_r8
    boundary = 0.0_r8
    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          dest(n) = dest(n) + src(k,j)
        case (:-1)
          boundary(-n) = boundary(-n) + src(k,j)
        end select
      end do
    end do

    call scatter_boundary_sum (new_mesh%node_ip, dest, boundary)

  end subroutine cn_sum_scatter_real64


  subroutine en_min_scatter_real64 (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: scatter_boundary_min

    real(r8), intent(out) :: dest(:)
    real(r8), intent(in)  :: src(:)

    integer  :: j, k, n
    real(r8) :: boundary(new_mesh%node_ip%offP_size())

    ASSERT(size(src) == ncells)
    ASSERT(size(dest) == nnodes)

    dest = huge(dest)
    boundary = huge(boundary)
    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          dest(n) = min(dest(n), src(j))
        case (:-1)
          boundary(-n) = min(boundary(-n), src(j))
        end select
      end do
    end do

    call scatter_boundary_min (new_mesh%node_ip, dest, boundary)

  end subroutine en_min_scatter_real64


  subroutine cn_max_scatter_real64 (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: scatter_boundary_max

    real(r8), intent(out) :: dest(:)
    real(r8), intent(in)  :: src(:,:)

    integer  :: j, k, n
    real(r8) :: boundary(new_mesh%node_ip%offP_size())

    ASSERT(size(src,1) == 8)
    ASSERT(size(src,2) == ncells)
    ASSERT(size(dest) == nnodes)

    dest = -huge(dest)
    boundary = -huge(boundary)
    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          dest(n) = max(dest(n), src(k,j))
        case (:-1)
          boundary(-n) = max(boundary(-n), src(k,j))
        end select
      end do
    end do

    call scatter_boundary_max (new_mesh%node_ip, dest, boundary)

  end subroutine cn_max_scatter_real64


  subroutine en_or_scatter_log (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: scatter_boundary_or

    logical, intent(out) :: dest(:)
    logical, intent(in)  :: src(:)

    integer :: j, k, n
    logical :: boundary(new_mesh%node_ip%offP_size())

    ASSERT(size(src) == ncells)
    ASSERT(size(dest) == nnodes)

    dest = .false.
    boundary = .false.
    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          dest(n) = (dest(n) .or. src(j))
        case (:-1)
          boundary(-n) = (boundary(-n) .or. src(j))
        end select
      end do
    end do

    call scatter_boundary_or (new_mesh%node_ip, dest, boundary)

  end subroutine en_or_scatter_log


  subroutine cn_or_scatter_log (dest, src)

    use mesh_impl, only: mesh
    use index_partitioning, only: scatter_boundary_or

    logical, intent(out) :: dest(:)
    logical, intent(in)  :: src(:,:)

    integer :: j, k, n
    logical :: boundary(new_mesh%node_ip%offP_size())

    ASSERT(size(src,1) == 8)
    ASSERT(size(src,2) == ncells)
    ASSERT(size(dest) == nnodes)

    dest = .false.
    boundary = .false.
    do j = 1, ncells
      do k = 1, 8
        n = mesh(j)%ngbr_vrtx(k)
        select case (n)
        case (1:)
          dest(n) = (dest(n) .or. src(k,j))
        case (:-1)
          boundary(-n) = (boundary(-n) .or. src(k,j))
        end select
      end do
    end do

    call scatter_boundary_or (new_mesh%node_ip, dest, boundary)

  end subroutine cn_or_scatter_log

end module en_gather_impl
