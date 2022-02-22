!!
!! VERTEX_IMPL
!!
!! Implements the legacy VERTEX structure array and associated procedures.
!!

#include "f90_assert.fpp"

module vertex_impl

  use kinds, only: r8
  implicit none
  private

  public :: vertex_data, vertex
  public :: init_vertex_impl, vertex_collate

  !! Copy of type from MESH_MODULE
  type :: vertex_data
     real(r8) :: coord(3)  = 0.0_r8  ! vertex coordinates
     real(r8) :: rsum_rvol = 0.0_r8  ! reciprocal sum of inverse volumes around vertices
  end type vertex_data

  type(vertex_data), allocatable :: vertex(:)

contains

  subroutine init_vertex_impl

    use common_impl, only: nnodes, mesh => new_mesh
    use cutoffs_module, only: alittle

    integer :: j
    real(r8) :: tmp(mesh%nnode)

    allocate(vertex(nnodes))

    do j = 1, nnodes
      vertex(j)%coord = mesh%x(:,j)
    end do

    tmp = 0.0_r8
    do j = 1, mesh%ncell
      associate (cnode => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
        tmp(cnode) = tmp(cnode) + 1.0_r8 / mesh%volume(j)
      end associate
    end do
    do j = 1, mesh%nlink
      if (mesh%link_cell_id(j) == 0) cycle ! not from a gap cell
      associate (lnode => mesh%lnode(mesh%xlnode(j):mesh%xlnode(j+1)-1))
        tmp(lnode) = tmp(lnode) + 1.0_r8 / (2*alittle)
      end associate
    end do
    vertex%rsum_rvol = 1.0_r8 / tmp(:mesh%nnode_onP)

  end subroutine init_vertex_impl

  !! Copy of function from MESH_MODULE
  function vertex_collate (vertex)
    use common_impl, only: nnodes_tot
    use parallel_communication, only: is_IOP, collate
    type(vertex_data), intent(in) :: vertex(:)
    type(vertex_data), pointer :: vertex_collate(:)
    integer :: n
    allocate(vertex_collate(merge(nnodes_tot,0,is_IOP)))
    do n = 1, 3
      call collate (vertex%coord(n), vertex_collate%coord(n))
    end do
    call collate (vertex%rsum_rvol, vertex_collate%rsum_rvol)
  end function vertex_collate

end module vertex_impl
