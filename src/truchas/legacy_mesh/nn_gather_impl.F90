!!
!! NN_GATHER_IMPL
!!
!! Implements the legacy node-node gather procedure.
!!
!! NOTES
!!
!! The public data and nn_gather_boundarydata procedure are only used by solid
!! mechanics, and thus only needs to work for that use case.  That means the
!! data only needs to account for links that came from gap elements, which it
!! does, but not the other links that do not correspond to any legacy mesh cell.
!!

#include "f90_assert.fpp"

module nn_gather_impl

  use kinds, only: r8
  use var_vector_module, only: int_var_vector
  implicit none
  private

  public :: init_nn_gather_impl
  public :: nn_gather_boundarydata

  type(int_var_vector), allocatable, public :: vertex_ngbr_all(:)
  type(int_var_vector), allocatable, public :: vertex_ngbr_all_orig(:)

contains

  subroutine init_nn_gather_impl
    use common_impl, only: nnodes
    allocate(vertex_ngbr_all(nnodes), vertex_ngbr_all_orig(nnodes))
    call init_vertex_ngbr_all(vertex_ngbr_all_orig, vertex_ngbr_all)
  end subroutine init_nn_gather_impl

  subroutine nn_gather_boundarydata (boundary, source)
    use common_impl, only: mesh => new_mesh
    real(r8), intent(in) :: source(:)
    real(r8), pointer :: boundary(:)
    ASSERT(size(source) == mesh%node_imap%onp_size)
    if (associated(boundary)) return  ! signal to do nothing
    allocate(boundary(mesh%node_imap%offp_size))
    call mesh%node_imap%gather_offp(source, boundary)
  end subroutine nn_gather_boundarydata

  subroutine init_vertex_ngbr_all (vertex_ngbr_all_orig, vertex_ngbr_all)

    use common_impl, only: nnodes, mesh => new_mesh
    use var_vector_module, only: create, sizes, flatten

    type(int_var_vector), intent(out) :: vertex_ngbr_all_orig(:), vertex_ngbr_all(:)

    integer, pointer :: container(:)

    ASSERT(size(vertex_ngbr_all) == nnodes)
    ASSERT(size(vertex_ngbr_all_orig) == nnodes)

    call init_vertex_ngbr_all_aux (mesh, vertex_ngbr_all)

    !! Generate VERTEX_NGBR_ALL_ORIG from VERTEX_NGBR_ALL
    call create (vertex_ngbr_all_orig, sizes(vertex_ngbr_all))
    container => flatten(vertex_ngbr_all_orig)
    container = mesh%node_imap%global_index(flatten(vertex_ngbr_all))

    !! Translate off-process node references to boundary buffer references.
    container => flatten(vertex_ngbr_all)
    where (container > nnodes) container = nnodes - container
    INSIST(all(container /= 0))

  end subroutine init_vertex_ngbr_all

  subroutine init_vertex_ngbr_all_aux (mesh, vertex_ngbr_all)

    use unstr_mesh_type
    use integer_set_type
    use var_vector_module, only: create
    use cell_topology, only: cell_edges, link_edges

    type(unstr_mesh), intent(in) :: mesh
    type(int_var_vector), intent(out) :: vertex_ngbr_all(:)

    integer :: j, k, n1, n2
    integer, pointer :: edges(:,:)
    type(integer_set) :: nset(mesh%nnode)

    ASSERT(size(vertex_ngbr_all) == mesh%nnode_onP)

    do j = 1, mesh%ncell
      associate (cnode => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
        edges => cell_edges(cnode)
        do k = 1, size(edges,dim=2)
          n1 = cnode(edges(1,k))
          n2 = cnode(edges(2,k))
          call nset(n1)%add(n2)
          call nset(n2)%add(n1)
        end do
      end associate
    end do

    do j = 1, mesh%nlink
      if (mesh%link_cell_id(j) == 0) cycle  ! not from a gap cell
      associate (lnode => mesh%lnode(mesh%xlnode(j):mesh%xlnode(j+1)-1))
        edges => link_edges(lnode)
        do k =1, size(edges,dim=2)
          n1 = lnode(edges(1,k))
          n2 = lnode(edges(2,k))
          call nset(n1)%add(n2)
          call nset(n2)%add(n1)
        end do
      end associate
    end do

    call create (vertex_ngbr_all, sizes=nset%size())

    do j = 1, mesh%nnode_onP
      call nset(j)%copy_to_array (vertex_ngbr_all(j)%v)
    end do

  end subroutine init_vertex_ngbr_all_aux

end module nn_gather_impl
