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
!! A consequence of not dealing with the latter kind of links is that the test
!! subroutine test_vertex_ngbr_all will produce output into the debug files
!! indicating differences in the vertex_ngbr_all structure array.  These can
!! be safely ignored.  The issue with such links is that their neighbors are
!! adjacent in the legacy mesh, but not in the new mesh.
!!

#include "f90_assert.fpp"

module nn_gather_impl

  use kinds, only: r8
  use var_vector_module, only: int_var_vector
  use index_partitioning, only: ip_desc
  implicit none
  private

  public :: init_nn_gather_impl
  public :: nn_gather_boundarydata

  type(int_var_vector), allocatable, public :: vertex_ngbr_all(:)
  type(int_var_vector), allocatable, public :: vertex_ngbr_all_orig(:)

  type(ip_desc) :: node_ip

contains

  subroutine init_nn_gather_impl
    use common_impl, only: nnodes
    allocate(vertex_ngbr_all(nnodes), vertex_ngbr_all_orig(nnodes))
    call init_vertex_ngbr_all_mapped (node_ip, vertex_ngbr_all_orig, vertex_ngbr_all)
  end subroutine init_nn_gather_impl

  subroutine nn_gather_boundarydata (boundary, source)
    use index_partitioning, only: gather_boundary
    real(r8), intent(in) :: source(:)
    real(r8), pointer :: boundary(:)
    ASSERT(size(source) == node_ip%onP_size())
    if (associated(boundary)) return  ! signal to do nothing
    allocate(boundary(node_ip%offP_size()))
    call gather_boundary (node_ip, source, boundary)
  end subroutine nn_gather_boundarydata

  subroutine init_vertex_ngbr_all_mapped (node_ip, vertex_ngbr_all_orig, vertex_ngbr_all)

    use common_impl, only: new_mesh, nnodes, pnode_new_to_old, pnode_old_to_new, gnmap
    use parallel_permutations, only: rearrange
    use index_partitioning, only: gather_boundary, localize_index_struct
    use var_vector_module, only: create, destroy, sizes, flatten
    use sort_module, only: heapsort
    use parallel_communication, only: global_maxval

    type(ip_desc), intent(out) :: node_ip
    type(int_var_vector), intent(out) :: vertex_ngbr_all_orig(:), vertex_ngbr_all(:)

    integer :: j, n
    integer :: old_gid(nnodes), new_gid(new_mesh%nnode)
    integer, allocatable :: new_sizes(:), old_sizes(:), src(:,:), dest(:,:), offP_index(:)
    type(int_var_vector), allocatable :: new_vertex_ngbr_all(:)
    integer, pointer :: container(:)

    ASSERT(size(vertex_ngbr_all) == nnodes)

    !! Push the array of old-mesh global node IDs to the new mesh.
    call node_ip%init (nnodes)
    do j = 1, nnodes
      old_gid(j) = node_ip%global_index(j)
    end do
    call rearrange (pnode_new_to_old, new_gid(:new_mesh%nnode_onP), old_gid, default=0)
    call gnmap%copy_from_parent (new_gid)
    call gather_boundary (new_mesh%node_ip, new_gid)

    !! Generate the VERTEX_NGBR_ALL data relative to the new mesh.
    allocate(new_vertex_ngbr_all(new_mesh%nnode_onP))
    call init_vertex_ngbr_all_aux (new_mesh, new_vertex_ngbr_all)

    !! Map the neighbor data to the legacy mesh global node numbering.
    do j = 1, size(new_vertex_ngbr_all)
      associate (list => new_vertex_ngbr_all(j)%v)
        INSIST(all(list > 0))
        list = new_gid(list)
        call heapsort (list)
      end associate
    end do

    !! Move this data to the old partition.  This is awkward because there is
    !! a variable amount of data per cell and no easy way to do it with PGSLib.
    !! So we unpack into a sufficiently large rank-2 array and move it instead.
    new_sizes = sizes(new_vertex_ngbr_all)
    allocate(old_sizes(nnodes))
    call rearrange (pnode_old_to_new, old_sizes, new_sizes) ! drops new gap node data
    n = global_maxval(new_sizes)
    allocate(src(n,size(new_sizes)), dest(n,nnodes))

    !! Unpack the data into a regular rank-2 array.
    src = 0
    do j = 1, size(new_vertex_ngbr_all)
      n = size(new_vertex_ngbr_all(j)%v)
      src(1:n,j) = new_vertex_ngbr_all(j)%v
    end do
    call destroy (new_vertex_ngbr_all)

    !! Move the data to the old mesh; drops data at new gap nodes.
    call rearrange (pnode_old_to_new, dest, src)

    !! Pack the data into an int_var_vector
    call create (vertex_ngbr_all, old_sizes)
    do j = 1, size(vertex_ngbr_all)
      n = size(vertex_ngbr_all(j)%v)
      vertex_ngbr_all(j)%v = dest(1:n,j)
    end do

    !! Save the neighbor data (global IDs) into another int_var_vector
    !! before it is modified to point to local IDs below.
    call create (vertex_ngbr_all_orig, old_sizes)
    container => flatten(vertex_ngbr_all_orig)
    container = flatten(vertex_ngbr_all)

    !! Localize the VERTEX_NGBR_ALL neighbor indices and augment the node index
    !! partition with off-process node info.
    container => flatten(vertex_ngbr_all)
    call localize_index_struct (node_ip, old_sizes, container, offP_index)
    call node_ip%add_offP_index (offP_index)

    !! Translate off-process node references to boundary buffer references.
    where (container > nnodes) container = nnodes - container
    INSIST(all(container /= 0))

  end subroutine init_vertex_ngbr_all_mapped

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
