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

    use common_impl, only: nnodes, pnode_old_to_new, gnmap, mesh => new_mesh
    use parallel_permutations, only: rearrange
    use cutoffs_module, only: alittle

    integer :: j, k
    real(r8) :: src(mesh%nnode)

    allocate(vertex(nnodes))

    do k = 1, 3
      call rearrange (pnode_old_to_new, vertex%coord(k), mesh%x(k,:mesh%nnode_onP))
    end do

    src = 0.0_r8
    do j = 1, mesh%ncell
      associate (cnode => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
        src(cnode) = src(cnode) + 1.0_r8 / mesh%volume(j)
      end associate
    end do
    do j = 1, mesh%nlink
      if (mesh%link_cell_id(j) == 0) cycle ! not from a gap cell
      associate (lnode => mesh%lnode(mesh%xlnode(j):mesh%xlnode(j+1)-1))
        src(lnode) = src(lnode) + 1.0_r8 / (2*alittle)
      end associate
    end do
    call gnmap%sum_to_parent(src)
    src = 1.0_r8 / src
    call rearrange (pnode_old_to_new, vertex%rsum_rvol, src(:mesh%nnode_onP))

    !call test_vertex ! If used, see remarks below (will fail for gap cells)

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
      call collate (vertex_collate%coord(n), vertex%coord(n))
    end do
    call collate (vertex_collate%rsum_rvol, vertex%rsum_rvol)
  end function vertex_collate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TESTING CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Testing of RSUM_RVOL values is difficult for several reasons.  First, we
  !! are computing cell volumes using an entirely different method and so they
  !! alone differ by some small hard-to-quantify amount.  But a more serious
  !! problem occurs when gap elements are present.  These are assigned an
  !! arbitrarily small volume of 2*alittle, which completely dominates the
  !! computation of RSUM_RVOL and lead to a result comparable to 2*alittle
  !! which is nigh impossible to match to the required relative precision.
  !! This is almost certainly *not* what was intended -- gap cells should
  !! be ignored in this case -- but it is what Truchas currently does.

  subroutine test_vertex

    use mesh_module, only: old_vertex => vertex
    use parallel_communication, only: global_all
    use truchas_logging_services

    integer :: j, k, unit
    logical :: error

    do k = 1, 3
      INSIST(global_all(vertex%coord(k) == old_vertex%coord(k)))
    end do

#define DEBUG
#ifdef DEBUG
#undef DEBUG
    unit = TLS_debug_unit()
    error = .false.
    do j = 1, size(vertex)
      if (abs(vertex(j)%rsum_rvol - old_vertex(j)%rsum_rvol) <= 128*spacing(old_vertex(j)%rsum_rvol)) cycle
      if (.not.error) then
        write(unit, '(/,a)') 'VERTEX%RSUM_RVOL **********'
        error = .true.
      end if
      write(unit,'(a,i0,a,2(1x,g0))') 'rsum_rvol[', j, ']=', vertex(j)%rsum_rvol, old_vertex(j)%rsum_rvol
    end do
    if (error) call TLS_fatal ('error validating vertex%rsum_rvol')
#else
    INSIST(global_all(abs(vertex%rsum_rvol - old_vertex%rsum_rvol) <= 128*spacing(old_vertex%rsum_rvol)))
#endif

  end subroutine test_vertex

end module vertex_impl
