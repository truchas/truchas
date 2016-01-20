#include "f90_assert.fpp"

module gap_node_map_type

  use kinds, only: r8
  use pgslib_module
  implicit none
  private

  type, public :: gap_node_map
    private
    logical :: nop = .true.
    integer :: nnode = 0
    integer, allocatable :: gnode(:)  ! list of (added) gap nodes
    integer, allocatable :: pnode(:)  ! corresponding list of parent nodes
    integer, allocatable :: dup_index(:)
    type(pgslib_gs_trace), pointer :: trace => null()
  contains
    procedure :: init
    procedure, private :: copy_from_parent_real64
    procedure, private :: copy_from_parent_int32
    procedure, private :: copy_from_parent_log32
    generic :: copy_from_parent => copy_from_parent_real64, copy_from_parent_int32, copy_from_parent_log32
    procedure, private :: sum_to_parent_real64
    generic :: sum_to_parent => sum_to_parent_real64
    procedure :: or_to_parent
    procedure, private :: min_to_parent_real64
    generic :: min_to_parent => min_to_parent_real64
    procedure, private :: max_to_parent_real64
    generic :: max_to_parent => max_to_parent_real64
    final :: gap_node_map_delete
  end type gap_node_map

contains

  subroutine gap_node_map_delete (this)
    type(gap_node_map), intent(inout) :: this
    if (associated(this%trace)) call pgslib_deallocate_trace (this%trace)
  end subroutine gap_node_map_delete


  subroutine init (this, mesh)

    use unstr_mesh_type
    use parallel_communication, only: is_IOP, collate, distribute
    use permutations, only: invert_perm

    class(gap_node_map), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh

    integer :: j, n
    integer, allocatable :: perm(:), pnode_g(:), pnode_l(:)

    if (.not.allocated(mesh%parent_node)) return

!    !! PARENT_NODE IDs are external IDs; we need internal IDs :-/
!    !! It is most straightforward to do the conversion in serial.
!    allocate(perm(merge(mesh%node_ip%global_size(),0,is_IOP)))
!    call collate (perm, mesh%xnode(:mesh%nnode_onP))
!    if (is_IOP) call invert_perm (perm) ! external-to-internal node ID
!    allocate(pnode_g(size(perm)))
!    call collate (pnode_g, mesh%parent_node(:mesh%nnode_onP))
!    do j = 1, size(pnode_g)
!      pnode_g(j) = perm(pnode_g(j))
!    end do
!    allocate(pnode_l(mesh%nnode_onP))
!    call distribute (pnode_l, pnode_g)
!    !! PNODE_L now has the internal IDs we need; what PARENT_NODE should be.
    pnode_l = mesh%parent_node
    !! Generate list of gap nodes
    n = 0
    do j = 1, mesh%nnode_onP
      if (pnode_l(j) /= mesh%node_ip%global_index(j)) n = n + 1
    end do
    allocate(this%gnode(n))
    n = 0
    do j = 1, mesh%nnode_onP
      if (pnode_l(j) /= mesh%node_ip%global_index(j)) then
        n = n + 1
        this%gnode(n) = j
      end if
    end do

    !! The parents of the gap nodes.
    this%pnode = pnode_l(this%gnode)

    !! Generate the PGSLib communication trace for accessing the parents
    this%trace => pgslib_setup_trace(this%pnode, mesh%nnode_onP)
    this%dup_index = pgslib_dup_index(this%trace)

    this%nnode = mesh%nnode ! for checking that data arrays are the right size
    this%nop = .false.

  end subroutine init

  subroutine copy_from_parent_real64 (this, array)

    class(gap_node_map), intent(in) :: this
    real(r8), intent(inout) :: array(:)

    integer :: j
    real(r8), allocatable :: sup_buffer(:), dup_buffer(:)

    if (this%nop) return

    ASSERT(size(array) == this%nnode)

    dup_buffer = array(this%dup_index)
    sup_buffer = pgslib_gather_buffer(dup_buffer, this%trace)

    do j = 1, size(this%gnode)
      if (this%pnode(j) > 0) then
        array(this%gnode(j)) = array(this%pnode(j))
      else
        array(this%gnode(j)) = sup_buffer(-this%pnode(j))
      end if
    end do

  end subroutine copy_from_parent_real64

  subroutine copy_from_parent_log32 (this, array)

    class(gap_node_map), intent(in) :: this
    logical, intent(inout) :: array(:)

    integer :: j
    logical, allocatable :: sup_buffer(:), dup_buffer(:)

    if (this%nop) return

    ASSERT(size(array) == this%nnode)

    dup_buffer = array(this%dup_index)
    sup_buffer = pgslib_gather_buffer(dup_buffer, this%trace)

    do j = 1, size(this%gnode)
      if (this%pnode(j) > 0) then
        array(this%gnode(j)) = array(this%pnode(j))
      else
        array(this%gnode(j)) = sup_buffer(-this%pnode(j))
      end if
    end do

  end subroutine copy_from_parent_log32

  subroutine copy_from_parent_int32 (this, array)

    class(gap_node_map), intent(in) :: this
    integer, intent(inout) :: array(:)

    integer :: j
    integer, allocatable :: sup_buffer(:), dup_buffer(:)

    if (this%nop) return

    ASSERT(size(array) == this%nnode)

    dup_buffer = array(this%dup_index)
    sup_buffer = pgslib_gather_buffer(dup_buffer, this%trace)

    do j = 1, size(this%gnode)
      if (this%pnode(j) > 0) then
        array(this%gnode(j)) = array(this%pnode(j))
      else
        array(this%gnode(j)) = sup_buffer(-this%pnode(j))
      end if
    end do

  end subroutine copy_from_parent_int32

  subroutine sum_to_parent_real64 (this, array)

    class(gap_node_map), intent(in) :: this
    real(r8), intent(inout) :: array(:)

    integer :: j, n
    real(r8), allocatable :: sup_buffer(:), dup_buffer(:)

    if (this%nop) return

    ASSERT(size(array) == this%nnode)

    allocate(sup_buffer(pgslib_size_of_sup(this%trace)))
    sup_buffer = 0.0_r8
    do j = 1, size(this%gnode)
      n = this%pnode(j)
      if (n > 0) then
        array(n) = array(n) + array(this%gnode(j))
      else
        sup_buffer(-n) = sup_buffer(-n) + array(this%gnode(j))
      end if
    end do

    dup_buffer = pgslib_scatter_buffer(sup_buffer, this%trace)

    do j = 1, size(this%dup_index)
      n = this%dup_index(j)
      array(n) = array(n) + dup_buffer(j)
    end do

  end subroutine sum_to_parent_real64

  subroutine min_to_parent_real64 (this, array)

    class(gap_node_map), intent(in) :: this
    real(r8), intent(inout) :: array(:)

    integer :: j, n
    real(r8), allocatable :: sup_buffer(:), dup_buffer(:)

    if (this%nop) return

    ASSERT(size(array) == this%nnode)

    allocate(sup_buffer(pgslib_size_of_sup(this%trace)))
    sup_buffer = huge(1.0_r8)
    do j = 1, size(this%gnode)
      n = this%pnode(j)
      if (n > 0) then
        array(n) = min(array(n), array(this%gnode(j)))
      else
        sup_buffer(-n) = min(sup_buffer(-n), array(this%gnode(j)))
      end if
    end do

    dup_buffer = pgslib_scatter_buffer(sup_buffer, this%trace)

    do j = 1, size(this%dup_index)
      n = this%dup_index(j)
      array(n) = min(array(n), dup_buffer(j))
    end do

  end subroutine min_to_parent_real64

  subroutine max_to_parent_real64 (this, array)

    class(gap_node_map), intent(in) :: this
    real(r8), intent(inout) :: array(:)

    integer :: j, n
    real(r8), allocatable :: sup_buffer(:), dup_buffer(:)

    if (this%nop) return

    ASSERT(size(array) == this%nnode)

    allocate(sup_buffer(pgslib_size_of_sup(this%trace)))
    sup_buffer = -huge(1.0_r8)
    do j = 1, size(this%gnode)
      n = this%pnode(j)
      if (n > 0) then
        array(n) = max(array(n), array(this%gnode(j)))
      else
        sup_buffer(-n) = max(sup_buffer(-n), array(this%gnode(j)))
      end if
    end do

    dup_buffer = pgslib_scatter_buffer(sup_buffer, this%trace)

    do j = 1, size(this%dup_index)
      n = this%dup_index(j)
      array(n) = max(array(n), dup_buffer(j))
    end do

  end subroutine max_to_parent_real64

  subroutine or_to_parent (this, array)

    class(gap_node_map), intent(in) :: this
    logical, intent(inout) :: array(:)

    integer :: j, n
    logical, allocatable :: sup_buffer(:), dup_buffer(:)

    if (this%nop) return

    ASSERT(size(array) == this%nnode)

    allocate(sup_buffer(pgslib_size_of_sup(this%trace)))
    sup_buffer = .false.
    do j = 1, size(this%gnode)
      n = this%pnode(j)
      if (n > 0) then
        array(n) = array(n) .or. array(this%gnode(j))
      else
        sup_buffer(-n) = sup_buffer(-n) .or. array(this%gnode(j))
      end if
    end do

    dup_buffer = pgslib_scatter_buffer(sup_buffer, this%trace)

    do j = 1, size(this%dup_index)
      n = this%dup_index(j)
      array(n) = array(n) .or. dup_buffer(j)
    end do

  end subroutine or_to_parent

end module gap_node_map_type
