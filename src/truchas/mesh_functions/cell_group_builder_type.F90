module cell_group_builder_type

  use base_mesh_class
  implicit none
  private

  type, public :: cell_group_builder
    private
    class(base_mesh), pointer :: mesh => null() ! reference only -- not owned
    integer :: ngroup = 0
    integer, allocatable :: tag(:)
    logical, allocatable :: mask(:) ! work space for add_group
  contains
    procedure :: init
    procedure :: add_cell_group
    procedure :: get_cell_groups
  end type cell_group_builder

contains

  subroutine init(this, mesh)
    class(cell_group_builder), intent(out) :: this
    class(base_mesh), target :: mesh
    this%mesh => mesh
    this%ngroup = 0
    allocate(this%tag(mesh%ncell), this%mask(mesh%ncell))
    this%tag = 0
  end subroutine init

  subroutine add_cell_group(this, setids, stat, errmsg)

    use parallel_communication, only: global_count
    use string_utilities, only: i_to_c

    class(cell_group_builder), intent(inout) :: this
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    integer(kind(this%mesh%cell_set_mask)) :: bitmask

    !! Identify the cells specified by SETIDS.
    call this%mesh%get_cell_set_bitmask(setids, bitmask, stat, errmsg)
    if (stat /= 0) return
    this%mask = (popcnt(iand(bitmask, this%mesh%cell_set_mask)) /= 0)

    !! Check that these cells don't overlap those from preceding calls.
    n = global_count(this%mask .and. this%tag /= 0)
    if (n /= 0) then
      stat = 1
      errmsg = i_to_c(n) // ' cells belong to an existing group'
      return
    end if

    !! Set the tag array.
    this%ngroup = 1 + this%ngroup
    where (this%mask) this%tag = this%ngroup

  end subroutine add_cell_group

  subroutine get_cell_groups(this, ngroup, xgroup, index)

    class(cell_group_builder), intent(in) :: this
    integer, intent(out) :: ngroup
    integer, allocatable, intent(out) :: xgroup(:), index(:)

    integer :: n, j

    ngroup = this%ngroup
    n = count(this%tag > 0)
    allocate(index(n), xgroup(ngroup+1))

    !! Prepare XGROUP: indices in group N will be INDEX(XGROUP(N):XGROUP(N+1)-1).
    xgroup(1) = 1
    do n = 1, ngroup
      xgroup(n+1) = xgroup(n) + count(this%tag == n)
    end do

    !! Fill the INDEX array; XGROUP(N) stores the next free location for group N.
    do j = 1, size(this%tag)
      n = this%tag(j)
      if (n == 0) cycle
      index(xgroup(n)) = j
      xgroup(n) = 1 + xgroup(n)
    end do

    !! Restore XGROUP; XGROUP(N) is now the start of group N+1 instead of N.
    do n = ngroup, 1, -1
      xgroup(n+1) = xgroup(n)
    end do
    xgroup(1) = 1

  end subroutine get_cell_groups

end module cell_group_builder_type
