!!
!! SCALAR_CELL_FUNC1_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scalar_cell_func1_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use base_mesh_class
  use scalar_mesh_func_class
  use scalar_func_containers
  implicit none
  private

  type :: array_list
    type(array_item), pointer :: first => null()
  contains
    procedure :: append => array_list_append
    final :: array_list_delete
  end type

  type :: array_item
    real(r8), allocatable :: array(:)
    type(array_item), pointer :: next => null()
  end type

  type, extends(scalar_mesh_func), public :: scalar_cell_func1
    class(base_mesh), pointer :: mesh => null() ! reference only -- not owned
    type(scalar_func_box), allocatable :: f(:)
    real(r8), allocatable :: array(:,:), w(:)
    ! construction phase temporaries
    integer, allocatable :: perm(:)
    real(r8), allocatable :: global_array(:)
    type(scalar_func_list) :: flist
    type(array_list) :: alist
  contains
    procedure :: init
    procedure :: add
    procedure :: assemble
    procedure :: compute
  end type scalar_cell_func1

contains

  subroutine init(this, mesh)
    use parallel_communication, only: is_IOP, global_sum, collate
    class(scalar_cell_func1), intent(out) :: this
    class(base_mesh), intent(in), target :: mesh
    integer :: n
    this%mesh => mesh
    n = global_sum(mesh%ncell_onP)
    allocate(this%perm(merge(n,0,is_IOP)))
    call collate(mesh%xcell(:mesh%ncell_onP), this%perm)
    allocate(this%global_array(size(this%perm)))
  end subroutine

  subroutine add(this, file, f, stat, errmsg)

    use parallel_communication, only: is_IOP, broadcast
    use permutations, only: reorder

    class(scalar_cell_func1), intent(inout) :: this
    character(*), intent(in) :: file
    class(scalar_func), allocatable, intent(inout) :: f
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: lun
    character(255) :: iom
    real(r8), allocatable :: array(:)

    !! Open the file for unformatted stream reading
    if (is_IOP) open(newunit=lun,file=file,access='stream',status='old',action='read',iostat=stat,iomsg=iom)
    call broadcast(stat)
    if (stat /= 0) then
      call broadcast(iom)
      errmsg = 'error opening file "' // file // '": ' // trim(iom)
      return
    end if

    !! Read the cell-based data array
    if (is_IOP) then
      read(lun,iostat=stat,iomsg=iom) this%global_array
      close(lun)
    end if
    call broadcast(stat)
    if (stat /= 0) then
      call broadcast(iom)
      errmsg = 'error reading file "' // file // '": ' // trim(iom)
      return
    end if

    !! Permute data to internal cell ordering, distribute, and sync ghosts
    if (is_IOP) call reorder(this%global_array, this%perm)
    allocate(array(this%mesh%ncell))
    call this%mesh%cell_imap%distribute(this%global_array, array)
    call this%mesh%cell_imap%gather_offp(array)

    call this%flist%append(f)
    call this%alist%append(array)

  end subroutine add


  subroutine assemble(this)

    use scalar_func_tools, only: is_const

    class(scalar_cell_func1), intent(inout) :: this

    integer :: i
    logical :: all_const
    real(r8) :: c
    type(array_item), pointer :: rest

    call scalar_func_list_to_box_array(this%flist, this%f)

    all_const = .true.
    do i = 1, size(this%f)
      all_const = all_const .and. is_const(this%f(i)%f)
    end do

    if (all_const) then ! go ahead and compute the final value array
      allocate(this%value(this%mesh%ncell))
      this%value = 0
      do i = 1, size(this%f)
        c = this%f(i)%eval([real(r8)::])
        this%value = this%value + c * this%alist%first%array
        rest => this%alist%first%next
        deallocate(this%alist%first)
        this%alist%first => rest
      end do
      deallocate(this%f) ! no longer needed; indicates no computation needed
    else ! move array list into a packed rank-2 array
      allocate(this%array(size(this%f),this%mesh%ncell))
      do i = 1, size(this%array,dim=1)
        this%array(i,:) = this%alist%first%array
        rest => this%alist%first%next
        deallocate(this%alist%first)
        this%alist%first => rest
      end do
    end if

  end subroutine assemble


  subroutine compute(this, t)

    class(scalar_cell_func1), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i, j
    real(r8), allocatable :: w(:)

    if (.not.allocated(this%f)) return ! nothing to do

    !! Evaluate the weights at the given time
    allocate(w(size(this%f)))
    do i = 1, size(this%f)
      w(i) = this%f(i)%eval([t])
    end do

    if (allocated(this%value)) then
      if (all(w == this%w)) return  ! no change; nothing to do
    else
      allocate(this%value(size(this%array,dim=2)))
    end if

    do j = 1, size(this%array,dim=2)
      this%value(j) = dot_product(w, this%array(:,j))
    end do

    this%w = w  ! save for comparison in next call

  end subroutine compute

  !! Final subroutine for ARRAY_LIST objects
  subroutine array_list_delete(this)
    type(array_list), intent(inout) :: this
    type(array_item), pointer :: rest
    do while (associated(this%first))
      rest => this%first%next
      deallocate(this%first)
      this%first => rest
    end do
  end subroutine

  subroutine array_list_append(this, array)
    class(array_list), intent(inout) :: this
    real(r8), allocatable, intent(inout) :: array(:)
    type(array_item), pointer :: last
    if (associated(this%first)) then
      last => this%first
      do while (associated(last%next))
        last => last%next
      end do
      allocate(last%next)
      last => last%next
    else
      allocate(this%first)
      last => this%first
    end if
    call move_alloc(array, last%array)
  end subroutine

end module scalar_cell_func1_type
