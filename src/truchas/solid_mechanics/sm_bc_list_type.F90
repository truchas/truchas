!!
!! SM_BC_LIST_TYPE
!!
!! Solid mechanics BCs are very coupled; they can't be independently applied
!! such as those in other physics. This is because despite being
!! defined along face sets, they are applied at nodes. The specific condition
!! applied at a node is determined by that node's neighboring faces. If those
!! neighboring faces require multiple different conditions, then that node
!! gets a unique condition determined by that combination of neighboring BCs.
!! This is not just a matter of adding or sequentially applying BCs; there is
!! unique logic depending on the combination.
!!
!! Thus there are multiple stages to setting up SM BC entities. We need to
!! read all BCs into one data structure that identifies and indexes them all
!! (this one). Then we build a structure which identifies, for each boundary
!! face, which BCs are applied there (there may be multiple on each face). Then
!! a collection of unique node-based BC objects, each designed for a particular
!! combination of face-BCs, identify the nodes for their own intersection set
!! and hold pointers to any functions read into this SM_BC_LIST_TYPE.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! March 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_list_type

  use truchas_logging_services
  use scalar_func_class
  implicit none
  private

  type :: smbcl_displacement
    integer :: setid, type
    class(scalar_func), allocatable :: f
  end type smbcl_displacement

  type :: smbcl_contact
    integer :: setid
  end type smbcl_contact

  type, public :: sm_bc_list
    integer :: dsz ! displacement size for indexing contact
    type(smbcl_displacement), allocatable :: displacement(:)
    type(smbcl_contact), allocatable :: contact(:)
    integer, allocatable :: bc_type(:)
    !type(smbcl_displacement), allocatable :: traction(:)
  contains
    procedure :: init
  end type sm_bc_list

  integer, parameter, public :: SMBCL_N = 0
  integer, parameter, public :: SMBCL_X = 1
  integer, parameter, public :: SMBCL_Y = 2
  integer, parameter, public :: SMBCL_Z = 3

contains

  subroutine init(this, params, stat, errmsg)

    use parameter_list_type
    use scalar_func_class
    use scalar_func_factories, only: alloc_scalar_func
    use string_utilities, only: lower_case

    class(sm_bc_list), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n

    n = count_entries('displacement-n', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('displacement-x', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('displacement-y', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('displacement-z', stat, errmsg)
    if (stat /= 0) return
    allocate(this%displacement(n))
    this%dsz = n
    n = count_entries('gap-contact', stat, errmsg)
    if (stat /= 0) return
    allocate(this%contact(n))

    n = 0
    ! do t = SMBCL_N, SMBCL_Z
    !   call init_displacement(t, n)
    !   if (stat /= 0) return
    ! end do
    call init_displacement(SMBCL_N, n)
    if (stat /= 0) return
    call init_displacement(SMBCL_X, n)
    if (stat /= 0) return
    call init_displacement(SMBCL_Y, n)
    if (stat /= 0) return
    call init_displacement(SMBCL_Z, n)
    if (stat /= 0) return
    call init_contact
    if (stat /= 0) return

    allocate(this%bc_type(size(this%displacement) + size(this%contact)))
    do n = 1, size(this%bc_type)
      if (n <= this%dsz) then
        this%bc_type(n) = this%displacement(n)%type
      else
        this%bc_type(n) = SMBCL_N
      end if
    end do

  contains

    subroutine init_displacement(type, index)

      integer, intent(in) :: type
      integer, intent(inout) :: index

      character(1), parameter :: dirstr(4) = ['n','x','y','z']
      character(:), allocatable :: this_type, type_string, data_label
      integer, allocatable :: setids(:)
      integer :: s
      class(scalar_func), allocatable :: f
      type(parameter_list_iterator) :: piter
      type(parameter_list), pointer :: plist

      data_label = 'displacement'
      type_string = data_label // '-' // dirstr(type+1)
      stat = 0
      piter = parameter_list_iterator(params, sublists_only=.true.)
      do while (.not.piter%at_end())
        plist => piter%sublist()
        call plist%get('type', this_type, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        if (lower_case(this_type) == type_string) then  ! use this sublist
          call TLS_info('  using SM_BC[' // piter%name() // ']')
          call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
          if (stat /= 0) exit
          call alloc_scalar_func(plist, data_label, f, stat, errmsg)
          if (stat /= 0) exit

          do s = 1, size(setids)
            index = index + 1
            this%displacement(index)%type = type
            this%displacement(index)%setid = setids(s)
            this%displacement(index)%f = f
          end do

        end if
        call piter%next
      end do
      if (stat /= 0) errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg

    end subroutine init_displacement


    subroutine init_contact()

      character(:), allocatable :: this_type, type_string
      integer, allocatable :: setids(:)
      type(parameter_list_iterator) :: piter
      type(parameter_list), pointer :: plist

      type_string = 'gap-contact'
      stat = 0
      piter = parameter_list_iterator(params, sublists_only=.true.)
      do while (.not.piter%at_end())
        plist => piter%sublist()
        call plist%get('type', this_type, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        if (lower_case(this_type) == type_string) then  ! use this sublist
          call TLS_info('  using SM_BC[' // piter%name() // ']')
          call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
          if (stat /= 0) exit
        end if
        call piter%next
      end do
      if (stat /= 0) errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg

    end subroutine init_contact


    !! Count the number of BCs of the given type. Each input face-set-id counts
    !! an independent BC. This is so that, for instance, normal-direction
    !! displacement BCs on different face sets are applied appropriately at any
    !! intersection between those face sets.
    integer function count_entries(type_string, stat, errmsg)

      character(*), intent(in) :: type_string
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      character(:), allocatable :: this_type
      integer, allocatable :: setids(:)
      type(parameter_list_iterator) :: piter
      type(parameter_list), pointer :: plist

      count_entries = 0
      stat = 0
      piter = parameter_list_iterator(params, sublists_only=.true.)
      do while (.not.piter%at_end())
        plist => piter%sublist()
        call plist%get('type', this_type, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        if (lower_case(this_type) == type_string) then  ! use this sublist
          !call TLS_info('  using SM_BC[' // piter%name() // ']')
          call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
          if (stat /= 0) exit
          count_entries = count_entries + size(setids)
        end if
        call piter%next
      end do
      if (stat /= 0) errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg

    end function count_entries

  end subroutine init

end module sm_bc_list_type
