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

  type :: smbcl_func
    integer :: setid, type
    logical :: nodeset = .false.
    class(scalar_func), allocatable :: f
  end type smbcl_func

  type :: smbcl_contact
    integer :: setid
  end type smbcl_contact

  type, public :: sm_bc_list
    integer :: xcontact ! offset for indexing traction
    integer :: xtraction ! offset for indexing traction
    type(smbcl_func), allocatable :: displacement(:)
    type(smbcl_contact), allocatable :: contact(:)
    type(smbcl_func), allocatable :: traction(:)
    integer, allocatable :: bc_type(:)
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

    integer :: n, t

    n = count_entries('displacement-n', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('displacement-x', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('displacement-y', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('displacement-z', stat, errmsg)
    if (stat /= 0) return
    allocate(this%displacement(n))
    n = count_entries('gap-contact', stat, errmsg)
    if (stat /= 0) return
    allocate(this%contact(n))
    n = count_entries('traction-n', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('traction-x', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('traction-y', stat, errmsg)
    if (stat /= 0) return
    n = n + count_entries('traction-z', stat, errmsg)
    if (stat /= 0) return
    allocate(this%traction(n))

    n = 0
    do t = SMBCL_N, SMBCL_Z
      call init_bc('displacement', t, this%displacement, n)
      if (stat /= 0) return
    end do
    n = 0
    do t = SMBCL_N, SMBCL_Z
      call init_bc('traction', t, this%traction, n)
      if (stat /= 0) return
    end do
    call init_contact
    if (stat /= 0) return

    this%xcontact = 1 + size(this%displacement)
    this%xtraction = this%xcontact + size(this%contact)
    allocate(this%bc_type(size(this%displacement) + size(this%contact) + size(this%traction)))
    this%bc_type(1:this%xcontact-1) = this%displacement(:)%type
    this%bc_type(this%xcontact:this%xtraction-1) = SMBCL_N
    this%bc_type(this%xtraction:) = this%traction(:)%type

    call check_consistency
    if (stat /= 0) return

  contains

    subroutine init_bc(data_label, type, bcdata, index)

      character(*), intent(in) :: data_label
      integer, intent(in) :: type
      type(smbcl_func), intent(inout) :: bcdata(:)
      integer, intent(inout) :: index

      character(1), parameter :: dirstr(4) = ['n','x','y','z']
      character(:), allocatable :: this_type, type_string
      integer, allocatable :: setids(:)
      logical :: nodeset
      integer :: s
      class(scalar_func), allocatable :: f
      type(parameter_list_iterator) :: piter
      type(parameter_list), pointer :: plist

      type_string = data_label // '-' // dirstr(type+1)
      stat = 0
      piter = parameter_list_iterator(params, sublists_only=.true.)
      do while (.not.piter%at_end())
        plist => piter%sublist()
        call plist%get('type', this_type, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        if (lower_case(this_type) == type_string) then  ! use this sublist
          call TLS_info('  using SM_BC[' // piter%name() // ']')
          nodeset = plist%is_parameter('node-set-ids')
          if (nodeset) then
            if (data_label /= 'displacement') then
              stat = 1
              errmsg = 'Given node-set-ids, incompatible with traction type BC.'
              exit
            end if
            call plist%get('node-set-ids', setids, stat=stat, errmsg=errmsg)
          else
            call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
          end if
          if (stat /= 0) exit
          call alloc_scalar_func(plist, data_label, f, stat, errmsg)
          if (stat /= 0) exit

          do s = 1, size(setids)
            index = index + 1
            bcdata(index)%type = type
            bcdata(index)%setid = setids(s)
            bcdata(index)%nodeset = nodeset
            allocate(bcdata(index)%f, source=f)
          end do

        end if
        call piter%next
      end do
      if (stat /= 0) errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg

    end subroutine init_bc


    subroutine init_contact()

      character(:), allocatable :: this_type, type_string
      integer, allocatable :: setids(:)
      type(parameter_list_iterator) :: piter
      type(parameter_list), pointer :: plist
      integer :: s, index

      type_string = 'gap-contact'
      index = 0
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

          do s = 1, size(setids)
            index = index + 1
            this%contact(index)%setid = setids(s)
          end do
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
          if (plist%is_parameter('node-set-ids')) then
            call plist%get('node-set-ids', setids, stat=stat, errmsg=errmsg)
          else
            call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
          end if
          if (stat /= 0) exit
          count_entries = count_entries + size(setids)
        end if
        call piter%next
      end do
      if (stat /= 0) errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg

    end function count_entries


    !! Perform a basic consistency check. Rules & neglected rules:
    !!   - Traction & displacement BCs at the same face set must be
    !!     applied in orthogonal directions. Here we only check if they
    !!     are of different types -- e.g. we don't check that the normal
    !!     of a face set is orthogonal to xhat for all faces in that
    !!     set.
    !!   - Gap contact can't be applied together with normal traction
    !!     or normal displacement.
    !!   - We don't check whether different face sets overlap.
    !!   - We don't check whether face sets and node sets overlap.
    subroutine check_consistency()

      integer :: i, j
      character(256) :: msg

      do i = 1, size(this%displacement)
        if (this%displacement(i)%nodeset) cycle

        ! check if any traction & displacement BCs are applied
        do j = 1, size(this%traction)
          if (this%displacement(i)%setid == this%traction(j)%setid &
              .and. this%displacement(i)%type == this%traction(j)%type) then
            stat = 1
            write(msg,'(a,i5)') "Tangential displacement and traction BCs applied on the same face set. face_set_ids = ", this%displacement(i)%setid
            errmsg = trim(msg)
            return
          end if
        end do
      end do

      ! check if any contact conditions overlap with normal traction or
      ! normal displacement.
      do i = 1, size(this%contact)
        do j = 1, size(this%displacement)
          if (this%contact(i)%setid == this%displacement(j)%setid &
              .and. this%displacement(j)%type == SMBCL_N) then
            stat = 1
            write(msg,'(a,i5)') "Gap contact BC applied to same face set as displacement BC. face_set_ids = ", this%contact(i)%setid
            errmsg = trim(msg)
            return
          end if
        end do

        do j = 1, size(this%traction)
          if (this%contact(i)%setid == this%traction(j)%setid &
              .and. this%traction(j)%type == SMBCL_N) then
            stat = 1
            write(msg,'(a,i5)') "Gap contact BC applied to same face set as traction BC. face_set_ids = ", this%contact(i)%setid
            errmsg = trim(msg)
            return
          end if
        end do
      end do

    end subroutine check_consistency

  end subroutine init

end module sm_bc_list_type
