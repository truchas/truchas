!!
!! STRING_SET_TYPE
!!
!! This module defines a set container that stores unique character string
!! values. This is very rudimentary with only the capabilities needed for the
!! moment. Only to be used for small sets and not in perfomance sensitive
!! contexts.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module string_set_type

  implicit none
  private

  type, public :: string_set
    private
    type(list_item), pointer :: first => null()
  contains
    procedure :: add
    procedure :: has
    final :: string_set_delete
  end type

  type :: list_item
    character(:), allocatable :: string
    type(list_item), pointer :: next => null()
  contains
    final :: list_item_delete
  end type

contains

  subroutine string_set_delete(this)
    type(string_set), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine

  recursive subroutine list_item_delete(this)
    type(list_item), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  subroutine add(this, string)
    class(string_set), intent(inout) :: this
    character(*), intent(in) :: string
    type(list_item), pointer :: item
    item => find_list_item(this, string)
    if (associated(item)) return
    allocate(item)
    item%string = string
    item%next => this%first
    this%first => item
  end subroutine

  logical function has(this, string)
    class(string_set), intent(in) :: this
    character(*), intent(in) :: string
    has = associated(find_list_item(this, string))
  end function

  function find_list_item(this, string) result(item)
    class(string_set), intent(in) :: this
    character(*), intent(in) :: string
    type(list_item), pointer :: item
    item => this%first
    do while (associated(item))
      if (item%string == string) return
      item => item%next
    end do
  end function

end module string_set_type
