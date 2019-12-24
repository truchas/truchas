!!
!! MATERIAL_DATABASE_TYPE
!!
!! This module provides a container for a collection of named MATERIAL class
!! objects. The implementation is very basic. MATERIAL class objects can be
!! added to the container, the container queried whether it contains a named
!! material, and requested for a pointer to a named material.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module material_database_type

  use material_class
  implicit none
  private

  type, public :: material_database
    private
    type(list_item), pointer :: first => null()
  contains
    procedure :: has_matl
    procedure :: matl_ref
    procedure :: add_matl
    final :: material_database_delete
  end type

  type :: list_item
    class(material), allocatable :: matl
    type(list_item), pointer :: next => null()
  contains
    final :: list_item_delete
  end type

contains

  !! Return true if the map contains the NAMEd material.
  logical function has_matl(this, name)
    class(material_database), intent(in) :: this
    character(*), intent(in) :: name
    has_matl = associated(find_list_item(this, name))
  end function

  !! Return a reference to the NAMEd material if it exists; otherwise NULL().
  function matl_ref(this, name)
    class(material_database), intent(in) :: this
    character(*), intent(in) :: name
    class(material), pointer :: matl_ref
    type(list_item), pointer :: item
    matl_ref => null()
    item => find_list_item(this, name)
    if (associated(item)) matl_ref => item%matl
  end function

  !! Add MATL to the map. The object is moved, not copied.
  !! The map must not already contain a material with the same name.
  subroutine add_matl(this, matl)
    class(material_database), intent(inout) :: this
    class(material), allocatable, intent(inout) :: matl
    type(list_item), pointer :: tail
    INSIST(.not.has_matl(this, matl%name))
    tail => this%first
    allocate(this%first)
    call move_alloc(matl, this%first%matl)
    this%first%next => tail
  end subroutine

  !! Return a pointer to the LIST_ITEM containing the NAMEd material; or NULL().
  function find_list_item(this, name) result(item)
    class(material_database), intent(in) :: this
    character(*), intent(in) :: name
    type(list_item), pointer :: item
    item => this%first
    do while (associated(item))
      if (item%matl%name == name) exit
      item => item%next
    end do
  end function

  !! Final subroutine for MATERIAL_DATABASE objects.
  subroutine material_database_delete(this)
    type(material_database), intent(inout) :: this
    if (associated(this%first)) deallocate(this%first)
  end subroutine

  !! Final subroutine for LIST_ITEM objects. This recursively follows the
  !! NEXT pointer. To deallocate a linked-list only the root needs to be
  !! be explicitly deallocated. To deallocate a single LIST_ITEM object,
  !! first nullify its NEXT pointer.
  recursive subroutine list_item_delete(this)
    type(list_item), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

end module material_database_type
