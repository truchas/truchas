!!
!! MATERIAL_TABLE
!!
!! Maintains a user-defined collection of material systems.
!!
!! PROGRAMMING INTERFACE
!!
!! The module defines a single instance of a material system table that is held
!! as private module data.  All interaction with the table is performed using
!! the procedures described below.  The table is initially empty.  At any point
!! in time the table will contain a user-defined collection of named material
!! systems, implemented as MAT_SYSTEM objects.  Associated with each material
!! name is a unique positive integer ID assigned by the system.  Many of the
!! procedures take the ID associated with the material system name as an
!! argument rather that the name itself.
!!
!! In the following procedures the name arguments are character variables/
!! values whose (trimmed) length does not exceed the value of the public
!! parameter MT_MAX_NAME_LEN.
!!
!!  MT_MATERIAL_ID(NAME) returns the ID assigned to the specified material
!!    system name.  It is an error if NAME is not a defined name.
!!
!!  MT_MATERIAL_NAME(ID) returns the material name corresponding to ID.
!!    It is an error if ID is not assigned to a material.
!!
!!  MT_HAS_MATERIAL(NAME) returns true if NAME is a defined material name;
!!    otherwise it returns false.
!!
!!  MT_VALID_MATERIAL(ID) returns true if ID is a valid material system ID;
!!    otherwise it returns false.
!!
!!  MT_NUM_MATERIAL() returns the number of materials defined.
!!
!!  CALL MT_GET_MATERIAL_IDS (ID) fills the user-supplied rank-1 integer
!!    array ID with all the defined material IDs.  The size of the array
!!    must be at least equal to the number of materials.  The value of any
!!    additional array elements is left unchanged.
!!
!!  CALL MT_GET_MATERIAL_NAMES (NAME) returns an array of all the material
!!    names.  NAME is a pointer to a rank-1 character array that is allocated
!!    by the subroutine with a size equal to the number of materials.  The
!!    character length must be at least large enough to hold the largest
!!    trimmed length of the names, which is no greater than MT_MAX_NAME_LEN.
!!
!!  MT_GET_MATERIAL(ID) returns a pointer to the MAT_SYSTEM object associated
!!    to the specified ID, or a null pointer if the ID is invalid.
!!
!!  CALL MT_ADD_MATERIAL (NAME, MS, ID) adds the material system MS with the
!!    specified NAME to the table and returns the assigned ID.  It is an error
!!    if the material name is already defined.  Note that the table stores an
!!    independent copy of the actual argument MS, which, as a consequence may
!!    be destroyed or modified after making the call if desired.
!!
!!  CALL MT_DELETE_MATERIAL (ID) deletes the material corresponding to ID from
!!    the table.  If the ID is invalid the table is not modified.  Any pointer
!!    associated to the deleted material system is rendered invalid.
!!
!!  CALL MT_RESET_TABLE () deletes all the materials from the table, returning
!!    it to its initial empty state.  Any pointers to the MAT_SYSTEM objects
!!    that were stored in the table are rendered invalid.
!!
!!  CALL MT_DUMP_TABLE () writes the contents of the table to the specified
!!    logical unit in a readable form.  Intended for debugging use only.
!!
!! IMPLEMENTATION NOTES
!!
!! The material systems are expected to number a few dozen at most, and the
!! table is not intended to be accessed intensively during the course of a
!! calculation.  The intent is for application code to cache pointers to the
!! MAT_SYSTEM objects returned by MT_GET_MATERIAL during an initialization phase
!! of a calculation and thereafter make few, if any, accesses to the table.
!! Thus speed of the implementation is not a concern; the table is implemented
!! as an unsorted linked-list of name-material system pairs -- the simplest
!! possible implementation having no a priori size limitation.
!!

#include "f90_assert.fpp"

module material_table

  use material_system
  implicit none
  private

  integer, parameter, public :: MT_MAX_NAME_LEN = 31

  public :: mt_material_id, mt_material_name, mt_has_material, mt_valid_material, mt_get_material
  public :: mt_num_material, mt_get_material_names, mt_get_material_ids
  public :: mt_add_material, mt_delete_material, mt_reset_table
  public :: mt_dump_table

  !! Single instance of the table held as private module data;
  !! All procedures operate on this data.
  type, private :: node
    integer :: id = 0
    character(len=MT_MAX_NAME_LEN) :: name = ''
    type(mat_system) :: ms
    type(node), pointer :: next => null()
  end type node
  type(node), save :: head

contains

  integer function mt_material_id (name) result (id)
    character(len=*), intent(in) :: name
    type(node), pointer :: p
    p => scan_for_name(head, name)
    INSIST(associated(p))
    id = p%id
  end function mt_material_id

  function mt_material_name (id) result (name)
    integer, intent(in) :: id
    character(len=MT_MAX_NAME_LEN) :: name
    type(node), pointer :: p
    p => scan_for_id(head, id)
    INSIST(associated(p))
    name = p%name
  end function mt_material_name

  logical function mt_has_material (name)
    character(len=*), intent(in) :: name
    mt_has_material = associated(scan_for_name(head, name))
  end function mt_has_material

  logical function mt_valid_material (id)
    integer, intent(in) :: id
    mt_valid_material = associated(scan_for_id(head, id))
  end function mt_valid_material

  integer function mt_num_material ()
    type(node), pointer :: p
    mt_num_material = 0
    p => head%next
    do while (associated(p))
      mt_num_material = mt_num_material + 1
      p => p%next
    end do
  end function mt_num_material

  subroutine mt_get_material_names (name)
    character(len=*), pointer :: name(:)
    integer :: j
    type(node), pointer :: p
    allocate(name(mt_num_material()))
    p => head%next
    do j = size(name), 1, -1
      ASSERT(associated(p))
      INSIST(len_trim(p%name) <= len(name))
      name(j) = p%name
      p => p%next
    end do
  end subroutine mt_get_material_names

  subroutine mt_get_material_ids (id)
    integer, intent(inout) :: id(:)
    integer :: j
    type(node), pointer :: p
    ASSERT(size(id) >= mt_num_material())
    p => head%next
    do j = mt_num_material(), 1, -1
      ASSERT(associated(p))
      id(j) = p%id
      p => p%next
    end do
  end subroutine mt_get_material_ids

  function mt_get_material (id) result (ms)
    integer, intent(in) :: id
    type(mat_system), pointer :: ms
    type(node), pointer :: p
    ms => null()
    p => scan_for_id(head, id)
    if (associated(p)) ms => p%ms
  end function mt_get_material

  subroutine mt_add_material (name, ms, id)

    character(len=*), intent(in) :: name
    type(mat_system), intent(in) :: ms
    integer, intent(out) :: id

    type(node), pointer :: p

    INSIST(len_trim(name) > 0 .and. len_trim(name) <= MT_MAX_NAME_LEN)
    INSIST(.not.associated(scan_for_name(head, name)))

    allocate(p)
    head%id = head%id + 1 ! generate new ID
    p%id = head%id
    p%name = name
    p%ms = ms !!! THIS MUST BE A DEEP COPY !!!
    p%next => head%next
    head%next => p

    id = p%id

  end subroutine mt_add_material

  subroutine mt_delete_material (id)
    integer, intent(in) :: id
    type(node), pointer :: p
    p => scan_for_id(head, id)
    call unlink (head, p)
    call deallocate_node (p)
  end subroutine mt_delete_material

  subroutine mt_reset_table ()
    type(node), pointer :: first
    type(node) :: default
    do while (associated(head%next)) ! unlink first list element and deallocate it
      first => head%next
      head%next => first%next
      call deallocate_node (first)
    end do
    head = default ! assign default initialization values to the table root node
  end subroutine mt_reset_table

  subroutine mt_dump_table (unit)

    use string_utilities, only: i_to_c

    integer, intent(in) :: unit

    type(node), pointer :: p

    write(unit,'(/,a)') 'TABLE HEAD: id=' // i_to_c(head%id) // &
                               ', name="' // trim(head%name) // '"'

    p => head%next
    if (associated(p)) then
      do while (associated(p))
        write(unit,'(2x,a)') 'matID=' // i_to_c(p%id) // ', name="' // trim(p%name) // '"'
        !! dump the contents of p%ms
        p => p%next
      end do
    else
      write(unit,'(2x,a)') '(none)'
    end if

  end subroutine mt_dump_table

 !!
 !! AUXILLARY ROUTINES
 !!
 !! SCAN_FOR_NAME/SCAN_FOR_ID searches TABLE for a node having the specified
 !!   name or ID and returns a pointer to the node if it is found; otherwise
 !!   a null pointer is returned.
 !!
 !! UNLINK searches TABLE for a node that is associated with the specified
 !!   node pointer P.  If it is found (there is at most one in a well-formed
 !!   table) it is unlinked from the list.  TABLE must be the head node of the
 !!   table.  Note that nothing is done with the pointer P itself; it is merely
 !!   unlinked from the list.
 !!
 !! DEALLOCATE_NODE deallocates all allocated storage associated with the data
 !!   components (principally the MAT_SYSTEM component) of the node pointer
 !!   THIS, and then deallocates the pointer itself.  If the pointer is
 !!   unassociated, the subroutine silently does nothing.
 !!

  function scan_for_name (head, name) result (p)
    type(node), intent(in) :: head
    character(len=*), intent(in) :: name
    type(node), pointer :: p
    p => head%next
    do while (associated(p))
     if (p%name == name) exit
     p => p%next
    end do
  end function scan_for_name

  function scan_for_id (head, id) result (p)
    type(node), intent(in) :: head
    integer, intent(in) :: id
    type(node), pointer :: p
    p => head%next
    do while (associated(p))
      if (p%id == id) exit
      p => p%next
    end do
  end function scan_for_id

  subroutine unlink (head, p)
    type(node), intent(inout) :: head
    type(node), pointer :: p
    type(node), pointer :: t  ! trailing pointer
    if (associated(head%next, p)) then
      head%next => p%next
    else
      t => head%next
      do while (associated(t%next))
        if (associated(t%next, p)) exit
        t => t%next
      end do
      if (associated(t%next)) t%next => p%next
    end if
  end subroutine unlink

  subroutine deallocate_node (this)
    type(node), pointer :: this
    if (associated(this)) then
      call destroy (this%ms)
      deallocate(this)
    end if
  end subroutine deallocate_node

end module material_table
