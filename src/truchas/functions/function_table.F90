!!
!! FUNCTION_TABLE
!!
!! Maintains a user-defined table of SCAFUN function objects indexed by name.
!!
!! PROGRAMMING INTERFACE
!!
!! The module defines a single instance of a function table that is held as
!! private module data.  All interaction with the table is performed using the
!! procedures described below.  The table is initially empty.  At any point in
!! time the table will contain a user-defined collection of SCAFUN function
!! objects that can be indexed by name -- a dictionary in other words.
!!
!! In the following procedures the NAME arguments are character variables with
!! values whose trimmed length does not exceed the value of the public parameter
!! FT_MAX_NAME_LEN.
!!
!!  FT_HAS_FUNCTION(NAME) returns true if the table contains the key NAME;
!!    otherwise it returns false.
!!
!!  FT_GET_FUNCTION(NAME) returns a pointer to the SCAFUN object associated
!!    with the key NAME.  If there is no such key in the table a null pointer
!!    is returned.
!!
!!  CALL FT_ADD_FUNCTION (NAME, FUNC) adds the SCAFUN object FUNC to the table
!!    and associates it to the key NAME.  It is an error if the table contains
!!    the key NAME.  Note that the table stores an independent copy of the
!!    actual argument FUNC, which, as a consequence may be destroyed or modified
!!    after making the call if desired without effecting the copy stored int
!!    the table.
!!
!!  CALL FT_DELETE_FUNCTION (NAME) deletes the function corresponding to the
!!    key NAME from the table.  If no such key exists, the table is unmodified.
!!    Note that any pointer associated with the deleted SCAFUN object is
!!    rendered invalid.
!!
!!  CALL FT_RESET_TABLE () deletes the all the functions from the table,
!!    returning it to its initial empty state.  Any pointers to the SCAFUN
!!    objects that were stored in the table are rendered invalid.
!!
!! IMPLEMENTATION NOTES
!!
!! The number of functions is expected to be limited, and the table is not
!! intended to be accessed intensively during the course of a calculation.
!! The intent is for application code to interrogate the table during the
!! initialization phase of a calculation and cache copies of the functions
!! (or perhaps pointers to them) as appropriate, and thereafter make few,
!! if any, accesses to the table.  Thus speed of the implementation is not
!! a concern; it is implemented as an unsorted linked-list of key-value
!! pairs -- the simplest possible implementation having no a priori size
!! limitation.
!!

#include "f90_assert.fpp"

module function_table

  use scalar_functions
  implicit none
  private
  
  public :: ft_has_function, ft_get_function
  public :: ft_add_function, ft_delete_function, ft_reset_table
  
  integer, parameter, public :: FT_MAX_NAME_LEN = 31

  type, private :: node
    character(len=FT_MAX_NAME_LEN) :: name = ''
    type(scafun) :: func
    type(node), pointer :: next => null()
  end type node
  type(node), save :: head

contains

  logical function ft_has_function (name)
    character(len=*), intent(in) :: name
    ft_has_function = associated(scan_for_name(head, name))
  end function ft_has_function
  
  function ft_get_function (name) result (func)
    character(len=*), intent(in) :: name
    type(scafun), pointer :: func
    type(node), pointer :: p
    func => null()
    p => scan_for_name(head, name)
    if (associated(p)) func => p%func
  end function ft_get_function
  
  subroutine ft_add_function (name, func)
  
    character(len=*), intent(in) :: name
    type(scafun), intent(in) :: func
    
    type(node), pointer :: p
    
    INSIST(len_trim(name) > 0 .and. len_trim(name) <= FT_MAX_NAME_LEN)
    INSIST(.not.associated(scan_for_name(head, name)))
    
    allocate(p)
    p%name = name
    p%func = func
    p%next => head%next
    head%next => p
    
  end subroutine ft_add_function
    
  subroutine ft_delete_function (name)
    character(len=*), intent(in) :: name
    type(node), pointer :: p
    p => scan_for_name(head, name)
    call unlink (head, p)
    call deallocate_node (p)
  end subroutine ft_delete_function

  subroutine ft_reset_table ()
    type(node), pointer :: first
    type(node) :: default
    do while (associated(head%next))
      first => head%next
      head%next => first%next
      call deallocate_node (first)
    end do
    head = default  ! assign default initialization values to the table head node
  end subroutine ft_reset_table
  
 !!
 !! AUXILLARY ROUTINES
 !!
 !! SCAN_FOR_NAME searches the table for a node having the specified name and
 !! returns a pointer to the node if it is found or a null pointer otherwise.
 !!
 !! UNLINK searches the table for a node that is associated with the specified
 !!   node pointer P.  If it is found (there is at most one in a well-formed
 !!   table) it is unlinked from the list.  HEAD must be the head node of the
 !!   table.  Note that nothing is done with the pointer P itself; it is merely
 !!   unlinked from the list.
 !!
 !! DEALLOCATE_NODE deallocates all allocated storage associated with the data
 !!   components (principally the SCAFUN component) of the node pointer THIS,
 !!   and then deallocates the pointer itself.  If the pointer is unassociated,
 !!   the subroutine silently does nothing.
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
      call destroy (this%func)
      deallocate(this)
    end if
  end subroutine deallocate_node

end module function_table
