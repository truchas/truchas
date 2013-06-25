!!
!! PHASE_PROPERTY_TABLE
!!
!! Maintains user-defined collections of phases and properties, and the
!! function objects that are assigned to certain phase-property pairs.
!!
!! Neil Carlson <nnc@lanl.gov>
!! 13 Nov 2008
!!
!! PROGRAMMING INTERFACE
!!
!! The module defines a single instance of a phase property table that is held
!! as private module data.  All interaction with the table is performed using
!! the procedures described below.  The table is initially empty.  At any point
!! in time the table will contain user-defined collections of phase names and
!! property names.  Associated with each phase or property name is a unique
!! positive integer ID assigned by the system.  Many of the procedures take the
!! IDs associated with the phase and property names as arguments rather than
!! the names themselves.  Each phase-property pair may have a scalar-valued
!! function associated with it, implemented as a SCAFUN object; this function
!! provides the value of the property for the phase.
!!
!! In the following procedures, IDs are integer variables/values and names are
!! character variables/values whose length does not exceed the value of the
!! public parameter PPT_MAX_NAME_LEN.
!!
!!  PPT_PHASE_ID(NAME) returns the ID assigned to the specified phase name.
!!    It is an error if NAME is not a defined phase name.
!!
!!  PPT_PROPERTY_ID(NAME) returns the ID assigned to the specified property
!!    name.  It is an error if NAME is not a defined property name.
!!
!!  PPT_PHASE_NAME(ID) returns the phase name corresponding to ID.  It is an
!!    error if ID is not assigned to a phase name.
!!
!!  PPT_PROPERTY_NAME(ID) returns the property name corresponding to ID.
!!    It is an error if ID is not assigned to a property name.
!!
!!  PPT_HAS_PHASE(NAME) returns true if NAME is a defined phase name; otherwise
!!    it returns false.
!!
!!  PPT_HAS_PROPERTY(NAME) returns true if NAME is a defined property name;
!!    otherwise it returns false.
!!
!!  PPT_VALID_PHASE(ID) returns true if ID is a valid phase ID; otherwise
!!    it returns false.
!!
!!  PPT_VALID_PROPERTY(ID) returns true if ID is a valid property ID;
!!    otherwise it returns false.
!!
!!  PPT_HAS_PHASE_PROPERTY(PHASE_ID, PROPERTY_ID) returns true if a function
!!    is assigned to the phase-property pair with the specified IDs;
!!    otherwise it returns false.  It is an error if either of the IDs is
!!    invalid.
!!
!!  PPT_NUM_PHASE() returns the number of defined phases.
!!
!!  PPT_NUM_PROPERTY() returns the number of defined properties.
!!
!!  PPT_GET_PHASE_IDS(ID) returns an array of all the phase IDs in the rank-1
!!    integer pointer array ID.  The array is allocated by the subroutine with
!!    a size equal to the number of phases.  The caller is responsible to
!!    deallocating the array.
!!
!!  PPT_GET_PROPERTY_IDS(ID) returns an array of all the property IDs in the
!!    rank-1 integer pointer array ID. The array is allocated by the subroutine
!!    with a size equal to the number of phases.  The caller is responsible to
!!    deallocating the array.
!!
!!  CALL PPT_GET_PHASE_NAMES (NAME) returns an array of all the phase names.
!!    NAME is a pointer to a rank-1 character array that is allocated by
!!    the subroutine with a size equal to the number of phases.  The character
!!    length must be at least large enough to hold the largest trimmed length
!!    of the names, which is no greater than PPT_MAX_NAME_LEN.
!!
!!  CALL PPT_GET_PROPERTY_NAMES (NAME) returns an array of all the property
!!    names.  NAME is a pointer to a rank-1 character array that is allocated
!!    by the subroutine with a size equal to the number of properties.  The
!!    character length must be at least large enough to hold the largest
!!    trimmed length of the names, which is no greater than PPT_MAX_NAME_LEN.
!!
!!  CALL PPT_GET_PHASE_PROPERTY (PHASE_ID, PROPERTY_ID, F) returns a pointer F
!!    to the SCAFUN object assigned to the phase-property pair with the
!!    specified IDs, or a null pointer if no function is assigned.  It is an
!!    error if either of the IDs is invalid.
!!
!!  CALL PPT_ADD_PHASE (NAME, ID) adds a phase with the specified NAME to the
!!    table and returns the assigned ID.  If the phase name already exists its
!!    assigned ID is simply returned and the table is not modified.
!!
!!  CALL PPT_ADD_PROPERTY (NAME, ID) adds a property with the specified NAME to
!!    the table and returns the assigned ID.  If the property name already
!!    exists its assigned ID is simply returned and the table is not modified.
!!
!!  CALL PPT_ASSIGN_PHASE_PROPERTY (PHASE_ID, PROPERTY_ID, F) assigns the
!!    SCAFUN object F to the phase-property pair with the specified IDs.
!!    It is an error if either of the IDs is invalid or if a function has
!!    already been assigned to the pair.  In the latter case use
!!    PPT_REASSIGN_PHASE_PROPERTY instead.  Note that the table stores a
!!    completely independent copy of the actual argument F, which, as a
!!    consequence, may be destroyed after making the call if desired.
!!
!!  CALL PPT_REASSIGN_PHASE_PROPERTY (PHASE_ID, PROPERTY_ID, F) reassigns the
!!    SCAFUN object F to the phase-property pair with the specified IDs.
!!    It is an error if either of the IDs is invalid or if a function is not
!!    already assigned to the pair.  Note that the original function is
!!    destroyed, rendering any pointers to it, such as were returned by
!!    PPT_GET_PHASE_PROPERTY, invalid.
!!
!!  CALL PPT_DELETE_PHASE (ID) deletes the phase with the specified ID,
!!    including all of its property functions, from the table.  It is an error
!!    if ID is not a valid phase ID.  Any pointers to the deleted SCAFUN
!!    function objects are rendered invalid.
!!
!!  CALL PPT_DELETE_PROPERTY (ID) deletes the property with the specified ID,
!!    including its functions for all the phases, from the table.  It is an
!!    error if ID is not a valid property ID.  Any pointers to the deleted
!!    SCAFUN function objects are rendered invalid.
!!
!!  CALL PPT_DELETE_PHASE_PROPERTY (PHASE_ID, PROPERTY_ID) deletes the SCAFUN
!!    object assigned to the phase-property pair with the specified IDs
!!    from the table.  It is an error if either of the IDs is invalid.  Any
!!    pointer to the deleted SCAFUN function object is rendered invalid.
!!
!!  CALL PPT_RESET_TABLE () deletes all the phases, properties, and assigned
!!    functions from the phase property table, returning it to its initial
!!    empty state.  Any pointers to the SCAFUN function objects that were
!!    stored in the table are rendered invalid.
!!
!!  CALL PPT_DUMP_TABLE (UNIT) writes the contents of the table to the
!!    specified logical unit in a moderately readable form.  Intended for
!!    debugging use.
!!
!! IMPLEMENTATION NOTES
!!
!! The phases and properties are expected to number a few dozen each at most,
!! and the table is not intended to be accessed intensively during the course
!! of a calculation.  The intent is for application code to cache pointers to
!! the SCAFUN function objects returned by PPT_GET_PHASE_PROPERTY during an
!! initialization phase of a calculation and thereafter make few, if any,
!! accesses to the table.  Thus speed of the implementation is not a concern,
!! and so, for example, we do not bother with such things as quickly searchable
!! lists; all lists are unsorted with insertion at the head -- the simplest
!! possible implementation (with no preset size limitation).
!!
!! With that as the context, the table is structured as a 'matrix' of doubly-
!! linked nodes -- across rows and down columns.  A node consists of data
!! components -- row and column IDs, name string, and a SCAFUN object -- plus
!! two node pointers -- one to the next node in the row and the other to the
!! next node in the column that the node belongs to.  Thus we obtain two
!! families of linear linked lists, 'rows' and 'columns', depending on
!! whether we are following the row or column pointers.  The nodes in a row
!! all have the same row ID data value, and the nodes in a column all have
!! the same column ID data value.  The entry point into the table is the
!! table head node.  The column linked-list that it heads forms the list
!! of phases: name stores the phase name, and the row ID the assigned ID.
!! Each of the nodes in this initial column (except for the starting table
!! head node) also serves as the head node for a linked-list row of the table
!! by following the row pointer.  Similarly, the row linked-list headed by
!! the table head node forms the list of properties: name stores the property
!! name and the column ID the assigned ID.  Each of the nodes in this initial
!! row (except for the starting table head node) also serves as the head node
!! of a linked-list column of the table by following the column pointer.
!! Proper table nodes -- that is nodes other than the table head node and
!! the row/column header nodes just described -- correspond to phase-property
!! pairs having an assigned function.  Each of these nodes store the phase
!! and property IDs as the row and column IDs and the assigned function.
!! Moreover they appear in the row linked list headed by the phase node
!! with the same row ID, and the column linked list headed by the property
!! node with the same column ID.
!!
!! The data components of a node are used (or not used) depending on how
!! the node is used.  The name component is currently unused by proper table
!! nodes, though this could be assigned the 'name' of the function if that
!! information was available.  Row header nodes also use the column ID to store
!! the length of the row linked list; that would be the number of assigned
!! property functions for a particular phase.  Similarly column header nodes
!! use the row ID to store the length of the column linked list; that would
!! be the number of phases that have a function assigned to a particular
!! property.  None of the row, column, and table header nodes use the SCAFUN
!! object component.  The table head node does not use the name, though if
!! we were to move to multiple table instances this could hold the name of
!! the table.  The row ID component of the table head node is used to dole
!! out the IDs assigned to phases or properties, being incremented each time
!! a phase or property is added to the table.  (Note that we could do
!! separate ID pools for phases and properties by using both the row and
!! column ID components of the table head node.)  It was tempting to use
!! the row/column ID components of the table head node to store the number
!! of phases and properties, but this would force us to use an external
!! method of assigning IDS -- the current implementation is self-contained
!! within the node structure.
!!

#include "f90_assert.fpp"

module phase_property_table

  use scalar_functions
  implicit none
  private

  integer, parameter, public :: PPT_MAX_NAME_LEN = 31

  public :: ppt_phase_id, ppt_property_id, ppt_phase_name, ppt_property_name
  public :: ppt_has_phase, ppt_has_property, ppt_valid_phase, ppt_valid_property
  public :: ppt_has_phase_property, ppt_get_phase_property
  public :: ppt_add_phase, ppt_add_property
  public :: ppt_assign_phase_property, ppt_reassign_phase_property
  public :: ppt_delete_phase, ppt_delete_property, ppt_delete_phase_property, ppt_reset_table
  public :: ppt_num_phase, ppt_num_property, ppt_get_phase_names, ppt_get_property_names
  public :: ppt_get_phase_ids, ppt_get_property_ids
  public :: ppt_dump_table

  !! Single instance of the table held as private module data;
  !! All procedures operate on this data.
  type, private :: node
    integer :: row_id = 0
    integer :: col_id = 0
    character(len=PPT_MAX_NAME_LEN) :: name = ''
    type(scafun) :: f
    type(node), pointer :: row_next => null()
    type(node), pointer :: col_next => null()
  end type
  type(node), save :: table

contains

  integer function ppt_phase_id (name)
    character(len=*), intent(in) :: name
    type(node), pointer :: p
    p => scan_col_for_name(table, name)
    INSIST(associated(p))
    ppt_phase_id = p%row_id
  end function ppt_phase_id

  integer function ppt_property_id (name)
    character(len=*), intent(in) :: name
    type(node), pointer :: p
    p => scan_row_for_name(table, name)
    INSIST(associated(p))
    ppt_property_id = p%col_id
  end function ppt_property_id

  function ppt_phase_name (id) result (name)
    integer, intent(in) :: id
    character(len=PPT_MAX_NAME_LEN) :: name
    type(node), pointer :: p
    p => scan_col_for_row_id(table, id)
    INSIST(associated(p))
    name = p%name
  end function ppt_phase_name

  function ppt_property_name (id) result (name)
    integer, intent(in) :: id
    character(len=PPT_MAX_NAME_LEN) :: name
    type(node), pointer :: p
    p => scan_row_for_col_id(table, id)
    INSIST(associated(p))
    name = p%name
  end function ppt_property_name

  logical function ppt_has_phase (name)
    character(len=*), intent(in) :: name
    ppt_has_phase = associated(scan_col_for_name(table, name))
  end function ppt_has_phase

  logical function ppt_has_property (name)
    character(len=*), intent(in) :: name
    ppt_has_property = associated(scan_row_for_name(table, name))
  end function ppt_has_property

  logical function ppt_valid_phase (id)
    integer, intent(in) :: id
    ppt_valid_phase = valid_phase_id(table, id)
  end function ppt_valid_phase

  logical function ppt_valid_property (id)
    integer, intent(in) :: id
    ppt_valid_property = valid_property_id(table, id)
  end function ppt_valid_property

  subroutine ppt_add_phase (name, id) !, stat, errmsg)

    character(len=*), intent(in) :: name
    integer, intent(out) :: id
    !integer, intent(out), optional :: stat
    !character(len=*), intent(out), optional :: errmsg

    type(node), pointer :: p

    INSIST( len_trim(name) > 0 .and. len_trim(name) <= PPT_MAX_NAME_LEN )

    p => scan_col_for_name(table, name)

    if (.not.associated(p)) then
      !! Insert a row header for the new phase.
      allocate(p)
      table%row_id = table%row_id + 1 ! generate a new ID
      p%row_id = table%row_id
      p%name = name
      p%col_next => table%col_next
      table%col_next => p
    end if

    id = p%row_id

  end subroutine ppt_add_phase

  subroutine ppt_add_property (name, id) !, stat, errmsg)

    character(len=*), intent(in) :: name
    integer, intent(out) :: id
    !integer, intent(out), optional :: stat
    !character(len=*), intent(out), optional :: errmsg

    type(node), pointer :: p

    INSIST( len_trim(name) > 0 .and. len_trim(name) <= PPT_MAX_NAME_LEN )

    p => scan_row_for_name(table, name)

    if (.not.associated(p)) then
      !! Insert a column header for the new property.
      allocate(p)
      table%row_id = table%row_id + 1 ! generate new ID
      p%col_id = table%row_id
      p%name = name
      p%row_next => table%row_next
      table%row_next => p
    end if

    id = p%col_id

  end subroutine ppt_add_property

  logical function ppt_has_phase_property (phase_id, property_id)
    integer, intent(in) :: phase_id, property_id
    INSIST(valid_phase_id(table,phase_id))
    INSIST(valid_property_id(table,property_id))
    ppt_has_phase_property = associated(scan_table(table, phase_id, property_id))
  end function ppt_has_phase_property

  subroutine ppt_assign_phase_property (phase_id, property_id, f) !, stat, errmsg)

    integer, intent(in) :: phase_id, property_id
    type(scafun), intent(in) :: f
    !integer, intent(out), optional :: stat
    !character(len=*), intent(out), optional errmsg

    type(node), pointer :: p, new

    INSIST(valid_phase_id(table,phase_id))
    INSIST(valid_property_id(table,property_id))
    INSIST(.not.associated(scan_table(table,phase_id,property_id)))

    !! Create and initialize a new phase property node for the table.
    allocate(new)
    new%row_id = phase_id
    new%col_id = property_id
    new%f = f

    !! Insert the new node into the row list.
    p => scan_col_for_row_id(table, phase_id)
    ASSERT(associated(p))
    new%row_next => p%row_next
    p%row_next => new
    p%col_id = p%col_id + 1 ! update the row list length

    !! Insert the new node into the column list.
    p => scan_row_for_col_id(table, property_id)
    ASSERT(associated(p))
    new%col_next => p%col_next
    p%col_next => new
    p%row_id = p%row_id + 1 ! update the column list length

  end subroutine ppt_assign_phase_property

  subroutine ppt_reassign_phase_property (phase_id, property_id, f) !, stat, errmsg)

    integer, intent(in) :: phase_id, property_id
    type(scafun), intent(in) :: f
    !integer, intent(out), optional :: stat
    !character(len=*), intent(out), optional errmsg

    type(node), pointer :: p

    p => scan_table(table, phase_id, property_id)
    INSIST(associated(p))

    !! Reassign the function.
    call destroy (p%f)
    p%f = f

  end subroutine ppt_reassign_phase_property

  subroutine ppt_get_phase_property (phase_id, property_id, f)  !, stat, errmsg)

    integer, intent(in) :: phase_id, property_id
    type(scafun), pointer :: f
    !integer, intent(out), optional :: stat
    !character(len=*), intent(out), optional errmsg

    type(node), pointer :: p

    INSIST(valid_phase_id(table,phase_id))
    INSIST(valid_property_id(table,property_id))

    p => scan_table(table, phase_id, property_id)
    if (associated(p)) then
      f => p%f
    else
      f => null()
    end if

  end subroutine ppt_get_phase_property

  subroutine ppt_delete_phase (id)
    integer, intent(in) :: id
    INSIST(valid_phase_id(table,id))
    call delete_row (table, id)
  end subroutine ppt_delete_phase

  subroutine ppt_delete_property (id)
    integer, intent(in) :: id
    INSIST(valid_property_id(table,id))
    call delete_column (table, id)
  end subroutine ppt_delete_property

  subroutine delete_row (table, id)

    type(node), intent(inout), target :: table
    integer, intent(in) :: id

    type(node), pointer :: row_head, col_head, p

    row_head => scan_col_for_row_id(table, id)
    if (.not.associated(row_head)) return

    !! Unlink all the row entries from their respective columns and
    !! then deallocate the row list entries.
    p => row_head%row_next
    do while (associated(p))
      col_head => scan_row_for_col_id(table, p%col_id)
      ASSERT(associated(col_head))
      call unlink_from_column (col_head, p)
      col_head%row_id = col_head%row_id - 1  ! update column list length
      p => p%row_next
    end do
    call deallocate_row (row_head)

    !! Finally unlink the row header and deallocate it.
    call unlink_from_column (table, row_head)
    call deallocate_node (row_head)

  end subroutine delete_row

  subroutine delete_column (table, id)

    type(node), intent(inout), target :: table
    integer, intent(in) :: id

    type(node), pointer :: row_head, col_head, p

    col_head => scan_row_for_col_id(table, id)
    if (.not.associated(col_head)) return

    !! Unlink all the column entries from their respective rows and
    !! then deallocate the column list entries.
    p => col_head%col_next
    do while (associated(p))
      row_head => scan_col_for_row_id(table, p%row_id)
      ASSERT(associated(row_head))
      call unlink_from_row (row_head, p)
      row_head%col_id = row_head%col_id - 1  ! update row list length
      p => p%col_next
    end do
    call deallocate_column (col_head)

    !! Finally unlink the column header and deallocate it.
    call unlink_from_row (table, col_head)
    call deallocate_node (col_head)

  end subroutine delete_column

  subroutine ppt_delete_phase_property (phase_id, property_id)
    integer, intent(in) :: phase_id, property_id
    INSIST(valid_phase_id(table,phase_id))
    INSIST(valid_property_id(table,property_id))
    call delete_entry (table, phase_id, property_id)
  end subroutine ppt_delete_phase_property

  subroutine delete_entry (table, row_id, col_id)

    type(node), intent(inout) :: table
    integer, intent(in) :: row_id, col_id

    type(node), pointer :: row_head, col_head, p

    row_head => scan_col_for_row_id(table, row_id)
    if (.not.associated(row_head)) return

    col_head => scan_row_for_col_id(table, col_id)
    if (.not.associated(col_head)) return

    p => scan_row_for_col_id(row_head, col_id)
    if (.not.associated(p)) return

    call unlink_from_row (row_head, p)
    row_head%col_id = row_head%col_id - 1 ! update row list length

    call unlink_from_column (col_head, p)
    col_head%row_id = col_head%row_id - 1 ! update column list length

    call deallocate_node (p)

  end subroutine delete_entry

  subroutine ppt_reset_table ()

    type(node), pointer :: first
    type(node) :: default

    do while (associated(table%col_next))
      first => table%col_next
      table%col_next => first%col_next  ! drop first row from the table
      call deallocate_row (first)
      call deallocate_node (first)
    end do

    call deallocate_row (table)

    table = default ! assign default initialization values to the table root node

  end subroutine ppt_reset_table

  integer function ppt_num_phase ()
    type(node), pointer :: p
    ppt_num_phase = 0
    p => table%col_next
    do while (associated(p))
      ppt_num_phase = ppt_num_phase + 1
      p => p%col_next
    end do
  end function ppt_num_phase

  integer function ppt_num_property ()
    type(node), pointer :: p
    ppt_num_property = 0
    p => table%row_next
    do while (associated(p))
      ppt_num_property = ppt_num_property + 1
      p => p%row_next
    end do
  end function ppt_num_property

  subroutine ppt_get_phase_names (name)
    character(len=*), pointer :: name(:)
    integer :: j
    type(node), pointer :: p
    allocate(name(ppt_num_phase()))
    p => table%col_next
    do j = size(name), 1, -1
      ASSERT(associated(p))
      INSIST(len_trim(p%name) <= len(name))
      name(j) = p%name
      p => p%col_next
    end do
  end subroutine ppt_get_phase_names

  subroutine ppt_get_property_names (name)
    character(len=*), pointer :: name(:)
    integer :: j
    type(node), pointer :: p
    allocate(name(ppt_num_property()))
    p => table%row_next
    do j = size(name), 1, -1
      ASSERT(associated(p))
      INSIST(len_trim(p%name) <= len(name))
      name(j) = p%name
      p => p%row_next
    end do
  end subroutine ppt_get_property_names

  subroutine ppt_get_phase_ids (id)
    integer, pointer :: id(:)
    integer :: j
    type(node), pointer :: p
    allocate(id(ppt_num_phase()))
    p => table%col_next
    do j = size(id), 1, -1
      ASSERT(associated(p))
      id(j) = p%row_id
      p => p%col_next
    end do
  end subroutine ppt_get_phase_ids

  subroutine ppt_get_property_ids (id)
    integer, pointer :: id(:)
    integer :: j
    type(node), pointer :: p
    allocate(id(ppt_num_property()))
    p => table%row_next
    do j = size(id), 1, -1
      ASSERT(associated(p))
      id(j) = p%col_id
      p => p%row_next
    end do
  end subroutine ppt_get_property_ids

  subroutine ppt_dump_table (unit)

    use string_utilities, only: i_to_c

    integer, intent(in) :: unit

    type(node), pointer :: p, q

    write(unit,'(/,a)') 'TABLE HEAD: row_id=' // i_to_c(table%row_id) // &
                                  ', col_id=' // i_to_c(table%col_id) // &
                                  ', name="'  // trim(table%name) // '"'

    !! Dump the rows of the table
    write(unit,'(a)') 'PHASES:'
    p => table%col_next
    if (associated(p)) then
      do while (associated(p))
        write(unit,'(2x,a)') 'phaseID=' // i_to_c(p%row_id) // &
                             ', name="' // trim(p%name) // &
                             '", functions=' // i_to_c(p%col_id)
        write(unit,'(4x,a)') 'FUNCTIONS:'
        q => p%row_next
        if (associated(q)) then
        do while (associated(q))
          write(unit,'(6x,a)') 'phaseID=' // i_to_c(q%row_id) // &
                               ', propertyID=' // i_to_c(q%col_id) // &
                               ', name="' // trim(q%name) // '"'
          q => q%row_next
        end do
        else
          write(unit,'(6x,a)') '(none)'
        end if
        p => p%col_next
      end do
    else
      write(unit, '(2x,a)') '(none)'
    end if

    !! Dump the columns of the table.
    write(unit,'(a)') 'PROPERTIES:'
    p => table%row_next
    if (associated(p)) then
      do while (associated(p))
        write(unit,'(2x,a)') 'propertyID=' // i_to_c(p%col_id) // &
                             ', name="' // trim(p%name) // &
                             '", functions=' // i_to_c(p%row_id)
        write(unit,'(4x,a)') 'FUNCTIONS:'
        q => p%col_next
        if (associated(q)) then
        do while (associated(q))
          write(unit,'(6x,a)') 'phaseID=' // i_to_c(q%row_id) // &
                               ', propertyID=' // i_to_c(q%col_id) // &
                               ', name="' // trim(q%name) // '"'
          q => q%col_next
        end do
        else
          write(unit,'(6x,a)') '(none)'
        end if
        p => p%row_next
      end do
    else
      write(unit, '(2x,a)') '(none)'
    end if

  end subroutine ppt_dump_table

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SCAN_COL_FOR_NAME   / SCAN_ROW_FOR_NAME
 !! SCAN_COL_FOR_ROW_ID / SCAN_ROW_FOR_COL_ID
 !!
 !! These auxillary functions search a row or column for a node having the
 !! specified name or row/column ID and return a pointer to the node if it is
 !! found; otherwise a null pointer is returned.  HEAD is the head node of the
 !! row/column. Note that if the table head is passed then it is the row/column
 !! headers themselves that are being searched.
 !!

  function scan_col_for_name (head, name) result (p)
    type(node), intent(in) :: head
    character(len=*), intent(in) :: name
    type(node), pointer :: p
    p => head%col_next
    do while (associated(p))
      if (p%name == name) exit
      p => p%col_next
    end do
  end function scan_col_for_name

  function scan_row_for_name (head, name) result (p)
    type(node), intent(in) :: head
    character(len=*), intent(in) :: name
    type(node), pointer :: p
    p => head%row_next
    do while (associated(p))
      if (p%name == name) exit
      p => p%row_next
    end do
  end function scan_row_for_name

  function scan_col_for_row_id (head, id) result (p)
    type(node), intent(in) :: head
    integer, intent(in) :: id
    type(node), pointer :: p
    p => head%col_next
    do while (associated(p))
      if (p%row_id == id) exit
      p => p%col_next
    end do
  end function scan_col_for_row_id

  function scan_row_for_col_id (head, id) result (p)
    type(node), intent(in) :: head
    integer, intent(in) :: id
    type(node), pointer :: p
    p => head%row_next
    do while (associated(p))
      if (p%col_id == id) exit
      p => p%row_next
    end do
  end function scan_row_for_col_id

  logical function valid_phase_id (table, id)
    type(node), intent(in) :: table
    integer, intent(in) :: id
    valid_phase_id = associated(scan_col_for_row_id(table, id))
  end function valid_phase_id

  logical function valid_property_id (table, id)
    type(node), intent(in) :: table
    integer, intent(in) :: id
    valid_property_id = associated(scan_row_for_col_id(table, id))
  end function valid_property_id

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SCAN_TABLE
 !!
 !! Searches the table for the entry having the specified ROW_ID and COL_ID,
 !! and returns a pointer to the table node if it is found; otherwise a null
 !! pointer is returned.  TABLE must be the table header node.
 !!

  function scan_table (table, row_id, col_id) result (p)
    type(node), intent(in) :: table
    integer,    intent(in) :: row_id, col_id
    type(node), pointer    :: p
    p => scan_col_for_row_id(table, row_id)
    if (associated(p)) p => scan_row_for_col_id(p, col_id)
  end function scan_table

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! UNLINK_FROM_ROW / UNLINK_FROM_COLUMN
 !!
 !! These auxillary subroutines search a row/column for a node that is
 !! associated with the specified node pointer P.  If it is found (there is
 !! at most one in a well-formed table) it is unlinked from the row/column
 !! list.  HEAD is the head node of the row/column list.  Note that nothing
 !! is done with the pointer P itself; it is merely unlinked from the list.
 !!
 !! Implementation note: because HEAD is itself not a pointer, the initial
 !! item in the list must be treated specially; otherwise the general code
 !! would suffice for all cases after making the necessary adjustments to
 !! the initialization of the trailing pointer.  Why not declare HEAD to be
 !! a pointer?  Actually all the row/column head nodes *are* pointers, except
 !! for the table head node itself.  I just thought it was more transparent
 !! (at the higher level) to not make it a pointer, which would require
 !! allocation at some point.
 !!

  subroutine unlink_from_row (head, p)

    type(node), intent(inout) :: head
    type(node), pointer :: p

    type(node), pointer :: t  ! trailing pointer

    if (associated(head%row_next, p)) then
      head%row_next => p%row_next
    else
      t => head%row_next
      do while(associated(t%row_next))
        if (associated(t%row_next, p)) exit
        t => t%row_next
      end do
      if (associated(t%row_next)) t%row_next => p%row_next
    end if

  end subroutine unlink_from_row

  subroutine unlink_from_column (head, p)

    type(node), intent(inout) :: head
    type(node), pointer :: p

    type(node), pointer :: t  ! trailing pointer

    if (associated(head%col_next, p)) then
      head%col_next => p%col_next
    else
      t => head%col_next
      do while(associated(t%col_next))
        if (associated(t%col_next, p)) exit
        t => t%col_next
      end do
      if (associated(t%col_next)) t%col_next => p%col_next
    end if

  end subroutine unlink_from_column

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEALLOCATE_NODE
 !!
 !! This auxillary routine deallocates all allocated storage associated with
 !! the data components of the node pointer THIS, and then deallocates the
 !! pointer itself.  If the pointer is unassociated, the subroutine silently
 !! does nothing.  The main intent here is to destroy the function component.
 !!

  subroutine deallocate_node (this)
    type(node), pointer :: this
    if (associated(this)) then
      call destroy (this%f)
      deallocate(this)
    end if
  end subroutine deallocate_node

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEALLOCATE_ROW / DEALLOCATE_COLUMN
 !!
 !! These auxillary subroutines deallocate a row/column linked list and all
 !! embedded allocated storage associated with the nodes in the list.  HEAD
 !! is the head node of the row/column.  Note that the head node itself is
 !! not deallocated.
 !!
 !! N.B.  When deallocating a row, for example, the column links of the nodes
 !! are completely ignored, and thus if any of the nodes appear in the column
 !! linked lists the table will be corrupted.  The UNLINK_FROM_COLUMN routine
 !! should be used beforehand, unless the entire table is being deleted.
 !!

  subroutine deallocate_row (head)
    type(node), intent(inout) :: head
    type(node), pointer :: first
    do while (associated(head%row_next))
      first => head%row_next
      head%row_next => first%row_next ! drop first item from list
      call deallocate_node (first)
    end do
  end subroutine deallocate_row

  subroutine deallocate_column (head)
    type(node), intent(inout) :: head
    type(node), pointer :: first
    do while (associated(head%col_next))
      first => head%col_next
      head%col_next => first%col_next ! drop first item from list
      call deallocate_node (first)
    end do
  end subroutine deallocate_column

end module phase_property_table
