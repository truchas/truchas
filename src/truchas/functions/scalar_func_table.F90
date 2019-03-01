!!
!! SCALAR_FUNC_TABLE
!!
!! This module manages the global table of SCALAR_FUNC objects.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! This essentially defines a singleton of type SCALAR_FUNC_MAP.  The object
!! is private to the module and accessed via the following public procedures.
!!
!!  CALL INSERT_FUNC(NAME, F) adds the SCALAR_FUNC object F to the table and
!!    associates it with NAME.  If a function with this name already exists,
!!    its associated object is replaced by this one.  F is allocatable and its
!!    allocation is moved to the table and it is returned unallocated.
!!
!!  CALL LOOKUP_FUNC(NAME, F) returns the named function.  F is an allocatable
!!    class SCALAR_FUNC variable.  It is allocated by the subroutine and returns
!!    a copy (via sourced allocation) of the named function; if no such function
!!    exists, F is returned unallocated.
!!
!!  KNOWN_FUNC(NAME) returns the value .TRUE. if the table contains the function
!!    with the given NAME; otherwise it returns .FALSE.
!!

module scalar_func_table

  use scalar_func_class
  use scalar_func_map_type
  implicit none
  private

  public :: insert_func, lookup_func, known_func

  type(scalar_func_map) :: ftable

contains

  subroutine insert_func(name, f)
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(inout) :: f
    call ftable%insert(name, f)
  end subroutine insert_func

  subroutine lookup_func(name, f)
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(out) :: f
    call ftable%lookup (name, f)
  end subroutine lookup_func

  logical function known_func(name)
    character(*), intent(in) :: name
    known_func = ftable%mapped(name)
  end function known_func

end module scalar_func_table
