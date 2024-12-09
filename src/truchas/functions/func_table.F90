!!
!! FUNC_TABLE
!!
!! This module manages the global table of function objects.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2024; merged the SCALAR/VECTOR_FUNC_TABLE modules
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! This essentially defines a private singleton of type UMAP wrapped by the
!! following public procedures which limit and customize access to the various
!! function classes: SCALAR_FUNC and VECTOR_FUNC.
!!
!!  CALL INSERT_FUNC(NAME, F) adds the function object F to the table and
!!    associates it with NAME. If a function with this name already exists,
!!    its associated object is replaced by this one. F is allocatable and its
!!    allocation is moved to the table and it is returned unallocated.
!!
!!  CALL LOOKUP_FUNC(NAME, F) returns the named function. F is an allocatable
!!    function class variable. It is allocated by the subroutine and returns
!!    a copy (via sourced allocation) of the named function; if no such function
!!    of the class exists, F is returned unallocated.
!!
!!  KNOWN_SCALAR_FUNC(NAME)
!!  KNOWN_VECTOR_FUNC(NAME)
!!    Return the value .TRUE. if the table contains the function of the
!!    indicated class with the given NAME; otherwise it returns .FALSE.
!!

module func_table

  use umap_type
  use scalar_func_class
  use vector_func_class
  !use complex_scalar_func_class
  !use complex_vector_func_class
  implicit none
  private

  public :: insert_func, lookup_func
  public :: known_scalar_func, known_vector_func
  !public :: known_complex_scalar_func, known_complex_vector_func

  type(umap) :: ftable

  interface insert_func
    procedure insert_scalar_func, insert_vector_func
    !procedure insert_complex_scalar_func, insert_complex_vector_func
  end interface

  interface lookup_func
    procedure lookup_scalar_func, lookup_vector_func
    !procedure lookup_complex_scalar_func, lookup_complex_vector_func
  end interface

contains

  subroutine insert_scalar_func(name, f)
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(inout) :: f
    class(*), allocatable :: tmp
    call move_alloc(f, tmp)
    call ftable%insert(name, tmp) ! tmp returned unallocated
  end subroutine

  subroutine insert_vector_func(name, f)
    character(*), intent(in) :: name
    class(vector_func), allocatable, intent(inout) :: f
    class(*), allocatable :: tmp
    call move_alloc(f, tmp)
    call ftable%insert(name, tmp) ! tmp returned unallocated
  end subroutine

!  subroutine insert_complex_scalar_func(name, f)
!    character(*), intent(in) :: name
!    class(complex_scalar_func), allocatable, intent(inout) :: f
!    class(*), allocatable :: tmp
!    call move_alloc(f, tmp)
!    call ftable%insert(name, tmp) ! tmp returned unallocated
!  end subroutine
!
!  subroutine insert_complex_vector_func(name, f)
!    character(*), intent(in) :: name
!    class(complex_vector_func), allocatable, intent(inout) :: f
!    class(*), allocatable :: tmp
!    call move_alloc(f, tmp)
!    call ftable%insert(name, tmp) ! tmp returned unallocated
!  end subroutine

  subroutine lookup_scalar_func(name, f)
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(out) :: f
    class(*), pointer :: uptr
    call ftable%lookup(name, uptr)
    if (.not.associated(uptr)) return
    select type (uptr)
    class is (scalar_func)
      allocate(f, source=uptr)
    end select
  end subroutine

  subroutine lookup_vector_func(name, f)
    character(*), intent(in) :: name
    class(vector_func), allocatable, intent(out) :: f
    class(*), pointer :: uptr
    call ftable%lookup(name, uptr)
    if (.not.associated(uptr)) return
    select type (uptr)
    class is (vector_func)
      allocate(f, source=uptr)
    end select
  end subroutine

!  subroutine lookup_complex_scalar_func(name, f)
!    character(*), intent(in) :: name
!    class(complex_scalar_func), allocatable, intent(out) :: f
!    class(*), pointer :: uptr
!    call ftable%lookup(name, uptr)
!    if (.not.associated(uptr)) return
!    select type (uptr)
!    class is (complex_scalar_func)
!      allocate(f, source=uptr)
!    end select
!  end subroutine
!
!  subroutine lookup_vector_func(name, f)
!    character(*), intent(in) :: name
!    class(complex_vector_func), allocatable, intent(out) :: f
!    class(*), pointer :: uptr
!    call ftable%lookup(name, uptr)
!    if (.not.associated(uptr)) return
!    select type (uptr)
!    class is (complex_vector_func)
!      allocate(f, source=uptr)
!    end select
!  end subroutine

  logical function known_scalar_func(name) result(known)
    character(*), intent(in) :: name
    class(*), pointer :: uptr
    call ftable%lookup(name, uptr)
    known = associated(uptr)
    if (known) then
      select type (uptr)
      class is (scalar_func)
      class default
        known = .false.
      end select
    end if
  end function

  logical function known_vector_func(name) result(known)
    character(*), intent(in) :: name
    class(*), pointer :: uptr
    call ftable%lookup(name, uptr)
    known = associated(uptr)
    if (known) then
      select type (uptr)
      class is (vector_func)
      class default
        known = .false.
      end select
    end if
  end function

!  logical function known_complex_scalar_func(name) result(known)
!    character(*), intent(in) :: name
!    class(*), pointer :: uptr
!    call ftable%lookup(name, uptr)
!    known = associated(uptr)
!    if (known) then
!      select type (uptr)
!      class is (complex_scalar_func)
!      class default
!        known = .false.
!      end select
!    end if
!  end function
!
!  logical function known_complex_vector_func(name) result(known)
!    character(*), intent(in) :: name
!    class(*), pointer :: uptr
!    call ftable%lookup(name, uptr)
!    known = associated(uptr)
!    if (known) then
!      select type (uptr)
!      class is (complex_vector_func)
!      class default
!        known = .false.
!      end select
!    end if
!  end function

end module func_table
