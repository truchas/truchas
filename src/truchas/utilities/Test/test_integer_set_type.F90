!!
!! Neil N. Carlson <nnc@lanl.gov> March 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program test_integer_set_type

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use integer_set_type!, only: integer_set => iset
  implicit none

  integer :: stat = 0

  call test_add
  call test_remove
  call test_to_array
  call test_copy_to_array
  call test_elemental
  call test_add_set
  call test_iterator

  call exit (stat)

contains

  subroutine test_add

    type(integer_set) :: A

    if (.not.A%is_empty()) call write_fail ('test_add: set is not empty')
    if (A%size() /= 0) call write_fail ('test_add: size not 0')

    call A%add (3)  ! first element
    if (A%is_empty()) call write_fail ('test_add: set is empty')
    if (A%size() /= 1) call write_fail ('test_add: size not 1')

    call A%add (2)  ! insert at the beginning
    if (A%is_empty()) call write_fail ('test_add: set is empty')
    if (A%size() /= 2) call write_fail ('test_add: size not 2')

    call A%add (5)  ! insert at the end
    if (A%size() /= 3) call write_fail ('test_add: size not 3')

    call A%add (4)  ! insert in the middle
    if (A%size() /= 4) call write_fail ('test_add: size not 4')

    call A%add (4)  ! insert duplicate
    if (A%size() /= 4) call write_fail ('test_add: size not 4')

  end subroutine test_add

  subroutine test_remove

    type(integer_set) :: A

    call A%add (1)
    call A%add (2)
    call A%add (3)
    call A%add (4)
    call A%add (5)

    call A%remove (3)
    if (A%size() /= 4) call write_fail ('test_remove: middle, wrong size')

    call A%remove (3)
    if (A%size() /= 4) call write_fail ('test_remove: absent middle, wrong size')

    call A%remove (0)
    if (A%size() /= 4) call write_fail ('test_remove: absent before, wrong size')

    call A%remove (9)
    if (A%size() /= 4) call write_fail ('test_remove: absent after, wrong size')

    call A%remove (1)
    if (A%size() /= 3) call write_fail ('test_remove: first, wrong size')

    call A%remove (5)
    if (A%size() /= 2) call write_fail ('test_remove: last, wrong size')

    call A%remove (2)
    call A%remove (4)
    if (.not.A%is_empty()) call write_fail ('test_remove: not empty')

  end subroutine test_remove

  subroutine test_to_array

    type(integer_set) :: a
    integer, allocatable :: b(:)

    !! Assignment to unallocated array
    call a%add (1)
    call a%add (2)
    call a%add (3)
    b = a
    if (size(b) /= 3) then
      call write_fail ('test_to_array: unallocated, wrong size')
    else
      if (any(b /= [1, 2, 3])) call write_fail ('test_to_array: unallocated, wrong values')
    end if

    !! Assignment to allocated array with wrong size
    call a%add (0)
    b = a
    if (size(b) /= 4) then
      call write_fail ('test_to_array: reallocated, wrong size')
    else
      if (any(b /= [0, 1, 2, 3])) call write_fail ('test_to_array: reallocated, wrong values')
    end if

    !! Assignment to allocated array with right size
    b = -1  ! overwrite with wrong values
    b = a
    if (size(b) /= 4) then
      call write_fail ('test_to_array: allocated, wrong size')
    else
      if (any(b /= [0, 1, 2, 3])) call write_fail ('test_to_array: allocated, wrong values')
    end if

  end subroutine test_to_array


  subroutine test_copy_to_array

    type(integer_set) :: a
    integer :: b(3) = -1

    !! Assignment to unallocated array
    call a%add (1)
    call a%add (2)
    call a%copy_to_array (b)
    if (any(b /= [1, 2, -1])) call write_fail ('test_copy_to_array: wrong values')

  end subroutine test_copy_to_array


  subroutine test_elemental

    type(integer_set) :: A(2), B(1)
    integer, allocatable :: array(:), array2(:)

    if (any(.not.A%is_empty())) call write_fail ('test_elemental: not empty')
    call A%add (1)
    call A%add ([2,3])

    if (any(A%is_empty())) call write_fail ('test_elemental: empty')
    if (any(A%size() /= [2,2])) call write_fail ('test_elemental: wrong sizes')

    array = A(1)
    if (any(array /= [1,2])) call write_fail ('test_elemental: wrong values set 1')

    array = A(2)
    if (any(array /= [1,3])) call write_fail ('test_elemental: wrong values set 2')

    array2 = B(1)
    if (size(array2) /= 0) call write_fail ('test_elemental: wrong size for array2')

  end subroutine test_elemental

  subroutine test_add_set

    type(integer_set) :: a, b, c, d, e
    integer, allocatable :: array(:)

    call b%add (2)
    call b%add (4)
    call a%add (b)
    if (a%size() /= 2) then
      call write_fail ('test_set_add_set: empty, wrong size')
    else
      array = a
      if (any(array /= [2,4])) call write_fail ('test_set_add_set: empty, wrong values')
    end if

    call c%add (1)
    call c%add (2)
    call a%add (c)
    if (a%size() /= 3) then
      call write_fail ('test_set_add_set: before, dupl, wrong size')
    else
      array = a
      if (any(array /= [1,2,4])) call write_fail ('test_set_add_set: before, dupl, wrong values')
    end if

    call d%add (3)
    call d%add (4)
    call d%add (5)
    call a%add (d)
    if (a%size() /= 5) then
      call write_fail ('test_set_add_set: after, dupl, wrong size')
    else
      array = a
      if (any(array /= [1,2,3,4,5])) call write_fail ('test_set_add_set: after, dupl, wrong values')
    end if

    call a%add(b)
    if (a%size() /= 5) then
      call write_fail ('test_set_add_set: all dupl, wrong size')
    else
      array = a
      if (any(array /= [1,2,3,4,5])) call write_fail ('test_set_add_set: all dupl, wrong values')
    end if

    call a%add(e)
    if (a%size() /= 5) then
      call write_fail ('test_set_add_set: add empty, wrong size')
    else
      array = a
      if (any(array /= [1,2,3,4,5])) call write_fail ('test_set_add_set: add empty, wrong values')
    end if

  end subroutine test_add_set


  subroutine test_iterator

    type(integer_set) :: A
    type(integer_set_iterator) :: iter
    integer, parameter :: m = 4
    integer :: n, j
    logical :: tag(m) = .false.

    do j = 1, m
      call A%add (j)
    end do

    iter = A%begin()
    do j = 1, m
      if (iter%at_end()) then
        call write_fail ('test_iterator: iterator at end prematurely')
        exit
      end if
      n = iter%value()
      if (n < 1 .or. n > m) then
        call write_fail ('test_iterator: invalid value')
      else if (tag(n)) then
        call write_fail ('test_iterator: element already visited')
      else
        tag(n) = .true.
      end if
      call iter%next
    end do
    if (.not.iter%at_end()) call write_fail ('test_iterator: iterator not at end')

  end subroutine test_iterator


  subroutine write_fail (errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    stat = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_integer_set_type
