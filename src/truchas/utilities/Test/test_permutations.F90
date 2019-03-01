!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_permutations

  use permutations
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif

  integer :: stat = 0
  integer, allocatable :: p(:)
  integer :: p1(7) = [2, 5, 7, 1, 3, 4, 6] ! 1 cycle
  integer :: p2(7) = [3, 6, 5, 2, 7, 4, 1] ! 2 cycle
  integer :: p3(7) = [3, 2, 6, 5, 4, 1, 7] ! 2 cycle w/ fixed pt and swap

  real, allocatable :: array_in(:,:), array_out(:,:), array_tmp(:,:)
  integer, allocatable :: xrag1(:), rag1(:), xrag2(:), rag2(:)

  !! Test IS_PERM
  call pass_fail (1, is_perm([4,3,2,1]) .eqv. .true.)
  call pass_fail (2, is_perm([4,3,3,1]) .eqv. .false.)
  call pass_fail (3, is_perm([4,0,2,1]) .eqv. .false.)
  call pass_fail (4, is_perm([4,3,2,5]) .eqv. .false.)

  !! Test INVERSE_PERM, IDENTITY_PERM
  call pass_fail (6, all(identity_perm(7) == p1(inverse_perm(p1))))
  call pass_fail (7, all(p1 == inverse_perm(inverse_perm(p1))))
  call pass_fail (8, all(identity_perm(7) == p2(inverse_perm(p2))))
  call pass_fail (9, all(p2 == inverse_perm(inverse_perm(p2))))
  call pass_fail (10, all(identity_perm(7) == p3(inverse_perm(p3))))
  call pass_fail (11, all(p3 == inverse_perm(inverse_perm(p3))))
  call pass_fail (12, all(identity_perm(7) == inverse_perm(identity_perm(7))))

  !! Test INVERT_PERM
  allocate(p(7))
  call invert_perm (p1, p)
  call pass_fail (13, all(p == inverse_perm(p1)))
  call pass_fail (14, all(identity_perm(7) == p(p1)))
  call invert_perm (p2, p)
  call pass_fail (15, all(p == inverse_perm(p2)))
  call pass_fail (16, all(identity_perm(7) == p(p2)))
  call invert_perm (p3, p)
  call pass_fail (17, all(p == inverse_perm(p3)))
  call pass_fail (18, all(identity_perm(7) == p(p3)))

  p = p1
  call invert_perm (p)
  call pass_fail (19, all(p == inverse_perm(p1)))
  p = p2
  call invert_perm (p)
  call pass_fail (20, all(p == inverse_perm(p2)))
  p = p3
  call invert_perm (p)
  call pass_fail (21, all(p == inverse_perm(p3)))

  !! Test REORDER
  allocate(array_in(2,7), array_out(2,7), array_tmp(2,7))
  array_in(1,:) = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
  array_in(2,:) = array_in(1,:) + 1
  array_out = array_in

  call reorder (array_out, p1) ! everything is moved
  call pass_fail (22, all(array_out /= array_in))
  call pass_fail (23, maxval(array_out(2,:) - 1 - array_out(1,:)) < 2*epsilon(1.0))

  call reorder (array_out, p3)
  call reorder (array_out, inverse_perm(p3))
  call reorder (array_out, inverse_perm(p1))
  call pass_fail (24, all(array_out == array_in))

  array_out = array_in
  array_tmp = array_in
  call reorder (array_out, p3, forward=.true.)
  call reorder (array_tmp, inverse_perm(p3), forward=.false.)
  call pass_fail (25, all(array_out == array_tmp))

  !! Test ragged reordering.
  allocate(xrag1(8), rag1(14), xrag2(8), rag2(14))
  xrag1 = [1, 2, 4, 7, 8, 10, 13, 15]
  rag1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
  xrag2 = xrag1; rag2 = rag1
  call reorder (xrag1, rag1, p3)
  call pass_fail (26, all(xrag1 == [1,4,6,9,11,12,13,15]))
  call pass_fail (27, all(rag1 == [4, 5, 6, 2, 3, 10, 11, 12, 8, 9, 7, 1, 13, 14]))
  call reorder (xrag1, rag1, p3, forward=.true.)
  call pass_fail (28, all(xrag1 == xrag2))
  call pass_fail (29, all(rag1 == rag2))
 
  call exit (stat)

contains

  subroutine pass_fail (n, pass)
    use,intrinsic :: iso_fortran_env, only: error_unit
    integer, intent(in) :: n
    logical, intent(in) :: pass
    if (pass) then
      write(error_unit,'(a,i0,a)') ' Test ', n, ': PASS'
    else
      stat = 1
      write(error_unit,'(a,i0,a)') ' Test ', n, ': FAIL'
    end if
  end subroutine pass_fail

end program test_permutations
