!!
!! The PERMUTATIONS Module
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! Last revised 10 Mar 2004; initial version Nov 2002.
!!
!! This module provides procedures for manipulating permutations and doing
!! in-place reordering of arrays according to a provided permutation.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2002, 2004 Neil N. Carlson, Carlson Science Computing
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! All procedures are pure.
!!
!! IS_PERM(PERM) returns the value true when the rank-1 default integer array
!! PERM describes a valid permutation. That is, PERM is a 1-1 mapping of
!! the set of integers {1, 2, ..., SIZE(PERM)} onto itself. Otherwise it
!! returns the value false. Implicit in this definition is that
!! LBOUND(PERM,1) equals 1.
!!
!! CALL INVERT_PERM (PERM [, INVP]) returns the inverse permutation of the
!! permutation PERM. PERM is a rank-1 default integer array. If the rank-1
!! default integer array INVP is present, it has INTENT(OUT) and the inverse
!! permutation is returned through it. In this case PERM has INTENT(IN).
!! If INVP is not present, then PERM has INTENT(INOUT) and its value is
!! overwritten with the inverse permutation. IS_PERM(PERM) must be true,
!! and if present, INVP must have the same array bounds as PERM, otherwise
!! the results are indeterminant.
!!
!! INVERSE_PERM(PERM) returns a rank-1 default integer array holding the
!! inverse of the permutation described by the rank-1 default integer
!! array PERM. The size of the result equals the size of PERM. The
!! results are indeterminant unless IS_PERM(PERM) is true.
!!
!! IDENTITY_PERM(N) returns a rank-1 default integer array of size N that
!! describes the identity permutation. That is, the result is the array
!! (/ 1, 2, ..., N /). The default integer N must be positive.
!!
!! IS_IDENTITY_PERM(PERM) returns the value true when the rank-1 default
!! integer array PERM describes the identity permutation.  Otherwise it
!! returns the value false.
!!
!! CALL REORDER (ARRAY, PERM [,FORWARD]) reorders the values in the array
!! ARRAY according to the permutation described by the rank-1 default
!! integer array PERM. PERM has INTENT(IN) and ARRAY has INTENT(INOUT).
!! If ARRAY is a multi-dimensional array, the permutation is performed along
!! the last dimension, and the array bounds in that dimension must equal
!! those of PERM. If the optional argument FORWARD is present and has value
!! true, the array is pushed-forward by the permutation. In the case of a
!! two-dimensional array, this means that the returned values
!! ARRAY(:,PERM(j)) equal the input values ARRAY(:,j). Otherwise the array
!! is pulled-back by the permutation; that is, the returned values ARRAY(:,j)
!! equal the input values ARRAY(:,PERM(j)). The results are indeterminant
!! unless IS_PERM(PERM) is true and the array bounds in the last dimension
!! of ARRAY equal those of PERM.
!!
!!
!! IMPLEMENTATION NOTES
!!
!! o REORDER is currently defined for rank-1 and rank-2 arrays of default
!! logical, integer, and real kinds, and for double precision kinds. The
!! specific versions were generated using commented template code at the
!! end of this file. This can be used to generate additional specific
!! versions as needed. (I hate having the actual code in a separate include
!! file and then including it multiple times. Better to go to m4 macro
!! preprocessing which would permit the template code to exist in-line.)
!!
!! o It should be relatively straightforward to extend the code to handle a
!! more general definition of a permutation as a 1-1 mapping of
!! {LBOUND(PERM),... UBOUND(PERM)} onto itself. This might be useful.
!! For this reason I've refered to the array bounds when describing the'
!! constraints wherever possible.
!!
!! o Some of the routines require a tag array to mark values that have been
!! dealt with. The permutation vector could itself be used for this
!! purpose by toggling the sign of the value (assuming permutations of
!! {1,...,N}), and then restoring the sign at the end of the routine.
!! However, this would require making the PERM argument INTENT(INOUT)
!! instead of INTENT(IN) which may cause havoc with the caller from whose
!! perspective the argument ought to be INTENT(IN). For this reason I've'
!! chosen to allocate a temporary tag array instead.
!!
!! o The REORDER routines do their work in place; there is no copying of
!! the array to new storage. The procedure is based on the result that
!! any permutation can be decomposed as a sequence of independent cycles.
!! The array values are simply moved around each cycle. The pull-back
!! version of reorder (the default) is optimal, in that requires the
!! minimum amount of data movement possible. The push-forward version is
!! less efficient, requiring swapping data with a temporary at each step.
!! Instead of pushing-forward an array, in may be more efficient to do the
!! equivalent operation of pulling-back the array by the inverse of the
!! permutation.
!!

module permutations

  implicit none
  private

  public :: is_perm, invert_perm, inverse_perm, identity_perm, is_identity_perm, reorder

  interface invert_perm
    module procedure invert_perm_1, invert_perm_2
  end interface

  interface reorder
    module procedure reorder_l0, reorder_l1
    module procedure reorder_i0, reorder_i1
    module procedure reorder_r0, reorder_r1
    module procedure reorder_d0, reorder_d1
  end interface

contains

  logical function is_perm (perm)

    integer, intent(in) :: perm(:)

    integer :: j
    logical :: image_value(size(perm))

    if (size(perm) < 1) then

      is_perm = .true.

    else ! Verify that PERM is 1-1 onto {1, 2, ..., SIZE(PERM)}.

      is_perm = .true.
      image_value = .false.
      do j = 1, size(perm)
        if (perm(j) < 1 .or. perm(j) > size(perm)) then
          is_perm = .false.
          exit
        else if (image_value(perm(j))) then ! PERM is not 1-1.
          is_perm = .false.
          exit
        end if
        image_value(perm(j)) = .true.
      end do

    end if

  end function is_perm

  logical function is_identity_perm (perm)

    integer, intent(in) :: perm(:)

    integer :: j

    is_identity_perm = .false.
    if (size(perm) < 1) return

    do j = 1, size(perm)
      if (perm(j) /= j) return
    end do
    is_identity_perm = .true.

  end function is_identity_perm

  pure subroutine invert_perm_1 (perm)

    integer, intent(inout) :: perm(:)

    integer :: k, j, pk, ppk

    do j = 1, size(perm)
      if (perm(j) < 0 .or. perm(j) == j) cycle
      k = j
      pk = perm(k)
      do while (pk /= j)
        ppk = perm(pk)
        perm(pk) = -k
        k = pk
        pk = ppk
      end do
      perm(j) = k
    end do

    perm = abs(perm)

  end subroutine invert_perm_1

  pure subroutine invert_perm_2 (perm, invp)

    integer, intent(in) :: perm(:)
    integer, intent(out) :: invp(:)

    integer :: j

    do j = 1, size(perm)
      invp(perm(j)) = j
    end do

  end subroutine invert_perm_2

  pure function inverse_perm (perm) result (invp)

    integer, intent(in) :: perm(:)
    integer :: invp(size(perm))

    integer :: j

    do j = 1, size(perm)
      invp(perm(j)) = j
    end do

  end function inverse_perm

  pure function identity_perm (n) result (perm)
    integer, intent(in) :: n
    integer :: perm(n)
    integer :: j
    do j = 1, n
      perm(j) = j
    end do
  end function identity_perm


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! REORDER -- specific procedures
 !!

  !!
  !! Rank-1 logical array
  !!
  subroutine reorder_l0 (vector, perm, forward)

    logical, dimension(:), intent(inout) :: vector
    integer, intent(in) :: perm(:)
    logical, intent(in), optional :: forward

    logical :: vk, vpk, vj
    integer :: j, k, pk
    logical :: pullback, tagged(size(perm))

    pullback = .true.
    if (present(forward)) pullback = .not.forward

    if (pullback) then ! Pull VECTOR back by PERM.

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vj = vector(j)
        k = j
        do while (perm(k) /= j)
          vector(k) = vector(perm(k))
          tagged(k) = .true. ! Tag vector(k) as defined.
          k = perm(k) ! Advance.
        end do
        vector(k) = vj
        tagged(k) = .true.
      end do

    else ! Push VECTOR forward by PERM

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vk = vector(j)
        pk = perm(j)
        do while (pk /= j)
          vpk = vector(pk) ! Swap vk and vector(pk).
          vector(pk) = vk
          vk = vpk
          tagged(pk) = .true. ! Tag vector(pk) as moved.
          pk = perm(pk) ! Advance.
        end do
        vector(j) = vk
      end do

    end if

  end subroutine reorder_l0
  !!
  !! Rank-2 logical array
  !!
  subroutine reorder_l1 (vector, perm, forward)

    logical, dimension(:,:), intent(inout) :: vector
    integer, intent(in) :: perm(:)
    logical, intent(in), optional :: forward

    logical, dimension(size(vector,1)) :: vk, vpk, vj
    integer :: j, k, pk
    logical :: pullback, tagged(size(perm))

    pullback = .true.
    if (present(forward)) pullback = .not.forward

    if (pullback) then ! Pull VECTOR back by PERM.

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vj = vector(:,j)
        k = j
        do while (perm(k) /= j)
          vector(:,k) = vector(:,perm(k))
          tagged(k) = .true. ! Tag vector(k) as defined.
          k = perm(k) ! Advance.
        end do
        vector(:,k) = vj
        tagged(k) = .true.
      end do

    else ! Push VECTOR forward by PERM

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vk = vector(:,j)
        pk = perm(j)
        do while (pk /= j)
          vpk = vector(:,pk) ! Swap vk and vector(pk).
          vector(:,pk) = vk
          vk = vpk
          tagged(pk) = .true. ! Tag vector(pk) as moved.
          pk = perm(pk) ! Advance.
        end do
        vector(:,j) = vk
      end do

    end if

  end subroutine reorder_l1
  !!
  !! Rank-1 integer array
  !!
  subroutine reorder_i0 (vector, perm, forward)

    integer, dimension(:), intent(inout) :: vector
    integer, intent(in) :: perm(:)
    logical, intent(in), optional :: forward

    integer :: vk, vpk, vj
    integer :: j, k, pk
    logical :: pullback, tagged(size(perm))

    pullback = .true.
    if (present(forward)) pullback = .not.forward

    if (pullback) then ! Pull VECTOR back by PERM.

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vj = vector(j)
        k = j
        do while (perm(k) /= j)
          vector(k) = vector(perm(k))
          tagged(k) = .true. ! Tag vector(k) as defined.
          k = perm(k) ! Advance.
        end do
        vector(k) = vj
        tagged(k) = .true.
      end do

    else ! Push VECTOR forward by PERM

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vk = vector(j)
        pk = perm(j)
        do while (pk /= j)
          vpk = vector(pk) ! Swap vk and vector(pk).
          vector(pk) = vk
          vk = vpk
          tagged(pk) = .true. ! Tag vector(pk) as moved.
          pk = perm(pk) ! Advance.
        end do
        vector(j) = vk
      end do

    end if

  end subroutine reorder_i0
  !!
  !! Rank-2 integer array
  !!
  subroutine reorder_i1 (vector, perm, forward)

    integer, dimension(:,:), intent(inout) :: vector
    integer, intent(in) :: perm(:)
    logical, intent(in), optional :: forward

    integer, dimension(size(vector,1)) :: vk, vpk, vj
    integer :: j, k, pk
    logical :: pullback, tagged(size(perm))

    pullback = .true.
    if (present(forward)) pullback = .not.forward

    if (pullback) then ! Pull VECTOR back by PERM.

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vj = vector(:,j)
        k = j
        do while (perm(k) /= j)
          vector(:,k) = vector(:,perm(k))
          tagged(k) = .true. ! Tag vector(k) as defined.
          k = perm(k) ! Advance.
        end do
        vector(:,k) = vj
        tagged(k) = .true.
      end do

    else ! Push VECTOR forward by PERM

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vk = vector(:,j)
        pk = perm(j)
        do while (pk /= j)
          vpk = vector(:,pk) ! Swap vk and vector(pk).
          vector(:,pk) = vk
          vk = vpk
          tagged(pk) = .true. ! Tag vector(pk) as moved.
          pk = perm(pk) ! Advance.
        end do
        vector(:,j) = vk
      end do

    end if

  end subroutine reorder_i1
  !!
  !! Rank-1 real array
  !!
  subroutine reorder_r0 (vector, perm, forward)

    real, dimension(:), intent(inout) :: vector
    integer, intent(in) :: perm(:)
    logical, intent(in), optional :: forward

    real :: vk, vpk, vj
    integer :: j, k, pk
    logical :: pullback, tagged(size(perm))

    pullback = .true.
    if (present(forward)) pullback = .not.forward

    if (pullback) then ! Pull VECTOR back by PERM.

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vj = vector(j)
        k = j
        do while (perm(k) /= j)
          vector(k) = vector(perm(k))
          tagged(k) = .true. ! Tag vector(k) as defined.
          k = perm(k) ! Advance.
        end do
        vector(k) = vj
        tagged(k) = .true.
      end do

    else ! Push VECTOR forward by PERM

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vk = vector(j)
        pk = perm(j)
        do while (pk /= j)
          vpk = vector(pk) ! Swap vk and vector(pk).
          vector(pk) = vk
          vk = vpk
          tagged(pk) = .true. ! Tag vector(pk) as moved.
          pk = perm(pk) ! Advance.
        end do
        vector(j) = vk
      end do

    end if

  end subroutine reorder_r0
  !!
  !! Rank-2 real array
  !!
  subroutine reorder_r1 (vector, perm, forward)

    real, dimension(:,:), intent(inout) :: vector
    integer, intent(in) :: perm(:)
    logical, intent(in), optional :: forward

    real, dimension(size(vector,1)) :: vk, vpk, vj
    integer :: j, k, pk
    logical :: pullback, tagged(size(perm))

    pullback = .true.
    if (present(forward)) pullback = .not.forward

    if (pullback) then ! Pull VECTOR back by PERM.

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vj = vector(:,j)
        k = j
        do while (perm(k) /= j)
          vector(:,k) = vector(:,perm(k))
          tagged(k) = .true. ! Tag vector(k) as defined.
          k = perm(k) ! Advance.
        end do
        vector(:,k) = vj
        tagged(k) = .true.
      end do

    else ! Push VECTOR forward by PERM

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vk = vector(:,j)
        pk = perm(j)
        do while (pk /= j)
          vpk = vector(:,pk) ! Swap vk and vector(pk).
          vector(:,pk) = vk
          vk = vpk
          tagged(pk) = .true. ! Tag vector(pk) as moved.
          pk = perm(pk) ! Advance.
        end do
        vector(:,j) = vk
      end do

    end if

  end subroutine reorder_r1
  !!
  !! Rank-1 double precision array
  !!
  subroutine reorder_d0 (vector, perm, forward)

    double precision, dimension(:), intent(inout) :: vector
    integer, intent(in) :: perm(:)
    logical, intent(in), optional :: forward

    double precision :: vk, vpk, vj
    integer :: j, k, pk
    logical :: pullback, tagged(size(perm))

    pullback = .true.
    if (present(forward)) pullback = .not.forward

    if (pullback) then ! Pull VECTOR back by PERM.

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vj = vector(j)
        k = j
        do while (perm(k) /= j)
          vector(k) = vector(perm(k))
          tagged(k) = .true. ! Tag vector(k) as defined.
          k = perm(k) ! Advance.
        end do
        vector(k) = vj
        tagged(k) = .true.
      end do

    else ! Push VECTOR forward by PERM

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vk = vector(j)
        pk = perm(j)
        do while (pk /= j)
          vpk = vector(pk) ! Swap vk and vector(pk).
          vector(pk) = vk
          vk = vpk
          tagged(pk) = .true. ! Tag vector(pk) as moved.
          pk = perm(pk) ! Advance.
        end do
        vector(j) = vk
      end do

    end if

  end subroutine reorder_d0
  !!
  !! Rank-2 double precision array
  !!
  subroutine reorder_d1 (vector, perm, forward)

    double precision, dimension(:,:), intent(inout) :: vector
    integer, intent(in) :: perm(:)
    logical, intent(in), optional :: forward

    double precision, dimension(size(vector,1)) :: vk, vpk, vj
    integer :: j, k, pk
    logical :: pullback, tagged(size(perm))

    pullback = .true.
    if (present(forward)) pullback = .not.forward

    if (pullback) then ! Pull VECTOR back by PERM.

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vj = vector(:,j)
        k = j
        do while (perm(k) /= j)
          vector(:,k) = vector(:,perm(k))
          tagged(k) = .true. ! Tag vector(k) as defined.
          k = perm(k) ! Advance.
        end do
        vector(:,k) = vj
        tagged(k) = .true.
      end do

    else ! Push VECTOR forward by PERM

      tagged = .false.
      do j = 1, size(perm)
        if (tagged(j) .or. perm(j) == j) cycle
        vk = vector(:,j)
        pk = perm(j)
        do while (pk /= j)
          vpk = vector(:,pk) ! Swap vk and vector(pk).
          vector(:,pk) = vk
          vk = vpk
          tagged(pk) = .true. ! Tag vector(pk) as moved.
          pk = perm(pk) ! Advance.
        end do
        vector(:,j) = vk
      end do

    end if

  end subroutine reorder_d1
end module permutations

#ifdef UNIT_TEST

program test_permutations
  use permutations

  integer, allocatable :: p(:)
  integer :: p1(7) = (/ 2, 5, 7, 1, 3, 4, 6 /) ! 1 cycle
  integer :: p2(7) = (/ 3, 6, 5, 2, 7, 4, 1 /) ! 2 cycle
  integer :: p3(7) = (/ 3, 2, 6, 5, 4, 1, 7 /) ! 2 cycle w/ fixed pt and swap

  real, allocatable :: array_in(:,:), array_out(:,:), array_tmp(:,:)

  !! Test IS_PERM
  call pass_fail (1, is_perm((/4,3,2,1/)) .eqv. .true.)
  call pass_fail (2, is_perm((/4,3,3,1/)) .eqv. .false.)
  call pass_fail (3, is_perm((/4,0,2,1/)) .eqv. .false.)
  call pass_fail (4, is_perm((/4,3,2,5/)) .eqv. .false.)

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
  array_in(1,:) = (/ 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7 /)
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

contains

  subroutine pass_fail (n, pass)
    integer, intent(in) :: n
    logical, intent(in) :: pass
    if (pass) then
      print *, ' Test', n, ': PASS'
    else
      print *, ' Test', n, ': FAIL'
    end if
  end subroutine pass_fail

end program test_permutations

#endif

! Here is code before preprocessing:
!
!  !!
!  !! Rank-1 logical array
!  !!
!#define _ROUTINE_NAME_  reorder_l0
!#define _DATA_TYPE_     logical
!#define _DATA_RANK_     0
!#include "reorder.fpp"
!
!  !!
!  !! Rank-2 logical array
!  !!
!#define _ROUTINE_NAME_  reorder_l1
!#define _DATA_TYPE_     logical
!#define _DATA_RANK_     1
!#include "reorder.fpp"
!
!  !!
!  !! Rank-1 integer array
!  !!
!#define _ROUTINE_NAME_  reorder_i0
!#define _DATA_TYPE_     integer
!#define _DATA_RANK_     0
!#include "reorder.fpp"
!
!  !!
!  !! Rank-2 integer array
!  !!
!#define _ROUTINE_NAME_  reorder_i1
!#define _DATA_TYPE_     integer
!#define _DATA_RANK_     1
!#include "reorder.fpp"
!
!  !!
!  !! Rank-1 real array
!  !!
!#define _ROUTINE_NAME_  reorder_r0
!#define _DATA_TYPE_     real
!#define _DATA_RANK_     0
!#include "reorder.fpp"
!
!  !!
!  !! Rank-2 real array
!  !!
!#define _ROUTINE_NAME_  reorder_r1
!#define _DATA_TYPE_     real
!#define _DATA_RANK_     1
!#include "reorder.fpp"
!
!  !!
!  !! Rank-1 double precision array
!  !!
!#define _ROUTINE_NAME_  reorder_d0
!#define _DATA_TYPE_     double precision
!#define _DATA_RANK_     0
!#include "reorder.fpp"
!
!  !!
!  !! Rank-2 double precision array
!  !!
!#define _ROUTINE_NAME_  reorder_d1
!#define _DATA_TYPE_     double precision
!#define _DATA_RANK_     1
!#include "reorder.fpp"
!
!
! Here is the contents of the file "reorder.fpp":
!
!#ifndef _ROUTINE_NAME_
!#error  "_ROUTINE_NAME_ must be defined before including this file"
!#endif
!
!#ifndef _DATA_TYPE_
!#error  "_DATA_TYPE_ must be defined before including this file"
!#endif
!
!#ifndef _DATA_RANK_
!#error "_DATA_RANK_ must be defined before including this file"
!#endif
!
!#if _DATA_RANK_ == 0
!#define _VEC_DECL_    _DATA_TYPE_, dimension(:)
!#define _TMP_DECL_    _DATA_TYPE_
!#define _VEC_VAL_(i)  vector(i)
!#elif _DATA_RANK_ == 1
!#define _VEC_DECL_    _DATA_TYPE_, dimension(:,:)
!#define _TMP_DECL_    _DATA_TYPE_, dimension(size(vector,1))
!#define _VEC_VAL_(i)  vector(:,i)
!#elif _DATA_RANK_ == 2
!#define _VEC_DECL_    _DATA_TYPE_, dimension(:,:,:)
!#define _TMP_DECL_    _DATA_TYPE_, dimension(size(vector,1),size(vector,2))
!#define _VEC_VAL_(i)  vector(:,:,i)
!#else
!#error "_DATA_RANK_ must be 0, 1, or 2"
!#endif
!
!  subroutine _ROUTINE_NAME_ (vector, perm, forward)
!
!    _VEC_DECL_, intent(inout) :: vector
!    integer, intent(in) :: perm(:)
!    logical, intent(in), optional :: forward
!
!    _TMP_DECL_ :: vk, vpk, vj
!    integer :: j, k, pk
!    logical :: pullback, tagged(size(perm))
!
!    pullback = .true.
!    if (present(forward)) pullback = .not.forward
!
!    if (pullback) then ! Pull VECTOR back by PERM.
!
!      tagged = .false.
!      do j = 1, size(perm)
!        if (tagged(j) .or. perm(j) == j) cycle
!        vj = _VEC_VAL_(j)
!        k = j
!        do while (perm(k) /= j)
!          _VEC_VAL_(k) = _VEC_VAL_(perm(k))
!          tagged(k) = .true.    ! Tag vector(k) as defined.
!          k = perm(k)           ! Advance.
!        end do
!        _VEC_VAL_(k) = vj
!        tagged(k) = .true.
!      end do
!
!    else ! Push VECTOR forward by PERM
!
!      tagged = .false.
!      do j = 1, size(perm)
!        if (tagged(j) .or. perm(j) == j) cycle
!        vk = _VEC_VAL_(j)
!        pk = perm(j)
!        do while (pk /= j)
!          vpk = _VEC_VAL_(pk)    ! Swap vk and vector(pk).
!          _VEC_VAL_(pk) = vk
!          vk = vpk
!          tagged(pk) = .true.   ! Tag vector(pk) as moved.
!          pk = perm(pk)         ! Advance.
!        end do
!        _VEC_VAL_(j) = vk
!      end do
!
!    end if
!
!  end subroutine _ROUTINE_NAME_
!    
!#undef _ROUTINE_NAME_
!#undef _DATA_TYPE_
!#undef _DATA_RANK_
!#undef _VEC_DECL_
!#undef _TMP_DECL_
!#undef _VEC_VAL_
