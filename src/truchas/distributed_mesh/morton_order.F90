!!
!! MORTON_ORDER
!!
!! This module provides a subroutine for computing the Morton order
!! of an array of multi-dimensional points. The order is returned
!! as a permutation vector p, where p(1) is the index of the first
!! point, p(2) the index of the second point, etc.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! March 2026
!!

#include "f90_assert.fpp"

module morton_order

  use,intrinsic :: iso_fortran_env, only: r8 => real64, i8 => int64
  use sort_utilities, only: heap_sort
  implicit none
  private

  public get_morton_permutation

contains

  subroutine get_morton_permutation(x, p)

    real(r8), intent(in) :: x(:,:)
    integer, intent(out) :: p(:)

    integer :: dim, i, j
    real(r8), dimension(size(x,dim=1)) :: xmin, xmax, dx, xscale
    integer(i8) :: w(size(x,dim=1))
    integer(i8), allocatable :: keys(:)
    integer, parameter :: bits(2:3) = [31, 21]

    INSIST(size(x,dim=2) == size(p))

    dim = size(x,dim=1)
    INSIST(dim >= 2 .and. dim <= 3)

    !! Data bounding box
    xmin = [(minval(x(i,:)), i = 1, dim)]
    xmax = [(maxval(x(i,:)), i = 1, dim)]
    dx = xmax - xmin

    xscale = 0.0_r8
    where (dx > scale(1.0_r8, bits(dim))*tiny(1.0_r8))
      xscale = (scale(1.0_r8, bits(dim)) - 1.0_r8) / dx
    end where

    !! Generate the Morton keys
    allocate(keys(size(p)))
    select case (dim)
    case(2)
      do j = 1, size(x,dim=2)
        ! Normalize to 31-bit integers
        w = int((x(:,j)-xmin)*xscale, i8)
        ! Interleave bits: YX YX YX ...
        keys(j) = ior(spread_bits_2d(w(1)), ishft(spread_bits_2d(w(2)), 1))
      end do
    case(3)
      do j = 1, size(x,dim=2)
        ! Normalize to 21-bit integers
        w = int((x(:,j)-xmin)*xscale, i8)
        ! Interleave bits: ZYX ZYX ZYX ...
        keys(j) = ior(ior(spread_bits_3d(w(1)), ishft(spread_bits_3d(w(2)), 1)), &
                         ishft(spread_bits_3d(w(3)), 2))
      end do
    end select

    !! Index sort the keys
    call heap_sort(keys, p)

  end subroutine get_morton_permutation

  pure function spread_bits_2d(x_in) result(x)
    integer(i8), intent(in) :: x_in
    integer(i8) :: x
    x = iand(x_in, z'7FFFFFFF') ! Limit to 31 bits (Safe for signed int64)
    ! Magic bit-shifting sequence
    x = iand(ior(x, ishft(x, 16)), z'0000FFFF0000FFFF')
    x = iand(ior(x, ishft(x, 8)),  z'00FF00FF00FF00FF')
    x = iand(ior(x, ishft(x, 4)),  z'0F0F0F0F0F0F0F0F')
    x = iand(ior(x, ishft(x, 2)),  z'3333333333333333')
    x = iand(ior(x, ishft(x, 1)),  z'5555555555555555')
  end function

  pure function spread_bits_3d(x_in) result(x)
    integer(i8), intent(in) :: x_in
    integer(i8) :: x
    x = iand(x_in, z'1FFFFF') ! Limit to 21 bits (Safe for signed int64)
    ! Magic bit-shifting sequence
    x = iand(ior(x, ishft(x, 32)), z'1F00000000FFFF')
    x = iand(ior(x, ishft(x, 16)), z'1F0000FF0000FF')
    x = iand(ior(x, ishft(x, 8)),  z'100F00F00F00F00F')
    x = iand(ior(x, ishft(x, 4)),  z'10C30C30C30C30C3')
    x = iand(ior(x, ishft(x, 2)),  z'1249249249249249')
  end function

end module morton_order
