!!
!! EDGE_HASH_TYPE
!!
!! A hash function specialized for hashing edges of enclosure faces.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 7 May 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type EDGE_HASH which has the hash function as
!!  a type bound procedure.  This is a two-round xorshift-multiply-xorshift hash
!!  function discovered by the Hash Function Prospector project. See the
!!  project's GitHub repository github.com/skeeto/hash-prospector.
!!
!!  The key is a rank-1 integer array of length 2 that describes the edge.  Note
!!  that we require the hash to be invariant with respect to permutations of the
!!  key components.  Thus we apply the same hash function to each key component
!!  and reduce the component hashes using exclusive-or which is commutative and
!!  associative. We then hash the result to reduce the bias inherent to the
!!  exclusive-or reduction process.
!!
!!  The derived type EDGE_HASH has the following type bound procedures.  All
!!  integer arguments are 32-bit (INT32 from the intrinsic ISO_FORTRAN_ENV
!!  module), and the hash function is specifically designed for that integer
!!  size.  The hash would need to be reimplemented to handle 64-bit integers.
!!
!!  INIT (HSIZE, HSIZE) initializes the hash function.  HSIZE is the (minimum)
!!    desired size of the hash (one greater than the largest possible hash
!!    value).  This is rounded up to the nearest power of 2, and that value is
!!    returned in HSIZE.  HSIZE must be at least as large as the number of
!!    edges (or required number of hash table entries).  The larger HSIZE, the
!!    fewer collisions and subsequent table probes will be needed, but at the
!!    expense of increased memory use.
!!
!!  HASH (KEY, H1) computes the hash H1 of KEY.  H1 lies in the interval
!!    [0, HSIZE-1].  The computed H1 value is invariant with respect to
!!    permutations of the values in KEY.  INIT must be called before calling
!!    this method.
!!
!! EXAMPLES
!!
!!  To hash the edges of a hexahedral mesh, where the number of edges is a bit
!!  more than 2 times the number of faces and at worst 4 times the number of
!!  faces (a totally disconnected mesh), one could call INIT with HSIZE equal 4
!!  times the number of faces.  INIT will redefine HSIZE to a power of 2, and
!!  this returned value should be used for the size of the table.
!!
!!  To hash the edges of a tetrahedral mesh, one could use HSIZE equal to 2
!!  times the number of faces, as there is approximately 2 edges per face.
!!
!! NOTES
!!
!! (1) It may be possible to replace the internal function bit_width with
!!    one of the new F2008 bit functions.
!!

#include "f90_assert.fpp"

module edge_hash_type

  use,intrinsic :: iso_fortran_env, only: int32, int64
  implicit none
  private

  type, public :: edge_hash
    private
    integer :: hbits  ! Number of hash address bits
    integer :: hsize  ! Hash address size (2**hbits)
  contains
    procedure :: init
    procedure :: hash
  end type edge_hash

contains

  !! Source: https://github.com/skeeto/hash-prospector
  !! NB: c2 = 0x846ca68b but the standard (as of F2018) first converts BOZ
  !!     constants to signed integers, which c2 overflows.  The workaround
  !!     below works for machines with a two's-complement representation of
  !!     signed integers.
  pure integer(int32) function integer_hash (x) result(ret)

    integer(int32), intent(in) :: x
    integer(int32), parameter :: c1 =  int(Z'7feb352d',int32)
    integer(int32), parameter :: c2 = -int(Z'7b935975',int32)

    ret = x
    ret = ieor(ret, ishft(ret, -16))
    ret = ret * c1
    ret = ieor(ret, ishft(ret, -15))
    ret = ret * c2
    ret = ieor(ret, ishft(ret, -16))

  end function integer_hash


  subroutine hash (this, key, h1)

    class(edge_hash), intent(in) :: this
    integer(int32), intent(in) :: key(2)
    integer(int32), intent(out) :: h1

    integer :: j
    integer(int32) :: hc

    ASSERT(all(key >= 0))

    hc = 0
    do j = 1, size(key)
      hc = ieor(hc, integer_hash(key(j)))
    end do

    !! Hash again to reduce XOR bias
    hc = integer_hash(hc)

    !! Hash address in range [0,2**this%hbits-1]
    h1 = iand(hc, this%hsize - 1)

  end subroutine hash


  subroutine init (this, hsize)

    class(edge_hash), intent(out) :: this
    integer(int32), intent(inout)  :: hsize ! hash address size (minimum)

    ASSERT(hsize > 0)

    this%hbits = bit_width(hsize)  ! number of hash address bits

    !! Return the hash address size (power of 2)
    hsize = ibset(0_int32, this%hbits)
    this%hsize = hsize

  contains

    pure integer function bit_width (i)
      integer(int32), intent(in) :: i
      bit_width = digits(i)
      do while (.not.btest(i,bit_width-1))
        bit_width = bit_width - 1
        if (bit_width == 0) exit
      end do
    end function bit_width

  end subroutine init

end module edge_hash_type
