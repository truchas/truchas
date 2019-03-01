!!
!! FACET_HASH_TYPE
!!
!! A hash function specialized for hashing facets of mesh cells.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, Jan 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type FACET_HASH which has the hash function
!!  as a type bound procedure.  This is a variable-length key Fibonacci hash
!!  function. See Knuth's The Art of Computer Programming, Vol 3, section 6.4.
!!  The key is a rank-1 integer array that describes the facet; the three node
!!  numbers of the vertices of a triangular facet, for example.  Note that we
!!  require the hash to be invariant with respect to permutations of the key
!!  components.  Thus we apply the same hash function to each key component
!!  and reduce the component hashes using exclusive-or which is commutative
!!  and associative.
!!
!!  The derived type FACET_HASH has the following type bound procedures.  All
!!  integer arguments are 32-bit (INT32 from the intrinsic ISO_FORTRAN_ENV
!!  module), and the hash function is specifically designed for that integer
!!  size.  The hash would need to be reimplemented to handle 64-bit integers.
!!
!!  INIT (HSIZE, KMAX) initializes the hash function.  HSIZE is the (mininum)
!!    desired size of the hash (one greater than the largest possible hash
!!    value).  This is rounded up to the nearest power of 2, and that value
!!    is returned in HSIZE.  KMAX is the maximum value of any key element.
!!    HSIZE must be at least as large as the number of facets (or required
!!    number of hash table entries).  The larger HSIZE, the fewer collisions
!!    and subsequent table probes will be needed, but at the expense of
!!    increased memory use.
!!
!!  HASH (KEY, H1 [,H2]) computes the hash H1 of KEY.  H1 lies in the interval
!!    [0, HSIZE-1]. The optional argument H2 is the linear probe increment to
!     use in the case of a cache collision: locations H1, H1-H2, H1-2*H2, ...
!!    should be probed in sequence, wrapping when necessary to remain within
!!    [0, HSIZE-1] (see Note 1).  The computed H1 and H2 values are invariant
!!    with respect to permutations of the values in KEY.  INIT must be called
!!    before calling this method.
!!
!! EXAMPLES
!!
!!  To hash the quadrilateral faces of a hexahedral mesh, where the number of
!!  faces is a bit more than 3 times the number of cells and at worst 6 times
!!  the number of cells (a totally disconnected mesh), one could call INIT
!!  with HSIZE equal 6 times the number of cells, and KMAX equal to the number
!!  of nodes.  INIT will redefine HSIZE to a power of 2, and this returned
!!  value should be used for the size of the table.
!!
!!  To hash the edges of a tetrahedral mesh, one could use HSIZE equal to 8
!!  or 9 times the number of nodes, as there is approximately 7 edges per node
!!  (a hex subdivided into 6 tets => 3 on edge, 3 on face, 1 interior).  KMAX
!!  would again be the number of nodes.
!!
!! NOTES
!!
!! (1) Decrementing by H2 reflects actual usage.  I do not recall at this
!!    point whether that is significant, or whether incrementing by H2 is
!!    equally good.
!!
!! (2) It may be possible to replace the internal function bit_width with
!!    one of the new F2008 bit functions.
!!

#include "f90_assert.fpp"

module facet_hash_type

  use,intrinsic :: iso_fortran_env, only: int32, int64
  implicit none
  private

  type, public :: facet_hash
    private
    integer :: hbits  ! Number of hash address bits
    integer :: wbits  ! Number of bits in hash multiplier
    integer :: h1bit  ! Start bit for extracting hash address
    integer :: h2bit  ! Start bit for extracting hash address increment
    integer :: hsize  ! Hash address size (2**hbits)
    integer(int64) :: a   ! Hash function multiplier
    integer(int64) :: wmask   ! Multiplier mask
  contains
    procedure :: init
    procedure :: hash
  end type facet_hash

contains

  subroutine hash (this, key, h1, h2)

    class(facet_hash), intent(in) :: this
    integer(int32), intent(in) :: key(:)
    integer(int32), intent(out) :: h1
    integer(int32), intent(out), optional :: h2

    integer :: j
    integer(int64) :: p

    ASSERT(size(key) > 0)
    ASSERT(all(key >= 0))

    p = key(1)*this%a
    do j = 2, size(key)
      p = ieor(p, key(j)*this%a)
    end do
    p = iand(p, this%wmask)

    !! Hash address in range [0,2**this%hbits-1]
    h1 = ibits(p, this%h1bit, this%hbits)
    ASSERT(h1 >= 0 .and. h1 < this%hsize)

    if (present(h2)) then
      !! Linear probing increment in range [1,2**this%hbits-1]
      h2 = ior(1_int64, ibits(p, this%h2bit, this%hbits))
      ASSERT(h2 >= 1 .and. h2 < this%hsize .and. modulo(h2,2) == 1)
    end if

  end subroutine hash

  subroutine init (this, hsize, kmax)

    class(facet_hash), intent(out) :: this
    integer(int32), intent(inout)  :: hsize ! hash address size (minimum)
    integer(int32), intent(in)     :: kmax  ! max value of any key component

    integer(int64) :: ASRC  ! leading 63 bits of 1/GoldenRatio
    data ASRC/o'474335715627751237012'/

    ASSERT(hsize > 0)
    ASSERT(kmax > 0)

    this%hbits = bit_width(hsize)                     ! number of hash address bits
    this%wbits = digits(ASRC) - bit_width(kmax)       ! number of multiplier bits
    this%wmask = ibset(0_int64, this%wbits) - 1_int64 ! hash product mask
    this%a     = ishft(ASRC, -bit_width(kmax))        ! hash multiplier
    this%h1bit = this%wbits - this%hbits              ! start bit for extracting hash address
    this%h2bit = max(0_int32, this%h1bit - this%hbits)! start bit for extracting address increment
    !this%h2bit = 0 ! this seems to work better (NOT!)

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

end module facet_hash_type
