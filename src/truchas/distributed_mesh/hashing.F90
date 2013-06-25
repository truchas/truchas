!!
!!  HASHING
!!
!!  Neil N. Carlson <nnc@lanl.gov>
!!

#include "f90_assert.fpp"

module hashing

  implicit none
  private
  
  public :: initialize_hash_param, hash

  !! HASH FUNCTION OBJECTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Don't mess with these!  The bit sizes are critical to the hash function.
  integer, parameter :: i4 = selected_int_kind(9)   ! 4-byte integer, < 2^31
  integer, parameter :: i8 = selected_int_kind(18)  ! 8-byte integer, < 2^63

  !! Container for the hash function parameters.
  type, public :: hash_param
    private
    integer :: hbits  ! Number of hash address bits
    integer :: wbits  ! Number of bits in hash multiplier
    integer :: h1bit  ! Start bit for extracting hash address
    integer :: h2bit  ! Start bit for extracting hash address increment
    integer :: hsize  ! Hash address size (2**hbits)
    integer(kind=i8) :: a   ! Hash function multiplier
    integer(kind=i8) :: wmask   ! Multiplier mask
  end type hash_param

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Hash Function Procedures
 !!
 !! Variable-length key Fibonacci hash function.  See section 6.4 of Knuth's
 !! The Art of Computing Programming, Vol 3.  Note that we require the hash
 !! to be invariant with respect to permutations of the key components.  Thus
 !! we apply the *same* hash function to each key component and reduce the
 !! component hashes using exclusive-or which is commutative and associative.
 !!

  subroutine hash (hpar, key, h1, h2)

    type(hash_param), intent(in) :: hpar
    integer(kind=i4), intent(in) :: key(:)
    integer(kind=i4), intent(out) :: h1
    integer(kind=i4), intent(out), optional :: h2

    integer :: j
    integer(kind=i8) :: p
    
    ASSERT( size(key) > 0 )
    ASSERT( all(key >= 0) )

    p = key(1)*hpar%a
    do j = 2, size(key)
      p = ieor(p, key(j)*hpar%a)
    end do
    p = iand(p, hpar%wmask)

    !! Hash address in range [0,2**hpar%hbits-1]
    h1 = ibits(p, hpar%h1bit, hpar%hbits)
    ASSERT( h1 >= 0 .and. h1 < hpar%hsize )
    
    if (present(h2)) then
      !! Linear probing increment in range [1,2**hpar%hbits-1]
      h2 = ior(1_i8, ibits(p, hpar%h2bit, hpar%hbits))
      ASSERT( h2 >= 1 .and. h2 < hpar%hsize .and. modulo(h2,2) == 1 )
    end if

  end subroutine hash

  subroutine initialize_hash_param (hpar, hsize, kmax)

    type(hash_param), intent(out)   :: hpar
    integer(kind=i4), intent(inout) :: hsize ! Hash address size (minimum)
    integer(kind=i4), intent(in)    :: kmax  ! Max value of any key component

    integer(kind=i8) :: ASRC    ! leading 63 bits of 1/GoldenRatio
    data ASRC/o'474335715627751237012'/

    ASSERT( hsize > 0 )
    ASSERT( kmax > 0 )

    hpar%hbits = bit_width(hsize)                   ! Number of hash address bits
    hpar%wbits = digits(ASRC) - bit_width(kmax)     ! Number of multiplier bits
    hpar%wmask = ibset(0_i8, hpar%wbits) - 1_i8     ! Hash product mask
    hpar%a     = ishft(ASRC, -bit_width(kmax))      ! Hash multiplier
    hpar%h1bit = hpar%wbits - hpar%hbits            ! start bit for extracting hash address
    hpar%h2bit = max(0_i4, hpar%h1bit - hpar%hbits) ! start bit for extracting address increment
    !hpar%h2bit = 0 ! this seems to work better (NOT!)

    !! Return the hash address size (power of 2)
    hsize = ibset(0_i4, hpar%hbits)
    hpar%hsize = hsize

  contains

    pure integer function bit_width (i)
      integer(kind=i4), intent(in) :: i
      bit_width = digits(i)
      do while (.not.btest(i,bit_width-1))
        bit_width = bit_width - 1
        if (bit_width == 0) exit
      end do
    end function bit_width

  end subroutine initialize_hash_param
  
end module hashing
