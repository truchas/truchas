!!
!! BITFIELD_TYPE: For when you need more bits than fit into an integer.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 26 May 2010
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bitfield_type

  implicit none
  private

  !! A derived type that stores a large, but fixed, number of bits.
  !! As for integers, bits are indexed 0 through BIT_SIZE(<bitfield>)-1.
  !! Internally, the bits are stored in a fixed-size integer array.
  public :: bitfield

  !! Intrinsic bit manipulation procedures extended to the BITFIELD type.
  !! Procedures are elemental where the instrinsics are.
  public :: bit_size, btest, ibset, ibclr, iand, ieor, ior, not
  public :: popcnt, trailz

  !! Equality/inequality operators overloaded for BITFIELD operands.
  public :: operator(==), operator(/=)

  !! Constant bitfield value with all bits set to 0.
  public :: ZERO_BITFIELD

  !! Generics from PARALLEL_COMMUNICATION extended to the BITFIELD type.
  !! Implementations for rank-1 BITFIELD arrays only.
  public :: distribute, allocate_collated_array

  !! Generics from INDEX_PARTITIONING extended to the BITFIELD type.
  !! Implementations for rank-1 BITFIELD arrays only.
  public :: gather_boundary
  
  !! Parallel extensions of some intrinsic bit manipulation procedures.
  public :: global_ior

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: NUM_CHUNK = 4 ! increase this if you need more bits
  integer, parameter :: CHUNK_BITS = bit_size(1)  ! bits in a default integer
  integer, parameter :: BITFIELD_BITS = NUM_CHUNK * CHUNK_BITS

  type :: bitfield
    private
    integer :: chunk(0:NUM_CHUNK-1)
  end type bitfield

  type(bitfield), parameter :: ZERO_BITFIELD = bitfield(0)

  interface bit_size
    module procedure bit_size_bitfield, bit_size_bitfield1
  end interface

  interface ibset
    module procedure ibset_bitfield
  end interface

  interface ibclr
    module procedure ibclr_bitfield
  end interface

  interface btest
    module procedure btest_bitfield
  end interface

  interface iand
    module procedure iand_bitfield
  end interface

  interface ieor
    module procedure ieor_bitfield
  end interface

  interface ior
    module procedure ior_bitfield
  end interface

  interface not
    module procedure not_bitfield
  end interface

  interface popcnt
    module procedure popcnt_bitfield
  end interface

  interface trailz
    module procedure trailz_bitfield
  end interface

  interface operator(==)
    module procedure eq_bitfield
  end interface

  interface operator(/=)
    module procedure ne_bitfield
  end interface

  interface distribute
    module procedure distribute_bitfield
  end interface

  interface allocate_collated_array
    module procedure alloc_coll_array_bitfield
  end interface

  interface gather_boundary
    module procedure gather_boundary_bitfield1, gather_boundary_bitfield2
  end interface

contains

  integer function bit_size_bitfield (bf)
    type(bitfield), intent(in) :: bf
    bit_size_bitfield = BITFIELD_BITS
  end function bit_size_bitfield

  integer function bit_size_bitfield1 (bf)
    type(bitfield), intent(in) :: bf(:)
    bit_size_bitfield1 = BITFIELD_BITS
  end function bit_size_bitfield1

  elemental logical function btest_bitfield (bf, pos)
    type(bitfield), intent(in) :: bf
    integer, intent(in) :: pos
    integer :: n, p
    n = pos / CHUNK_BITS
    p = modulo(pos, CHUNK_BITS)
    btest_bitfield = btest(bf%chunk(n),p)
  end function btest_bitfield

  elemental function ibset_bitfield (bf, pos) result (bf_out)
    type(bitfield), intent(in) :: bf
    integer, intent(in) :: pos
    type(bitfield) :: bf_out
    integer :: n, p
    bf_out = bf
    n = pos / CHUNK_BITS
    p = modulo(pos, CHUNK_BITS)
    bf_out%chunk(n) = ibset(bf%chunk(n),p)
  end function ibset_bitfield

  elemental function ibclr_bitfield (bf, pos) result (bf_out)
    type(bitfield), intent(in) :: bf
    integer, intent(in) :: pos
    type(bitfield) :: bf_out
    integer :: n, p
    bf_out = bf
    n = pos / CHUNK_BITS
    p = modulo(pos, CHUNK_BITS)
    bf_out%chunk(n) = ibclr(bf%chunk(n),p)
  end function ibclr_bitfield

  elemental function iand_bitfield (bf1, bf2) result (bf_out)
    type(bitfield), intent(in) :: bf1, bf2
    type(bitfield) :: bf_out
    bf_out%chunk = iand(bf1%chunk, bf2%chunk)
  end function iand_bitfield

  elemental function ior_bitfield (bf1, bf2) result (bf_out)
    type(bitfield), intent(in) :: bf1, bf2
    type(bitfield) :: bf_out
    bf_out%chunk = ior(bf1%chunk, bf2%chunk)
  end function ior_bitfield

  elemental function ieor_bitfield (bf1, bf2) result (bf_out)
    type(bitfield), intent(in) :: bf1, bf2
    type(bitfield) :: bf_out
    bf_out%chunk = ieor(bf1%chunk, bf2%chunk)
  end function ieor_bitfield

  elemental function not_bitfield (bf) result (bf_out)
    type(bitfield), intent(in) :: bf
    type(bitfield) :: bf_out
    bf_out%chunk = not(bf%chunk)
  end function not_bitfield

  elemental integer function popcnt_bitfield (bf)
    type(bitfield), intent(in) :: bf
    popcnt_bitfield = sum(popcnt(bf%chunk))
  end function popcnt_bitfield

  elemental integer function trailz_bitfield (bf) result (p)
    type(bitfield), intent(in) :: bf
    integer :: n, q
    p = 0
    do n = 0, NUM_CHUNK-1
      q = trailz(bf%chunk(n))
      p = p + q
      if (q < CHUNK_BITS) exit
    end do
  end function trailz_bitfield

  elemental logical function eq_bitfield (bf1, bf2)
    type(bitfield), intent(in) :: bf1, bf2
    eq_bitfield = all(bf1%chunk == bf2%chunk)
  end function eq_bitfield

  elemental logical function ne_bitfield (bf1, bf2)
    type(bitfield), intent(in) :: bf1, bf2
    ne_bitfield = any(bf1%chunk /= bf2%chunk)
  end function ne_bitfield

  subroutine distribute_bitfield (vout, vin, bsize)

    use parallel_communication, only: distribute

    type(bitfield), intent(out) :: vout(:)
    type(bitfield), intent(in)  :: vin(:)
    integer, intent(in), optional :: bsize(:)

    integer :: n

    do n = 0, NUM_CHUNK-1
      call distribute (vout%chunk(n), vin%chunk(n), bsize)
    end do

  end subroutine distribute_bitfield

  subroutine alloc_coll_array_bitfield (array, size1, stat)

    use parallel_communication, only: is_IOP

    type(bitfield), pointer :: array(:)
    integer, intent(in) :: size1
    integer, intent(out), optional :: stat

    if (is_IOP) then
      if (present(stat)) then
        allocate(array(size1), stat=stat)
      else
        allocate(array(size1))
      end if
    else
      allocate(array(0))
    end if

  end subroutine alloc_coll_array_bitfield

  subroutine gather_boundary_bitfield1 (this, local_data)

    use index_partitioning, only: ip_desc, gather_boundary

    type(ip_desc),  intent(in)    :: this   ! partition descriptor
    type(bitfield), intent(inout) :: local_data(:)  ! local data array

    integer :: n

    do n = 0, NUM_CHUNK-1
      call gather_boundary (this, local_data%chunk(n))
    end do

  end subroutine gather_boundary_bitfield1

  subroutine gather_boundary_bitfield2 (this, onP_data, offP_data)

    use index_partitioning, only: ip_desc, gather_boundary

    type(ip_desc),  intent(in)  :: this         ! partition descriptor
    type(bitfield), intent(in)  :: onP_data(:)  ! on-process data array
    type(bitfield), intent(out) :: offP_data(:) ! off-process data array

    integer :: n

    do n = 0, NUM_CHUNK-1
      call gather_boundary (this, onP_data%chunk(n), offP_data%chunk(n))
    end do

  end subroutine gather_boundary_bitfield2

  function global_ior (bf) result (bf_out)

    use parallel_communication, only: nPE, is_IOP, collate, broadcast

    type(bitfield), intent(in) :: bf
    type(bitfield) :: bf_out

    integer :: n, p, tmp, array(nPE)

    do n = 0, NUM_CHUNK-1
      call collate (array, bf%chunk(n))
      if (is_IOP) then
        tmp = array(1)
        do p = 2, nPE
          tmp = ior(tmp, array(p))
        end do
        bf_out%chunk(n) = tmp
      end if
    end do

    call broadcast (bf_out%chunk)

  end function global_ior

end module bitfield_type
