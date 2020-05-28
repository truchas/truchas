!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! distribute -- include file for parallel_communication.F90

#ifdef _TYPE_
#undef _TYPE_
#endif

#ifdef _INT8_DATA_
#undef _INT8_DATA_
#define _TYPE_ integer(int8)
#define _PROC2_ distribute_int8_2
#endif

#ifdef _INTEGER_DATA_
#undef _INTEGER_DATA_
#define _TYPE_ integer
#define _PROC2_ distribute_I2
#endif

#ifdef _SINGLE_DATA_
#undef _SINGLE_DATA_
#define _TYPE_ real
#define _PROC2_ distribute_S2
#endif

#ifdef _DOUBLE_DATA_
#undef _DOUBLE_DATA_
#define _TYPE_ real(kind(1.0d0))
#define _PROC2_ distribute_D2
#endif

#ifdef _LOGICAL_DATA_
#undef _LOGICAL_DATA_
#define _TYPE_ logical
#define _PROC2_ distribute_L2
#endif

  subroutine _PROC2_ (vout, vin, bsize)

    _TYPE_,  intent(out) :: vout(:,:)
    _TYPE_,  intent(in)  :: vin(:,:)
    integer, intent(in), optional :: bsize(:)

    integer :: k

    ASSERT( size(vout,1) == size(vin,1) )
    if (present(bsize)) then
      ASSERT( size(bsize) == nPE )
      ASSERT( size(vout,2) == bsize(this_PE))
      if (is_IOP) then
        ASSERT( size(vin,2)  == sum(bsize))
      end if
    end if

    !! This is really inefficient; PGSLib ought to handle this directly.
    do k = 1, size(vout,1)
      call distribute (vout(k,:), vin(k,:), bsize)
    end do

  end subroutine _PROC2_

#undef _TYPE_
#undef _PROC2_
