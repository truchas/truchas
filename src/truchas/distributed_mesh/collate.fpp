# /* collate -- include file for parallel_communication.F90 */

#ifdef _TYPE_
#undef _TYPE_
#endif

#ifdef _INT8_DATA_
#undef _INT8_DATA_
#define _TYPE_ integer(int8)
#define _PROC2_ collate_int8_2
#endif

#ifdef _INTEGER_DATA_
#undef _INTEGER_DATA_
#define _TYPE_ integer
#define _PROC2_ collate_I2
#endif

#ifdef _SINGLE_DATA_
#undef _SINGLE_DATA_
#define _TYPE_ real
#define _PROC2_ collate_S2
#endif

#ifdef _DOUBLE_DATA_
#undef _DOUBLE_DATA_
#define _TYPE_ real(kind(1.0d0))
#define _PROC2_ collate_D2
#endif

#ifdef _LOGICAL_DATA_
#undef _LOGICAL_DATA_
#define _TYPE_ logical
#define _PROC2_ collate_L2
#endif

#ifndef _TYPE_
#error "one of LOGICAL_DATA, INTEGER_DATA, SINGLE_DATA, DOUBLE_DATA must be defined"
#endif

  subroutine _PROC2_ (vout, vin)
  
    _TYPE_, intent(out) :: vout(:,:)
    _TYPE_, intent(in)  :: vin(:,:)
    
    integer :: k
    
    ASSERT( size(vout,1) == size(vin,1) )
    
    !! This is really inefficient; PGSLib ought to handle this directly.
    do k = 1, size(vout,1)
      call collate (vout(k,:), vin(k,:))
    end do
    
  end subroutine _PROC2_

#undef _TYPE_
#undef _PROC2_
