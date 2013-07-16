! -*- Mode: f90 -*-

#ifdef aix
#define FLUSH FLUSH_
#endif
!*********************************************************************
! $Id: FlushUnit.F,v 1.5 2002/12/03 05:33:21 turner Exp $
!
! <para>
! Flush the contents of the buffer for logical unit lu.
! </para>
!
! <PARAMETER LIST>
!
!  Input:
!   lu - logical unit to flush
!
!  Output:
!   status - return status
!      0  ==>  success
!     <0  ==>  internal error
!
!*********************************************************************
module Ubik_FlushUnit_module

  implicit none

  private
  public :: FlushUnit

contains

  subroutine FlushUnit (lu, status)
    use ubik_types_module

    ! <arguments>
    integer(int_type), intent(in) :: lu
    integer(int_type), intent(out) :: status
    ! </arguments>

    status = 0
    flush lu

    return
  end subroutine FlushUnit

end module Ubik_FlushUnit_module
