! -*- Mode: f90 -*-

!*********************************************************************
! $Id: WriteString.F,v 1.3 2002/12/03 04:45:47 turner Exp $
! 
! <para>
! Write a string to a specified logical unit.
! </para>
! 
! <PARAMETER LIST>
! 
!  Input:
!   lu - logical unit to which to write
!   string - string to write out
!
!*********************************************************************
module Ubik_WriteString_module
  implicit none

  interface Spew
     module procedure WriteString
  end interface

  private
  public :: Spew

contains

  subroutine WriteString (lu, string)
    use ubik_types_module
    use Ubik_FlushUnit_module
    use Ubik_ParallelInfo_module, only: on_IO_PE

    ! <arguments>
    integer(int_type), intent(in) :: lu
    character(len=*), intent(in) :: string
    ! </arguments>
    integer(int_type) :: status

    ! only write out string if this PE is the IO PE
    status = 0
    if (on_IO_PE()) then
       write(lu,'(a)') TRIM(string)
       call FlushUnit (lu, status)
    end if

    return

  end subroutine WriteString

end module Ubik_WriteString_module
