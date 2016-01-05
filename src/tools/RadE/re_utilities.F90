!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module re_utilities

  use scl
  implicit none
  private
  
  public :: re_info, re_halt
  
contains

  subroutine re_info (mesg)
    character(len=*), intent(in) :: mesg
    if (scl_rank()==1) write(*,'(a)') mesg
  end subroutine re_info

  subroutine re_halt (errmsg)
    character(len=*), intent(in) :: errmsg
    if (scl_rank()==1) write(*,'(2a)') 'FATAL: ', errmsg
    call scl_finalize ()
    stop
  end subroutine re_halt

end module re_utilities
