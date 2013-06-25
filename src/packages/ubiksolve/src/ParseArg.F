! -*- Mode: f90 -*-

!*********************************************************************
! $Id: ParseArg.F,v 1.7 2002/12/03 04:58:11 turner Exp $
!
! Routines to parse command-line arguments of the form:
!
!   name=value
!
!*********************************************************************
module Ubik_ParseArg_module
#ifdef NAG_COMPILER
  use f90_unix
#endif
  use ubik_types_module
  implicit none

  interface ParseArg
     module procedure ParseArgChar
     module procedure ParseArgFloat
     module procedure ParseArgInt
  end interface

  private
  public :: ParseArg

contains
  
  !*********************************************************************
  !
  ! <para>
  ! Read a character value from a command line argument of the form:
  !
  !      name=value
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   name - character string used on command line
  !
  !  Output
  !   value - value of the variable
  !
  !*********************************************************************
  subroutine ParseArgChar (name, value)

    ! <arguments>
    character(len=*), intent(in) :: name
    character(len=*), intent(inout) :: value
    ! </arguments>

    ! Local:
    character(len=40) :: arg
    integer(int_type) :: i, ipos
#ifndef NAG_COMPILER
    integer(int_type) :: IARGC
#endif

    do i=1,IARGC()
       call GETARG (i, arg)
       ipos = INDEX(arg, '=')
       if (arg(1:ipos-1) == name) then
#if defined(aix) || defined(osf1)
          read(arg(ipos+1:),*) value
#else
          read(arg(ipos+1:),'(a)') value
#endif
          return
       end if
    end do

    return
  end subroutine ParseArgChar
  
  !*********************************************************************
  !
  ! <para>
  ! Read a real value from a command line argument of the form:
  !
  !   name=value
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   name - character string used on command line
  !
  !  Output
  !   value - value of the variable
  !
  !*********************************************************************
  subroutine ParseArgFloat (name, value)

    ! <arguments>
    character(len=*), intent(in) :: name
    real(real_type), intent(inout) :: value
    ! </arguments>

    ! Local:
    character(len=40) :: arg
    integer(int_type) :: i, ipos
#ifndef NAG_COMPILER
    integer(int_type) :: IARGC
#endif

    do i=1,IARGC()
       call GETARG (i, arg)
       ipos = INDEX(arg, '=')
       if (arg(1:ipos-1) == name) then
#if defined(osf1)
          read(arg(ipos+1:),'(g15.8)',end=9999) value
#else
          read(arg(ipos+1:),*,end=9999) value
#endif
          return
       end if
    end do

9999 continue
    return
  end subroutine ParseArgFloat
  
  !*********************************************************************
  !
  ! <para>
  ! Read an integer value from a command line argument of the form:
  !
  !   name=value
  ! </para>
  !
  ! <PARAMETER LIST>
  !
  !  Input:
  !   name - character string used on command line
  !
  !  Output
  !   value - value of the variable
  !
  !*********************************************************************
  subroutine ParseArgInt (name, value)

    ! <arguments>
    character(len=*), intent(in) :: name
    integer(int_type), intent(inout) :: value
    ! </arguments>

    ! Local:
    character(len=40) :: arg
    integer(int_type) :: i, ipos
#ifndef NAG_COMPILER
    integer(int_type) :: IARGC
#endif

    do i=1,IARGC()
       call GETARG (i, arg)
       ipos = INDEX(arg, '=')
       if (arg(1:ipos-1) == name) then
#if defined(osf1)
          read(arg(ipos+1:),'(i10)') value
#else
          read(arg(ipos+1:),*) value
#endif
          return
       end if
    end do

    return
  end subroutine ParseArgInt

end module Ubik_ParseArg_module
