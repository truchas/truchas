!!
!! TOOLPATH_DRIVER
!!
!! This module holds the optional TOOLPATH_FACTORY object and provides several
!! driver procedures that act upon the object. Conceptually this is part of the
!! Truchas multiphysics driver and functions as an adapter to the OO design of
!! the toolpath component.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module toolpath_driver

  use toolpath_type
  use toolpath_factory_type
  use parameter_list_type
  implicit none
  private

  public :: read_toolpath_namelists, known_toolpath, alloc_toolpath

  type(toolpath_factory), public :: tp_fac

contains

  subroutine read_toolpath_namelists(lun)
    use toolpath_namelist, only: read_namelists => read_toolpath_namelists
    integer, intent(in) :: lun
    call read_namelists(lun, tp_fac)
  end subroutine

  logical function known_toolpath(name)
    character(*), intent(in) :: name
    known_toolpath = tp_fac%known_toolpath(name)
  end function

  subroutine alloc_toolpath(path, name, stat, errmsg)

    type(toolpath), allocatable, intent(out) :: path
    character(*), intent(in) ::  name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist

    call tp_fac%alloc_toolpath(path, name, stat, errmsg)
    if (stat /= 0) return

    !! Write a plotfile of the toolpath if requested.
    plist => tp_fac%toolpath_plist(name)
    if (plist%is_parameter('plotfile')) then
      block
        use,intrinsic :: iso_fortran_env, only: r8 => real64
        use parallel_communication, only: is_IOP, broadcast
        use string_utilities, only: i_to_c
        integer :: lun, ios
        real(r8) :: dt
        character(:), allocatable :: plotfile
        call plist%get('plotfile', plotfile, stat, errmsg)
        if (stat /= 0) return
        call plist%get('plotfile-dt', dt, stat, errmsg)
        if (stat /= 0) return
        if (dt <= 0) then
          stat = 1
          errmsg = 'plotfile-dt must be > 0.0'
          return
        end if
        if (is_IOP) then
          open(newunit=lun,file=plotfile,action='write',status='replace',iostat=ios)
          if (ios == 0) then
            call path%write_plotfile(lun, dt)
            close(lun)
          endif
        endif
        call broadcast(ios)
        if (ios /= 0) then
          stat = 1
          errmsg = 'error opening ' // plotfile // ': iostat=' // i_to_c(ios)
          return
        end if
      end block
    end if

  end subroutine alloc_toolpath

end module toolpath_driver
