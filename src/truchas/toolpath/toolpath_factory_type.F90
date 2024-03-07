!!
!! TOOLPATH_FACTORY_TYPE
!!
!! This module defines the derived type TOOLPATH_FACTORY, which is a high
!! level wrapper around the fundamental ALLOC_TOOLPATH procedure provided by
!! the TOOLPATH_FACTORY module. A TOOLPATH_FACTORY object holds a parameter
!! list that consists of a collection of named sublists, each a parameter list
!! that defines a toolpath. Client code can then use a factory object to
!! instantiate a particular toolpath by referencing its name.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module toolpath_factory_type

  use toolpath_type
  use parameter_list_type
  implicit none
  private

  type, public :: toolpath_factory
    private
    type(parameter_list) :: params
  contains
    generic :: alloc_toolpath => alloc_by_name, alloc_by_plist
    procedure, private :: alloc_by_name
    procedure, nopass, private :: alloc_by_plist
    procedure :: toolpath_plist
    procedure :: known_toolpath
  end type

contains

  function toolpath_plist(this, name) result(plist)
    class(toolpath_factory), intent(inout) :: this
    character(*), intent(in) :: name
    type(parameter_list), pointer :: plist
    plist => this%params%sublist(name)
  end function

  logical function known_toolpath(this, name)
    class(toolpath_factory), intent(in) :: this
    character(*), intent(in) :: name
    known_toolpath = this%params%is_sublist(name)
  end function

  subroutine alloc_by_name(this, path, name, stat, errmsg)

    class(toolpath_factory), intent(inout) :: this
    type(toolpath), allocatable, intent(out) :: path
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist

    if (this%params%is_sublist(name)) then
      plist => this%params%sublist(name)
      call alloc_by_plist(path, plist, stat, errmsg)
    else
      stat = 1
      errmsg = 'unknown toolpath: ' // name
    end if

  end subroutine

  subroutine alloc_by_plist(path, plist, stat, errmsg)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use toolpath_type
    use toolpath_factory, only: alloc_tp => alloc_toolpath

    type(toolpath), allocatable, intent(out) :: path
    type(parameter_list), intent(inout) :: plist
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: ds

    call alloc_tp(path, plist, stat, errmsg)
    if (stat /= 0) return

    !! Add path partition data if requested.
    if (plist%is_parameter('partition-ds')) then
      call plist%get('partition-ds', ds, stat, errmsg)
      if (stat /= 0) return
      call path%set_partition(ds)
    end if

  end subroutine alloc_by_plist

end module toolpath_factory_type
