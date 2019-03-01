!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flow_input_utils
  use kinds, only: r8
  use input_utilities
  use parameter_list_type
  implicit none
  private

  public :: plist_set_if
  ! re-export input_utilities
  public :: seek_to_namelist, null_i, null_r, null_c

  interface plist_set_if
    module procedure set_i0, set_i1, set_r0, set_r1, set_c, set_l
  end interface plist_set_if

contains

  subroutine set_i0(p, name, v)
    type(parameter_list), pointer, intent(inout) :: p
    character(*), intent(in) :: name
    integer, intent(in) :: v
    if (v /= null_i) call p%set(name, v)
  end subroutine set_i0

  subroutine set_i1(p, name, v)
    type(parameter_list), pointer, intent(inout) :: p
    character(*), intent(in) :: name
    integer, intent(in) :: v(:)
    integer :: i, j

    j = -1
    do i = 1, size(v)
      if (v(i) == null_i) then
        j = i-1
        exit
      end if
    end do

    if (j == -1) then ! all value are valid
      call p%set(name, v)
    else if (j > 0) then ! at least on value is valid
      call p%set(name, v(1:j))
    end if
  end subroutine set_i1

  subroutine set_r0(p, name, v)
    type(parameter_list), pointer, intent(inout) :: p
    character(*), intent(in) :: name
    real(r8), intent(in) :: v
    if (v /= null_r) call p%set(name, v)
  end subroutine set_r0

  subroutine set_r1(p, name, v)
    type(parameter_list), pointer, intent(inout) :: p
    character(*), intent(in) :: name
    real(r8), intent(in) :: v(:)
    integer :: i, j

    j = -1
    do i = 1, size(v)
      if (v(i) == null_r) then
        j = i-1
        exit
      end if
    end do

    if (j == -1) then ! all value are valid
      call p%set(name, v)
    else if (j > 0) then ! at least on value is valid
      call p%set(name, v(1:j))
    end if
  end subroutine set_r1

  subroutine set_c(p, name, v)
    type(parameter_list), pointer, intent(inout) :: p
    character(*), intent(in) :: name
    character(*), intent(in) :: v
    if (len_trim(v) == 0) return
    if (v /= null_c) call p%set(name, v)
  end subroutine set_c

  subroutine set_l(p, name, v)
    type(parameter_list), pointer, intent(inout) :: p
    character(*), intent(in) :: name
    logical, intent(in) :: v
    if (v) call p%set(name, v)
  end subroutine set_l


end module flow_input_utils
