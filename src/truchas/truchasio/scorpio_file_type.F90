!!
!! SCORPIO_FILE_TYPE
!!
!! This module defines an objected-oriented interface for doing parallel HDF5
!! output via the Scorpio library.
!!
!! Ondrej Certik <certik@lanl.gov>
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module scorpio_file_type

  use,intrinsic :: iso_fortran_env, only: int8, int32, int64, real64
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_associated, c_null_char
  use scorpio_c_binding
  implicit none
  private

  type, public :: scorpio_file
    integer(c_int) :: fhandle = -1
    type(c_ptr) :: myIOgroup = c_null_ptr
  contains
    procedure :: open_file
    procedure :: reopen_file
    procedure :: close_file

    procedure :: create_group
    procedure :: close_group
    procedure :: create_link
    procedure, private :: write_attr_int32
    procedure, private :: write_attr_real64
    procedure, private :: write_attr_real64_r1
    procedure, private :: write_attr_string
    generic :: write_attr => &
        write_attr_int32, &
        write_attr_real64, &
        write_attr_real64_r1, &
        write_attr_string

    procedure, private :: write_dataset_char_r0
    procedure, private :: write_dataset_int8_r1
    procedure, private :: write_dataset_int8_r2
    procedure, private :: write_dataset_int32_r1
    procedure, private :: write_dataset_int32_r2
    procedure, private :: write_dataset_real64_r1
    procedure, private :: write_dataset_real64_r2
    generic :: write_dataset => write_dataset_int8_r1, write_dataset_int8_r2, &
                                write_dataset_int32_r1, write_dataset_int32_r2, &
                                write_dataset_real64_r1, write_dataset_real64_r2, &
                                write_dataset_char_r0
  end type scorpio_file

contains

  subroutine open_file(this, filename, io_group_size)
    class(scorpio_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: io_group_size
    call scorpio_open_file_ext(filename//c_null_char, io_group_size, this%fhandle, this%myIOgroup)
  end subroutine open_file

  subroutine reopen_file(this, filename, io_group_size)
    class(scorpio_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: io_group_size
    call scorpio_reopen_file_ext(filename//c_null_char, io_group_size, this%fhandle, this%myIOgroup)
  end subroutine reopen_file

  subroutine close_file(this)
    class(scorpio_file), intent(inout) :: this
    if (c_associated(this%myIOgroup)) then
      call scorpio_close_file_ext(this%fhandle, this%myIOgroup)
      this%myIOgroup = c_null_ptr
      this%fhandle = -1
    end if
  end subroutine close_file

  integer(int64) function create_group(this, name) result(groupid)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: name
    groupid = scorpio_create_dataset_group(name//c_null_char, this%fhandle, this%myIOgroup)
  end function create_group

  subroutine close_group(this, groupid)
    class(scorpio_file), intent(in) :: this
    integer(int64), intent(in) :: groupid
    call scorpio_close_dataset_group(groupid, this%fhandle, this%myIOgroup)
  end subroutine close_group

  subroutine create_link(this, link_target, link_loc_id, link_name)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: link_target
    integer(int64), intent(in) :: link_loc_id
    character(*), intent(in) :: link_name
    call scorpio_create_link(link_target//c_null_char, link_loc_id, link_name//c_null_char, &
                             this%fhandle, this%myIOgroup)
  end subroutine create_link

  subroutine write_dataset_char_r0(this, name, ldata, glen)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: ldata
    integer(int32), intent(in) :: glen
    call scorpio_write_dataset2_char(&
        ldata, 1, [glen], [len(ldata)], this%fhandle, name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_dataset_int8_r1(this, name, ldata, glen)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int8), intent(in) :: ldata(:)
    integer(int32), intent(in) :: glen
    call scorpio_write_dataset2_byte(&
        ldata, 1, [glen], shape(ldata), this%fhandle, name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_dataset_int8_r2(this, name, ldata, glen)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int8), intent(in) :: ldata(:,:)
    integer(int32), intent(in) :: glen
    integer(int32) :: gdims(2), ldims(2)
    ldims(1) = size(ldata,dim=2)
    ldims(2) = size(ldata,dim=1)
    gdims(1) = glen
    gdims(2) = ldims(2)
    call scorpio_write_dataset2_byte(&
        ldata, 2, gdims, ldims, this%fhandle, name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_dataset_int32_r1(this, name, ldata, glen)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: ldata(:)
    integer(int32), intent(in) :: glen
    call scorpio_write_dataset2_int(&
        ldata, 1, [glen], shape(ldata), this%fhandle, name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_dataset_int32_r2(this, name, ldata, glen)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: ldata(:,:)
    integer(int32), intent(in) :: glen
    integer(int32) :: gdims(2), ldims(2)
    ldims(1) = size(ldata,dim=2)
    ldims(2) = size(ldata,dim=1)
    gdims(1) = glen
    gdims(2) = ldims(2)
    call scorpio_write_dataset2_int(&
        ldata, 2, gdims, ldims, this%fhandle, name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_dataset_real64_r1(this, name, ldata, glen)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:)
    integer(int32), intent(in) :: glen
    call scorpio_write_dataset2_double(&
        ldata, 1, [glen], shape(ldata), this%fhandle, name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_dataset_real64_r2(this, name, ldata, glen)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: ldata(:,:)
    integer(int32), intent(in) :: glen
    integer(int32) :: gdims(2), ldims(2)
    ldims(1) = size(ldata,dim=2)
    ldims(2) = size(ldata,dim=1)
    gdims(1) = glen
    gdims(2) = ldims(2)
    call scorpio_write_dataset2_double(&
        ldata, 2, gdims, ldims, this%fhandle, name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_attr_int32(this, obj_name, attr_name, attr_data)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: obj_name, attr_name
    integer(int32), intent(in) :: attr_data
    call scorpio_write_simple_attr_int(attr_name//c_null_char, attr_data, &
        this%fhandle, obj_name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_attr_real64(this, obj_name, attr_name, attr_data)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: obj_name, attr_name
    real(real64), intent(in) :: attr_data
    call scorpio_write_simple_attr_double(attr_name//c_null_char, attr_data, &
        this%fhandle, obj_name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_attr_string(this, obj_name, attr_name, attr_data)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: obj_name, attr_name
    character(*), intent(in) :: attr_data
    call scorpio_write_simple_attr_string(attr_name//c_null_char, attr_data//c_null_char, &
        this%fhandle, obj_name//c_null_char, this%myIOgroup)
  end subroutine

  subroutine write_attr_real64_r1(this, obj_name, attr_name, attr_data)
    class(scorpio_file), intent(in) :: this
    character(*), intent(in) :: obj_name, attr_name
    real(real64), intent(in) :: attr_data(:)
    call scorpio_write_attr_double(attr_name//c_null_char, attr_data, 1, shape(attr_data), &
        this%fhandle, obj_name//c_null_char, this%myIOgroup)
  end subroutine

end module scorpio_file_type
