!!
!! SCORPIO_C_BINDING
!!
!! Bindings to a subset of the Scorpio C library and additions.
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

module scorpio_c_binding

  use,intrinsic :: iso_c_binding, only: c_ptr, c_char, c_double, c_int, c_int8_t, c_int64_t
  implicit none
  private

  public :: scorpio_open_file_ext
  public :: scorpio_close_file_ext
  public :: scorpio_reopen_file_ext
  public :: scorpio_create_dataset_group
  public :: scorpio_close_dataset_group
  public :: scorpio_create_link
  public :: scorpio_write_dataset2_char
  public :: scorpio_write_dataset2_byte
  public :: scorpio_write_dataset2_int
  public :: scorpio_write_dataset2_double
  public :: scorpio_write_simple_attr_int
  public :: scorpio_write_simple_attr_double
  public :: scorpio_write_attr_double
  public :: scorpio_write_simple_attr_string

  !! SCORPIO library functions
  interface
    function scorpio_create_dataset_group(group_name, fhandle, myIOgroup) result(groupid) bind(c)
      import c_int, c_char, c_ptr, c_int64_t
      character(kind=c_char), intent(in) :: group_name(*)
      integer(c_int), value :: fhandle
      type(c_ptr), value :: myIOgroup
      integer(c_int64_t) :: groupid
    end function
    !NB: really function with a non-collective return code
    subroutine scorpio_close_dataset_group(groupid, fhandle, myIOgroup) bind(c)
      import c_int, c_ptr, c_int64_t
      integer(c_int64_t), value :: groupid
      integer(c_int), value :: fhandle
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_create_link(target, link_loc_id, link_name, fhandle, myIOgroup) bind(c)
      import c_int, c_char, c_ptr, c_int64_t
      character(kind=c_char) :: target(*), link_name(*)
      integer(c_int64_t), value :: link_loc_id
      integer(c_int), value :: fhandle
      type(c_ptr), value :: myIOgroup
    end subroutine
  end interface

  !! SCORPIO_EXT -- our custom extensions to Scorpio
  interface
    subroutine scorpio_open_file_ext(filename, groupSize, fhandle, myIOgroup) bind(c)
      import c_char, c_int, c_ptr
      character(kind=c_char), intent(in) :: filename(*)
      integer(c_int), value :: groupSize
      integer(c_int), intent(out) :: fhandle
      type(c_ptr), intent(out) :: myIOgroup
    end subroutine
    subroutine scorpio_reopen_file_ext(filename, groupSize, fhandle, myIOgroup) bind(c)
      import c_char, c_int, c_ptr
      character(kind=c_char), intent(in) :: filename(*)
      integer(c_int), value :: groupSize
      integer(c_int), intent(out) :: fhandle
      type(c_ptr), intent(out) :: myIOgroup
    end subroutine
    subroutine scorpio_close_file_ext(fhandle, myIOgroup) bind(c)
      import c_int, c_ptr
      integer(c_int), value :: fhandle
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_write_dataset2_char(vector, ndims, globaldims, localdims, &
        fhandle, dset_name, myIOgroup) bind(c)
      import c_int, c_char, c_ptr
      character(kind=c_char), intent(in) :: vector(*)
      integer(c_int), value :: ndims, fhandle
      integer(c_int), intent(in) :: globaldims(*), localdims(*)
      character(kind=c_char), intent(in) :: dset_name(*)
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_write_dataset2_byte(vector, ndims, globaldims, localdims, &
        fhandle, dset_name, myIOgroup) bind(c)
      import c_int, c_int8_t, c_char, c_ptr
      integer(c_int8_t), intent(in) :: vector(*)
      integer(c_int), value :: ndims, fhandle
      integer(c_int), intent(in) :: globaldims(*), localdims(*)
      character(kind=c_char), intent(in) :: dset_name(*)
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_write_dataset2_int(vector, ndims, globaldims, localdims, &
        fhandle, dset_name, myIOgroup) bind(c)
      import c_int, c_char, c_ptr
      integer(c_int), intent(in) :: vector(*)
      integer(c_int), value :: ndims, fhandle
      integer(c_int), intent(in) :: globaldims(*), localdims(*)
      character(kind=c_char), intent(in) :: dset_name(*)
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_write_dataset2_double(vector, ndims, globaldims, localdims, &
        fhandle, dset_name, myIOgroup) bind(c)
      import c_int, c_double, c_char, c_ptr
      real(c_double), intent(in) :: vector(*)
      integer(c_int), value :: ndims, fhandle
      integer(c_int), intent(in) :: globaldims(*), localdims(*)
      character(kind=c_char), intent(in) :: dset_name(*)
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_write_simple_attr_int(attr_name, attr_data, &
        fhandle, obj_name, myIOgroup) bind(c)
      import c_int, c_char, c_ptr
      character(kind=c_char), intent(in) :: attr_name(*)
      integer(c_int), intent(in) :: attr_data
      integer(c_int), value :: fhandle
      character(kind=c_char), intent(in) :: obj_name(*)
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_write_simple_attr_double(attr_name, attr_data, &
        fhandle, obj_name, myIOgroup) bind(c)
      import c_int, c_char, c_double, c_ptr
      character(kind=c_char), intent(in) :: attr_name(*)
      real(c_double), intent(in) :: attr_data
      integer(c_int), value :: fhandle
      character(kind=c_char), intent(in) :: obj_name(*)
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_write_attr_double(attr_name, attr_data, ndims, adims, &
        fhandle, obj_name, myIOgroup) bind(c)
      import :: c_int, c_char, c_double, c_ptr
      character(kind=c_char), intent(in) :: attr_name(*)
      real(c_double), intent(in) :: attr_data(*)
      integer(c_int), value :: ndims
      integer(c_int), intent(in) :: adims(*)
      integer(c_int), value :: fhandle
      character(kind=c_char), intent(in) :: obj_name(*)
      type(c_ptr), value :: myIOgroup
    end subroutine
    subroutine scorpio_write_simple_attr_string(attr_name, attr_data, &
        fhandle, obj_name, myIOgroup) bind(c)
      import c_int, c_char, c_ptr
      character(kind=c_char), intent(in) :: attr_name(*), attr_data(*)
      integer(c_int), value :: fhandle
      character(kind=c_char), intent(in) :: obj_name(*)
      type(c_ptr), value :: myIOgroup
    end subroutine
  end interface

end module
