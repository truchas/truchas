!!
!! HL_HDF5
!!
!! A high-level layer over HDF5 which provides simplified procedures for
!! writing attributes, datasets, creating extendable datasets (of a specific
!! form) and incrementally appending to them. The provided procedures are
!! primarily limited to the needs of VTKHDF_FILE_TYPE. This depends on, and
!! supplements, the procedures from HDF5_C_BINDING which provides a limited
!! custom low-level binding to HDF5's C interface.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! March 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module hl_hdf5

  use,intrinsic :: iso_fortran_env
  use hdf5_c_binding
  implicit none
  private

  public :: h5_write_attr, h5_write_dataset, h5_create_unlimited_dataset, h5_append_to_dataset, h5_write_dataset_element

  interface h5_write_attr
    procedure h5_write_attr_int8
    procedure h5_write_attr_int32
    procedure H5_write_attr_real64
    procedure H5_write_attr_string
  end interface

  interface h5_write_dataset
    procedure write_dataset_int8
    procedure write_dataset_int32
    procedure write_dataset_real64
  end interface

  interface h5_create_unlimited_dataset
    procedure create_unlimited_dataset_int32
    procedure create_unlimited_dataset_real64
  end interface

  interface h5_append_to_dataset
    procedure append_to_dataset_int32
    procedure append_to_dataset_real64
  end interface

  interface h5_write_dataset_element
    procedure write_dataset_element_int32
  end interface

contains

  subroutine h5_write_attr_int8(obj_id, attr_name, attr_value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: attr_name
    integer(int8), intent(in) :: attr_value(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5T_STD_U8LE

    call get_attr_id(obj_id, attr_name, type_id, shape(attr_value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

    select rank (attr_value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [attr_value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, attr_value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, attr_value)
    rank default
      stat = 1
      errmsg = 'unsupported attribute value rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)
    !INSIST(istat == 0)

    if (stat /= 0) then
      errmsg = 'error writing attribute'
      return
    end if

  end subroutine h5_write_attr_int8


  subroutine h5_write_attr_int32(obj_id, attr_name, attr_value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: attr_name
    integer(int32), intent(in) :: attr_value(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5T_NATIVE_INTEGER

    call get_attr_id(obj_id, attr_name, type_id, shape(attr_value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

    select rank (attr_value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [attr_value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, attr_value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, attr_value)
    rank default
      stat = 1
      errmsg = 'unsupported attribute value rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute'
      return
    end if

  end subroutine h5_write_attr_int32


  subroutine h5_write_attr_real64(obj_id, attr_name, attr_value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: attr_name
    real(real64), intent(in) :: attr_value(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5T_NATIVE_DOUBLE

    call get_attr_id(obj_id, attr_name, type_id, shape(attr_value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

    select rank (attr_value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [attr_value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, attr_value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, attr_value)
    rank default
      stat = 1
      errmsg = 'unsupported attribute value rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute'
      return
    end if

  end subroutine h5_write_attr_real64


  subroutine h5_write_attr_string(obj_id, attr_name, attr_value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: attr_name
#ifdef INTEL_BUG20240327
    character(*), intent(in) :: attr_value
#else
    character(*), intent(in) :: attr_value(..)
#endif
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5Tcopy(H5T_NATIVE_CHARACTER)
    INSIST(type_id > 0)
    istat = H5Tset_size(type_id, len(attr_value, kind=c_size_t))
    istat = min(0, stat) ! >= 0 from H5Tset_size is successful
    INSIST(istat == 0)

    call get_attr_id(obj_id, attr_name, type_id, shape(attr_value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

#ifdef INTEL_BUG20240327
    stat = H5Awrite(attr_id, type_id, [attr_value])
#else
    select rank (attr_value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [attr_value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, attr_value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, attr_value)
    rank default
      stat = 1
      errmsg = 'unsupported attribute value rank'
      return
    end select
#endif
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)
    istat = H5Tclose(type_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute'
      return
    end if

  end subroutine h5_write_attr_string


  subroutine get_attr_id(obj_id, attr_name, type_id, dims, attr_id, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: attr_name
    integer(hid_t), intent(in) :: type_id
    integer(hsize_t), intent(in) :: dims(:) ! the Fortran shape
    integer(hid_t), intent(out) :: attr_id
    integer, intent(out) :: stat  !TODO: do we need this too?
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: space_id
    integer :: istat

    if (H5Aexists(obj_id, attr_name)) then  ! open the attribute

      attr_id = H5Aopen(obj_id, attr_name)
      if (attr_id < 0) then
        stat = 1
        errmsg = 'unable to open attribute'
        return
      end if

      !NB: we assume type_id and dims match the attribute

      stat = 0

    else ! create the attribute

      select case (size(dims))
      case (0)
        space_id = H5Screate()
      case (1)
        space_id = H5Screate(dims)
      case (2:)
        block
          integer(hsize_t) :: c_dims(size(dims))
          c_dims = dims(size(dims):1:-1) ! reverse order for the C shape
          space_id = H5Screate(c_dims)
        end block
      end select
      INSIST(space_id > 0)

      attr_id = H5Acreate(obj_id, attr_name, type_id, space_id)
      if (attr_id < 0) then
        stat = 1
        errmsg = 'unable to create attribute'
        return
      end if

      istat = H5Sclose(space_id)

      stat = 0

    end if

  end subroutine get_attr_id


  subroutine write_dataset_int8(loc_id, name, data, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int8), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, dset_id
    integer :: istat  ! ignored returned status value

    type_id = H5T_STD_U8LE

    call create_dataset(loc_id, name, type_id, shape(data, kind=hsize_t), dset_id, errmsg)
    if (dset_id < 0) then
      stat = 1
      return
    end if

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data])
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data)
    rank default
      stat = 1
      errmsg = 'unsupported dataset rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful
    istat = H5Dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing dataset "' // name // '"'
      return
    end if

  end subroutine write_dataset_int8


  subroutine write_dataset_int32(loc_id, name, data, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer, intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, dset_id
    integer :: istat  ! ignored returned status value

    type_id = H5T_NATIVE_INTEGER

    call create_dataset(loc_id, name, type_id, shape(data, kind=hsize_t), dset_id, errmsg)
    if (dset_id < 0) then
      stat = 1
      return
    end if

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data])
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data)
    rank default
      stat = 1
      errmsg = 'unsupported dataset rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful
    istat = H5Dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing dataset "' // name // '"'
      return
    end if

  end subroutine write_dataset_int32


  subroutine write_dataset_real64(loc_id, name, data, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, dset_id
    integer :: istat  ! ignored returned status value

    type_id = H5T_NATIVE_DOUBLE

    call create_dataset(loc_id, name, type_id, shape(data, kind=hsize_t), dset_id, errmsg)
    if (dset_id < 0) then
      stat = 1
      return
    end if

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data])
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data)
    rank default
      stat = 1
      errmsg = 'unsupported dataset rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful
    istat = H5Dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing dataset "' // name // '"'
      return
    end if

  end subroutine write_dataset_real64


  subroutine create_dataset(loc_id, name, type_id, dims, dset_id, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in) :: type_id
    integer(hsize_t), intent(in) :: dims(:)
    integer(hid_t), intent(out) :: dset_id
    character(:), allocatable :: errmsg

    integer(hid_t) :: space_id
    integer :: istat  ! ignored returned status value

    select case (size(dims))
    case (0)
      space_id = H5Screate()
    case (1)
      space_id = H5Screate(dims)
    case (2:)
      block
        integer(hsize_t) :: c_dims(size(dims))
        c_dims = dims(size(dims):1:-1) ! reverse order for the C shape
        space_id = H5Screate(c_dims)
      end block
    end select
    INSIST(space_id > 0)

    dset_id = H5Dcreate(loc_id, name, type_id, space_id)
    istat = H5Sclose(space_id)
    if (dset_id < 0) then
      errmsg = 'error creating dataset "' // name // '"'
      return
    end if

  end subroutine create_dataset


  subroutine create_unlimited_dataset_int32(loc_id, name, mold, chunk_size, stat, errmsg)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call create_unlimited_dataset_aux(loc_id, name, H5T_NATIVE_INTEGER, shape(mold,hsize_t), chunk_size, stat, errmsg)
  end subroutine

  subroutine create_unlimited_dataset_real64(loc_id, name, mold, chunk_size, stat, errmsg)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call create_unlimited_dataset_aux(loc_id, name, H5T_NATIVE_DOUBLE, shape(mold,hsize_t), chunk_size, stat, errmsg)
  end subroutine


  subroutine create_unlimited_dataset_aux(loc_id, name, type_id, dims, chunk_size, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in) :: type_id
    integer(hsize_t), intent(in) :: dims(:)
    integer, intent(in) :: chunk_size
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: istat
    integer(hid_t) :: space_id, dcpl_id, dset_id
    integer(hsize_t), allocatable :: c_dims(:), maxdims(:), chunk_dims(:)

    if (size(dims) == 0) then
      stat = 1
      errmsg = 'invalid rank for unlimited dataset'
      return
    end if

    c_dims = dims
    if (size(dims) > 1) c_dims = dims(size(dims):1:-1)

    maxdims = c_dims
    maxdims(1) = H5S_UNLIMITED
    chunk_dims = c_dims
    chunk_dims(1) = chunk_size
    c_dims(1) = 0

    space_id = H5Screate(c_dims, maxdims)
    INSIST(space_id > 0)

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE)
    INSIST(dcpl_id > 0)

    stat = H5Pset_chunk(dcpl_id, size(chunk_dims), chunk_dims)
    stat = min(0, stat) ! >= 0 from H5Pset_chunk is successful
    INSIST(stat == 0)

    dset_id = H5Dcreate(loc_id, name, type_id, space_id, dcpl_id=dcpl_id)
    if (dset_id < 0) then
      stat = 1
      errmsg = 'error creating unlimited dataset "' // name // '"'
    end if

    istat = h5sclose(space_id)
    istat = h5pclose(dcpl_id)
    istat = h5dclose(dset_id)

  end subroutine create_unlimited_dataset_aux


  subroutine append_to_dataset_aux(loc_id, name, dims, dset_id, mem_space_id, data_space_id, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hsize_t), intent(in) :: dims(:)
    integer(hid_t), intent(out) :: dset_id, data_space_id, mem_space_id
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hsize_t), allocatable :: mem_dims(:), data_dims(:), maxdims(:), start(:)
    integer :: ndims, istat

    dset_id = H5Dopen(loc_id, name)
    if (dset_id < 0) then
      stat = 1
      errmsg = 'unable to open dataset "' // name // '"'
      return
    end if

    if (size(dims) == 0) then ! scalar appended data
      mem_dims = [1]  ! compatible with a rank-1 extendable dataset
    else
      mem_dims = dims
      if (size(dims) > 1) mem_dims = dims(size(dims):1:-1) ! reverse order for the C shape
    end if

    mem_space_id = H5Screate(mem_dims)
    INSIST(mem_space_id > 0)

    data_space_id = H5Dget_space(dset_id)
    INSIST(data_space_id > 0)

    ndims = H5Sget_simple_extent_ndims(data_space_id)
    INSIST(ndims > 0)

    if (size(mem_dims) /= ndims) then
      stat = 1
      errmsg = 'rank of appended data is incompatible with dataset "' // name // '"'
      return
    end if

    allocate(data_dims, maxdims, start, mold=mem_dims)

    ndims = H5Sget_simple_extent_dims(data_space_id, data_dims, maxdims)
    INSIST(ndims > 0)
    istat = H5Sclose(data_space_id)

    if (any(mem_dims(2:) /= data_dims(2:))) then
      stat = 1
      errmsg = 'shape of appended data is incompatible with dataset "' // name // '"'
      return
    end if

    ! Resize the dataset to accomodate the appended data
    start = 0
    start(1) = data_dims(1)
    data_dims(1) = data_dims(1) + mem_dims(1)
    stat = H5Dset_extent(dset_id, data_dims)
    stat = min(0, stat) ! >= 0 from H5Dset_extent is successful
    if (stat /= 0) then
      errmsg = 'error resizing dataset "' // name // '"'
      return
    end if

    ! Create space corresponding to appended elements in dataset
    data_space_id = H5Dget_space(dset_id)
    stat = H5Sselect_hyperslab(data_space_id, H5S_SELECT_SET, start, mem_dims)
    stat = min(0, stat) ! >= 0 from H5Sselect_hyperslab is successful
    if (stat /= 0) then
      stat = 1
      errmsg = 'error creating hyperslab of dataset "' // name // '"'
      return
    end if

  end subroutine append_to_dataset_aux


  subroutine append_to_dataset_int32(loc_id, name, data, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id
    integer :: istat

    type_id = H5T_NATIVE_INTEGER

    call append_to_dataset_aux(loc_id, name, shape(data, hsize_t), dset_id, mem_space_id, data_space_id, stat, errmsg)
    if (stat /= 0) return

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data], mem_space_id, data_space_id)
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id)
    rank default
      stat = 1
      errmsg = 'unsupported rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful

    istat = h5sclose(mem_space_id)
    istat = h5sclose(data_space_id)
    istat = h5dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing to dataset "' // name // '"'
      return
    end if

  end subroutine append_to_dataset_int32


  subroutine append_to_dataset_real64(loc_id, name, data, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id
    integer :: istat

    type_id = H5T_NATIVE_DOUBLE

    call append_to_dataset_aux(loc_id, name, shape(data, hsize_t), dset_id, mem_space_id, data_space_id, stat, errmsg)
    if (stat /= 0) return

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data], mem_space_id, data_space_id)
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id)
    rank default
      stat = 1
      errmsg = 'unsupported rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful

    istat = h5sclose(mem_space_id)
    istat = h5sclose(data_space_id)
    istat = h5dclose(dset_id)

  end subroutine append_to_dataset_real64


  subroutine write_dataset_element_int32(loc_id, name, elem, val, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer, intent(in) :: elem
    integer(int32), intent(in) :: val
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id
    integer(hsize_t) :: coord(1)
    integer :: istat

    type_id = H5T_NATIVE_INTEGER

    dset_id = H5Dopen(loc_id, name)
    if (dset_id < 0) then
      stat = 1
      errmsg = 'unable to open dataset "' // name // '"'
      return
    end if

    mem_space_id = H5Screate()
    INSIST(mem_space_id > 0)

    data_space_id = H5Dget_space(dset_id)
    INSIST(data_space_id > 0)

    coord = [elem]
    stat = H5Sselect_elements(data_space_id, H5S_SELECT_SET, 1_c_size_t, coord)
    stat = min(0, stat) ! >= 0 from H5Sselect_elements is successful
    INSIST(stat == 0)

    stat = H5Dwrite(dset_id, type_id, [val], mem_space_id, data_space_id)
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful

    istat = h5sclose(mem_space_id)
    istat = h5sclose(data_space_id)
    istat = h5dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing to dataset "' // name // '"'
      return
    end if

  end subroutine write_dataset_element_int32

end module hl_hdf5
