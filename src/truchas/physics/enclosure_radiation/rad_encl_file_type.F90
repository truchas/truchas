!!
!! RAD_ENCL_FILE_TYPE
!!
!! This module provides the interface for reading and writing data from and to
!! the radiation enclosure file.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! October 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module rad_encl_file_type

  use,intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use,intrinsic :: iso_c_binding, only: c_size_t
  use netcdf_c_binding
  use netcdf_file_type
  implicit none
  private

  type, public :: rad_encl_file
    private
    type(netcdf_file) :: file
    !integer(c_int) :: num_nodes = 0, num_faces = 0, num_fnode = 0, num_groups = 0
  contains
    procedure :: create
    procedure :: open_ro
    procedure :: open_rw
    procedure :: close
    procedure :: put_surface
    procedure :: get_surface
    procedure :: put_symmetry_info
    procedure :: get_symmetry_info
    procedure :: put_group_info
    procedure :: get_group_info
    procedure :: put_source_info
    procedure :: get_source_info
    procedure :: init_patch
    procedure :: get_patch_dims
    procedure :: init_vf
    procedure :: get_vf_dims
    procedure :: put_vf_rowcount
    procedure :: get_vf_rowcount
    procedure :: put_vf_rows
    procedure :: get_vf_rows
    procedure :: put_ambient
    procedure :: get_ambient
    procedure :: put_f2p_map
    procedure :: get_f2p_map
    procedure :: put_area
    procedure :: get_area
    procedure :: put_face_weight
    procedure :: get_face_weight
    procedure :: has_vf_data
    procedure :: has_ambient
    procedure :: has_patches
    procedure :: has_area
    procedure :: has_face_weight
  end type rad_encl_file

contains

  subroutine create(this, path, stat, errmsg)
    class(rad_encl_file), intent(out) :: this
    character(*), intent(in) :: path
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    call this%file%create(path, ior(NF_CLOBBER,NF_NETCDF4), stat, errmsg)
  end subroutine create

  subroutine open_ro(this, path, stat, errmsg)
    class(rad_encl_file), intent(out) :: this
    character(*), intent(in) :: path
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    call this%file%open(path, NF_NOWRITE, stat, errmsg)
  end subroutine open_ro

  subroutine open_rw(this, path, stat, errmsg)
    class(rad_encl_file), intent(out) :: this
    character(*), intent(in) :: path
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    call this%file%open(path, NF_WRITE, stat, errmsg)
  end subroutine open_rw

  subroutine close(this, stat, errmsg)
    class(rad_encl_file), intent(inout) :: this
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    call this%file%close(stat, errmsg)
  end subroutine close

  !! Write the enclosure surface to the netCDF file.
  subroutine put_surface(this, xface, fnode, coord)

    class(rad_encl_file), intent(in) :: this
    integer(int32), intent(in) :: xface(:), fnode(:)
    real(real64), intent(in) :: coord(:,:)

    integer :: stat, dim1, dim2, coord_id, xface_id, fnode_id
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%put_surface: '

    !! Put the dataset into define mode.
    call this%file%redef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Define the coordinate variable.
    call this%file%def_dim('num_dim', size(coord,dim=1,kind=c_size_t), dim1, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_dim('num_nodes', size(coord,dim=2,kind=c_size_t), dim2, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('coord', NF_REAL64, [dim1,dim2], coord_id, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Define the face connectivity variables.
    call this%file%def_dim('num_faces', size(xface,kind=c_size_t)-1, dim1, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_dim('num_xface', size(xface,kind=c_size_t), dim1, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('xface', NF_INT32, dim1, xface_id, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_dim('num_fnode', size(fnode,kind=c_size_t), dim1, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('fnode', NF_INT32, dim1, fnode_id, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Take the dataset out of define mode.
    call this%file%enddef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Write the coordinate data.
    call this%file%put_var(coord_id, coord, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Write the face connectivity data.
    call this%file%put_var(xface_id, xface, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%put_var(fnode_id, fnode, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !this%num_nodes = size(coord,dim=2)
    !this%num_faces = size(xface) - 1
    !this%num_fnode = size(fnode)

  end subroutine put_surface

  subroutine get_surface(this, fsize, fnode, coord)

    class(rad_encl_file), intent(in) :: this
    integer(int32), allocatable, intent(out) :: fsize(:), fnode(:)
    real(real64), allocatable, intent(out) :: coord(:,:)

    integer :: stat, dimid, varid
    integer(c_size_t) :: n, ndim
    integer(int32), allocatable :: xface(:)
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%get_surface: '

    !! Read the fsize data.
    call this%file%inq_dimid('num_faces', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    allocate(fsize(n), xface(n+1))
    call this%file%inq_varid('xface', varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%get_var(varid, xface, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    fsize = xface(2:) - xface(:size(fsize))
    deallocate(xface)

    !! Read the fnode data.
    call this%file%inq_dimid('num_fnode', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    allocate(fnode(n))
    call this%file%inq_varid('fnode', varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%get_var(varid, fnode, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Read the coord data.
    call this%file%inq_dimid('num_dim', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, ndim, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimid('num_nodes', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    allocate(coord(ndim,n))
    call this%file%inq_varid('coord', varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%get_var(varid, coord, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine get_surface

  subroutine put_group_info (this, gnum, group_ids)

    class(rad_encl_file), intent(in) :: this
    integer(int32), intent(in) :: gnum(:), group_ids(:)

    integer :: stat, dimid, gnum_id, group_id
    integer(c_size_t) :: n
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%put_group_info: '

    !! Get the num_faces dimension from the database.
    call this%file%inq_dimid('num_faces', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    INSIST(size(gnum,kind=c_size_t) == n)
    INSIST(minval(gnum) > 0)
    INSIST(maxval(gnum) <= size(group_ids))

    !! Put the dataset into define mode.
    call this%file%redef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Define the gnum variable.
    call this%file%def_var('gnum', NF_INT32, dimid, gnum_id, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Define the group_ids variable.
    call this%file%def_dim('num_group', size(group_ids,kind=c_size_t), dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('group_ids', NF_INT32, dimid, group_id, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Take the dataset out of define mode.
    call this%file%enddef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Write the gnum data.
    call this%file%put_var(gnum_id, gnum, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Write the group_ids data.
    call this%file%put_var(group_id, group_ids, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine put_group_info

  subroutine get_group_info(this, gnum, group_ids)

    class(rad_encl_file), intent(in) :: this
    integer(int32), allocatable, intent(out) :: gnum(:), group_ids(:)

    integer :: stat, dimid, varid
    integer(c_size_t) :: n
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%get_group_info: '

    !! Read the gnum data.
    call this%file%inq_dimid('num_faces', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    allocate(gnum(n))
    call this%file%inq_varid('gnum', varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%get_var(varid, gnum, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Read the group_ids data.
    call this%file%inq_dimid('num_group', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    allocate(group_ids(n))
    call this%file%inq_varid('group_ids', varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%get_var(varid, group_ids, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine get_group_info

  subroutine put_symmetry_info(this, mirror, rot_axis, num_rot)

    class(rad_encl_file), intent(in) :: this
    logical, intent(in) :: mirror(:)
    integer, intent(in) :: rot_axis, num_rot

    integer :: stat, dimid, varid
    integer(int32) :: symmetry(5)
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%put_symmetry_info: '

    INSIST(size(mirror) == 3)
    INSIST(rot_axis >= 0 .and. rot_axis <=3)

    !! Define the symmetry variable.
    call this%file%redef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_dim('num_symmetry', size(symmetry,kind=c_size_t), dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('symmetry', NF_INT32, dimid, varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%enddef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Write the serialized symmetry info.
    symmetry = 0
    where (mirror) symmetry(1:3) = 1
    symmetry(4) = rot_axis
    if (rot_axis > 0) then
      INSIST(num_rot > 1)
      symmetry(5) = num_rot
    end if
    call this%file%put_var(varid, symmetry, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine put_symmetry_info

  subroutine get_symmetry_info(this, mirror, rot_axis, num_rot)

    class(rad_encl_file), intent(in) :: this
    logical, intent(out) :: mirror(:)
    integer, intent(out) :: rot_axis, num_rot

    integer :: stat, varid
    integer(int32) :: symmetry(5)
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%get_symmetry_info: '

    INSIST(size(mirror) == 3)

    !! Read the serialized symmetry info.
    call this%file%inq_varid('symmetry', varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%get_var(varid, symmetry, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    !! Unpack the serialized symmetry info.
    mirror  = (symmetry(1:3) /= 0)
    rot_axis = symmetry(4)
    num_rot  = symmetry(5)

  end subroutine get_symmetry_info

  subroutine put_source_info(this, src_elem, src_side)

    class(rad_encl_file), intent(in) :: this
    integer(int32), intent(in) :: src_elem(:), src_side(:)

    integer :: stat, dimid, src_elem_id, src_side_id
    integer(c_size_t) :: n
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%put_source_info: '

    !! Get num_faces from the database.
    call this%file%inq_dimid('num_faces', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    INSIST(size(src_elem) == n)
    INSIST(size(src_side) == n)

    !! Define the src_elem and src_side variables.
    call this%file%redef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('src_elem', NF_INT32, dimid, src_elem_id, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('src_side', NF_INT32, dimid, src_side_id, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%enddef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    call this%file%put_var(src_elem_id, src_elem, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%put_var(src_side_id, src_side, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine put_source_info

  subroutine get_source_info(this, src_elem, src_side)

    class(rad_encl_file), intent(in) :: this
    integer(int32), allocatable, intent(out) :: src_elem(:), src_side(:)

    integer :: stat, dimid, varid
    integer(c_size_t) :: n
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%get_source_info: '

    call this%file%inq_dimid('num_faces', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    allocate(src_elem(n), src_side(n))
    call this%file%inq_varid('src_elem', varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%get_var(varid, src_elem, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_varid('src_side', varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%get_var(varid, src_side, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine get_source_info

  subroutine init_patch(this, npatch, write_patch_data)

    class(rad_encl_file), intent(in) :: this
    integer(c_size_t), intent(in) :: npatch
    logical, intent(in) :: write_patch_data

    integer :: stat, dimid, varid
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%init_patch: '

    call this%file%redef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    call this%file%def_dim('num_patches', npatch, dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    if (write_patch_data) then
      call this%file%inq_dimid('num_faces', dimid, stat, errmsg)
      if (stat /= 0) call error_exit(proc//errmsg)
      call this%file%def_var('f2p_map', NF_INT32, dimid, varid, stat, errmsg)
      if (stat /= 0) call error_exit(proc//errmsg)
      call this%file%def_var('face_weight', NF_REAL64, dimid, varid, stat, errmsg)
      if (stat /= 0) call error_exit(proc//errmsg)
    end if

    call this%file%enddef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine init_patch

  subroutine get_patch_dims(this, nface, npatch)

    class(rad_encl_file), intent(in) :: this
    integer, intent(out) :: nface, npatch

    integer :: stat, dimid
    integer(c_size_t) :: n
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%get_patch_dims: '

    !! Get the number of faces.
    call this%file%inq_dimid('num_faces', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, n, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    nface = n

    !! Get the number of patches
    if (this%has_patches()) then
      call this%file%inq_dimid('num_patches', dimid, stat, errmsg)
      if (stat /= 0) call error_exit(proc//errmsg)
      call this%file%inq_dimlen(dimid, n, stat, errmsg)
      if (stat /= 0) call error_exit(proc//errmsg)
      npatch = n
    else
      npatch = nface
    end if

  end subroutine get_patch_dims

  subroutine init_vf(this, nnonz, write_ambient)

    class(rad_encl_file), intent(in) :: this
    integer(c_size_t), intent(in) :: nnonz
    logical, intent(in) :: write_ambient

    integer :: stat, dimid, varid
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%init_vf: '

    call this%file%redef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    call this%file%def_dim('num_nonzero', nnonz, dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('val', NF_REAL32, dimid, varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('ja', NF_INT32, dimid, varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

    call this%file%inq_dimid('num_patches', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('icount', NF_INT32, dimid, varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%def_var('area', NF_REAL64, dimid, varid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    if (write_ambient) then
      call this%file%def_var('ambient', NF_REAL32, dimid, varid, stat, errmsg)
      if (stat /= 0) call error_exit(proc//errmsg)
    end if

    call this%file%enddef(stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine init_vf

  subroutine get_vf_dims(this, nface, npatch, nnonz)

    class(rad_encl_file), intent(in) :: this
    integer, intent(out) :: nface, npatch
    integer(c_size_t), intent(out) :: nnonz

    integer :: stat, dimid
    character(:), allocatable :: errmsg
    character(*), parameter :: proc = 'rad_encl_file%get_vf_dims: '

    !! Get the number of faces and patches.
    call this%get_patch_dims(nface, npatch)

    !! Get the number of nonzeros in the VF matrix.
    call this%file%inq_dimid('num_nonzero', dimid, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)
    call this%file%inq_dimlen(dimid, nnonz, stat, errmsg)
    if (stat /= 0) call error_exit(proc//errmsg)

  end subroutine get_vf_dims

  subroutine put_vf_rowcount(this, rowcount)
    class(rad_encl_file), intent(in) :: this
    integer(int32), intent(in) :: rowcount(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('icount', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_vf_rowcount: '//errmsg)
    call this%file%put_var(varid, rowcount)
    if (stat /= 0) call error_exit('rad_encl_file%put_vf_rowcount: '//errmsg)
  end subroutine put_vf_rowcount

  subroutine get_vf_rowcount(this, rowcount)
    class(rad_encl_file), intent(in) :: this
    integer(int32), intent(out) :: rowcount(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('icount', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_vf_rowcount: '//errmsg)
    call this%file%get_var(varid, rowcount)
    if (stat /= 0) call error_exit('rad_encl_file%get_vf_rowcount: '//errmsg)
  end subroutine get_vf_rowcount

  subroutine put_vf_rows(this, val, ja, start)
    class(rad_encl_file), intent(in) :: this
    real(real32), intent(in) :: val(:)
    integer(int32), intent(in) :: ja(:)
    integer(c_size_t), intent(in) :: start
    integer :: stat, varid
    character(:), allocatable :: errmsg
    ASSERT(size(ja) == size(val))
    if (size(ja) == 0) return
    !! Should check that we're not trying to write too much data
    call this%file%inq_varid('val', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_vf_rows: '//errmsg)
    call this%file%put_var(varid, start, val, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_vf_rows: '//errmsg)
    call this%file%inq_varid('ja', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_vf_rows: '//errmsg)
    call this%file%put_var(varid, start, ja, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_vf_rows: '//errmsg)
  end subroutine put_vf_rows

  subroutine get_vf_rows(this, val, ja, start)
    class(rad_encl_file), intent(in) :: this
    real(real32), intent(out) :: val(:)
    integer(int32), intent(out) :: ja(:)
    integer(c_size_t), intent(in) :: start
    integer :: stat, varid
    character(:), allocatable :: errmsg
    ASSERT(size(ja) == size(val))
    if (size(ja) == 0) return
    !! Should check that we're not trying to read too much data
    call this%file%inq_varid('val', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_vf_rows: '//errmsg)
    call this%file%get_var(varid, start, val, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_vf_rows: '//errmsg)
    call this%file%inq_varid('ja', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_vf_rows: '//errmsg)
    call this%file%get_var(varid, start, ja, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_vf_rows: '//errmsg)
  end subroutine get_vf_rows

  subroutine put_ambient(this, ambient)
    class(rad_encl_file), intent(in) :: this
    real(real32), intent(in) :: ambient(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('ambient', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_ambient: '//errmsg)
    call this%file%put_var(varid, ambient, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_ambient: '//errmsg)
  end subroutine put_ambient

  subroutine get_ambient(this, ambient)
    class(rad_encl_file), intent(in) :: this
    real(real32), intent(out) :: ambient(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('ambient', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_ambient: '//errmsg)
    call this%file%get_var(varid, ambient, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_ambient: '//errmsg)
  end subroutine get_ambient

  subroutine put_f2p_map(this, f2p_map)
    class(rad_encl_file), intent(in) :: this
    integer(int32), intent(in) :: f2p_map(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('f2p_map', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_f2p_map: '//errmsg)
    call this%file%put_var(varid, f2p_map, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_f2p_map: '//errmsg)
  end subroutine put_f2p_map

  subroutine get_f2p_map(this, f2p_map)
    class(rad_encl_file), intent(in) :: this
    integer(int32), intent(out) :: f2p_map(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('f2p_map', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_f2p_map: '//errmsg)
    call this%file%get_var(varid, f2p_map, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_f2p_map: '//errmsg)
  end subroutine get_f2p_map

  subroutine put_area(this, area)
    class(rad_encl_file), intent(in) :: this
    real(real64), intent(in) :: area(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('area', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_area: '//errmsg)
    call this%file%put_var(varid, area, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_area: '//errmsg)
  end subroutine put_area

  subroutine get_area(this, area)
    class(rad_encl_file), intent(in) :: this
    real(real64), intent(out) :: area(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('area', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_area: '//errmsg)
    call this%file%get_var(varid, area, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_area: '//errmsg)
  end subroutine get_area

  subroutine put_face_weight(this, face_weight)
    class(rad_encl_file), intent(in) :: this
    real(real64), intent(in) :: face_weight(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('face_weight', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_face_weight: '//errmsg)
    call this%file%put_var(varid, face_weight, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%put_face_weight: '//errmsg)
  end subroutine put_face_weight

  subroutine get_face_weight(this, face_weight)
    class(rad_encl_file), intent(in) :: this
    real(real64), intent(out) :: face_weight(:)
    integer :: stat, varid
    character(:), allocatable :: errmsg
    !! Should check the length is correct
    call this%file%inq_varid('face_weight', varid, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_face_weight: '//errmsg)
    call this%file%get_var(varid, face_weight, stat, errmsg)
    if (stat /= 0) call error_exit('rad_encl_file%get_face_weight: '//errmsg)
  end subroutine get_face_weight

  logical function has_vf_data(this)
    class(rad_encl_file), intent(in) :: this
    integer :: varid, stat
    call this%file%inq_varid('val', varid, stat)
    has_vf_data = (stat == 0)
  end function has_vf_data

  logical function has_ambient(this)
    class(rad_encl_file), intent(in) :: this
    integer :: varid, stat
    call this%file%inq_varid('ambient', varid, stat)
    has_ambient = (stat == 0)
  end function has_ambient

  logical function has_patches(this)
    class(rad_encl_file), intent(in) :: this
    integer :: varid, stat
    call this%file%inq_varid('f2p_map', varid, stat)
    has_patches = (stat == 0)
  end function has_patches

  logical function has_area(this)
    class(rad_encl_file), intent(in) :: this
    integer :: varid, stat
    call this%file%inq_varid('area', varid, stat)
    has_area = (stat == 0)
  end function has_area

  logical function has_face_weight(this)
    class(rad_encl_file), intent(in) :: this
    integer :: varid, stat
    call this%file%inq_varid('face_weight', varid, stat)
    has_face_weight = (stat == 0)
  end function has_face_weight

  subroutine error_exit(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
#if defined(NAGFOR) | defined(NAG_COMPILER)
    use,intrinsic :: f90_unix, only: exit
#endif
    character(*), intent(in) :: errmsg
    write(error_unit,'(2a)') 'error: ', errmsg
    call exit(1)
  end subroutine error_exit

end module rad_encl_file_type
