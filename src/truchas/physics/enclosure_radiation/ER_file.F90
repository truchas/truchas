!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! IMPLEMENTATION NOTES
!!
!! With one exception, the present implementation limits the problem size to
!! that representable with default integer variables.  That means all positive
!! integer values are limited to 2**31 (about 2e9).  The exception is the size
!! of the view factor matrix which is represented with a 64 bit integer.
!!
!! 1. While the native C netCDF interface uses size_t values (typically 64 bit)
!! to represent dimension sizes, its Fortran interface uses default integer,
!! which are 32 bits.  Thus to define/retrieve a 64-bit dimension we have to
!! use the C functions.
!!
!! 2. Because the (global) view factor matrix size is limited by the largest
!! 64-bit integer, the IA array of the standard CSR representation would also
!! need to use 64-bit integers.  An alternative for the disk file, is to
!! instead store the number of non-zeros per row.  These are legitimately
!! limited to 32-bit integers.
!!
!! 3. The C interface (nc_*) apparently uses 0-based IDs for variables and
!! dimensions, whereas the Fortran interface (nf90_*) uses 1-based IDs.
!! Consequently an ID returned by an nc_* call must be passed as ID+1 to a
!! nf90_* call, and an ID returned by an nf90_* call must be passed as ID-1
!! to an nc_* call. (I frankly do not understand why there is this difference
!! in the first place; the IDs mean absolutely nothing to the user code.)
!!

#include "f90_assert.fpp"

module ER_file

  use,intrinsic :: iso_fortran_env, only: r8 => real64, i8 => int64
  use,intrinsic :: iso_c_binding, only: c_int, c_size_t, c_char, c_null_char, c_float
  use netcdf
  implicit none
  private

  public :: ERF_create, ERF_open_ro, ERF_open_rw, ERF_close, ERF_get_dims
  public :: ERF_put_surface, ERF_get_surface
  public :: ERF_put_symmetry_info, ERF_get_symmetry_info
  public :: ERF_put_group_info, ERF_get_group_info
  public :: ERF_put_source_info, ERF_get_source_info
  public :: ERF_init_vf, ERF_get_vf_dims
  public :: ERF_put_ia, ERF_get_ia, ERF_get_vf_rowcount
  public :: ERF_put_vf_rows, ERF_get_vf_rows
  public :: ERF_put_ambient, ERF_get_ambient
  public :: ERF_put_icount, ERF_get_icount

  !! Bindings to netCDF C interface functions. Needed for large dimensions; see Note 1.
  interface
    function nc_def_dim(ncid, name, len, idp) result(status) bind(c,name='nc_def_dim')
      import c_char, c_int, c_size_t
      integer(c_int), value :: ncid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_size_t), value :: len
      integer(c_int) :: idp
      integer(c_int) :: status
    end function
    function nc_inq_dimlen(ncid, dimid, lenp) result(status) bind(c,name='nc_inq_dimlen')
      import c_int, c_size_t
      integer(c_int), value :: ncid, dimid
      integer(c_size_t) :: lenp
      integer(c_int) :: status
    end function
    function nc_put_vara_int(ncid, varid, startp, countp, op) result(status) &
        bind(c,name='nc_put_vara_int')
      import c_int, c_size_t
      integer(c_int), value :: ncid, varid
      integer(c_size_t), intent(in) :: startp(*), countp(*)
      integer(c_int), intent(in) :: op(*)
      integer(c_int) :: status
    end function
    function nc_put_vara_float(ncid, varid, startp, countp, op) result(status) &
        bind(c,name='nc_put_vara_float')
      import c_int, c_size_t, c_float
      integer(c_int), value :: ncid, varid
      integer(c_size_t), intent(in) :: startp(*), countp(*)
      real(c_float), intent(in) :: op(*)
      integer(c_int) :: status
    end function
    function nc_get_vara_int(ncid, varid, startp, countp, ip) result(status) &
        bind(c,name='nc_get_vara_int')
      import c_int, c_size_t
      integer(c_int), value :: ncid, varid
      integer(c_size_t), intent(in) :: startp(*), countp(*)
      integer(c_int), intent(out) :: ip(*)
      integer(c_int) :: status
    end function
    function nc_get_vara_float(ncid, varid, startp, countp, ip) result(status) &
        bind(c,name='nc_get_vara_float')
      import c_int, c_size_t, c_float
      integer(c_int), value :: ncid, varid
      integer(c_size_t), intent(in) :: startp(*), countp(*)
      real(c_float), intent(out) :: ip(*)
      integer(c_int) :: status
    end function
  end interface

contains

  subroutine ERF_create (path, ncid)

    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid

    integer :: cmode, status

    cmode = ior(NF90_CLOBBER, NF90_HDF5)
    status = nf90_create(path, cmode, ncid)
    call handle_netcdf_error (status)

    !! Take the dataset out of define mode.
    status = nf90_enddef (ncid)
    call handle_netcdf_error (status)

  end subroutine ERF_create

  subroutine ERF_open_rw (path, ncid)
    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid
    integer :: status
    status = nf90_open(trim(path), NF90_WRITE, ncid)
    call handle_netcdf_error (status)
  end subroutine ERF_open_rw

  subroutine ERF_open_ro (path, ncid)
    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid
    integer :: status
    status = nf90_open(path, NF90_NOWRITE, ncid)
    call handle_netcdf_error (status)
  end subroutine ERF_open_ro

  subroutine ERF_close (ncid)
    integer, intent(in) :: ncid
    integer :: status
    status = nf90_close (ncid)
    call handle_netcdf_error (status)
  end subroutine ERF_close

  subroutine ERF_put_surface (ncid, xface, fnode, coord)

    integer,  intent(in) :: ncid
    integer,  intent(in) :: xface(:), fnode(:)
    real(r8), intent(in) :: coord(:,:)

    integer :: status, dim1, dim2
    integer :: coord_id, xface_id, fnode_id

    !! Put the dataset into define mode.
    status = nf90_redef (ncid)
    call handle_netcdf_error (status)

    !! Define the coordinate variable.
    status = nf90_def_dim(ncid, 'num_dim', size(coord,dim=1), dim1)
    call handle_netcdf_error (status)
    status = nf90_def_dim(ncid, 'num_nodes', size(coord,dim=2), dim2)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'coord', NF90_DOUBLE, (/dim1,dim2/), coord_id)
    call handle_netcdf_error (status)

    !! Define the face connectivity variables.
    status = nf90_def_dim(ncid, 'num_faces', size(xface)-1, dim1)
    call handle_netcdf_error (status)
    status = nf90_def_dim(ncid, 'num_xface', size(xface), dim1)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'xface', NF90_INT, (/dim1/), xface_id)
    call handle_netcdf_error (status)
    status = nf90_def_dim(ncid, 'num_fnode', size(fnode), dim1)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'fnode', NF90_INT, (/dim1/), fnode_id)
    call handle_netcdf_error (status)

    !! Take the dataset out of define mode.
    status = nf90_enddef (ncid)
    call handle_netcdf_error (status)

    !! Write the coordinate data.
    status = nf90_put_var(ncid, coord_id, coord)
    call handle_netcdf_error (status)

    !! Write the face connectivity data.
    status = nf90_put_var(ncid, xface_id, xface)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, fnode_id, fnode)
    call handle_netcdf_error (status)

  end subroutine ERF_put_surface

  subroutine ERF_get_dims (ncid, nface, nnode, nfnode, ngroup)

    integer, intent(in)  :: ncid
    integer, intent(out) :: nface, nnode, nfnode, ngroup

    integer :: status, dimid

    !! Get the number of faces.
    status = nf90_inq_dimid(ncid, 'num_faces', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=nface)
    call handle_netcdf_error (status)

    !! Get the number of nodes.
    status = nf90_inq_dimid(ncid, 'num_nodes', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=nnode)
    call handle_netcdf_error (status)

    !! Get the size of the fnode array.
    status = nf90_inq_dimid(ncid, 'num_fnode', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=nfnode)
    call handle_netcdf_error (status)

    !! Get the number of face groups.
    status = nf90_inq_dimid(ncid, 'num_group', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=ngroup)
    call handle_netcdf_error (status)

  end subroutine ERF_get_dims

  subroutine ERF_get_surface (ncid, fsize, fnode, coord)

    integer,  intent(in)  :: ncid
    integer,  intent(out) :: fsize(:), fnode(:)
    real(r8), intent(out) :: coord(:,:)

    integer :: status, dimid, n, varid
    integer, allocatable :: xface(:)

    status = nf90_inq_dimid(ncid, 'num_dim', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    INSIST(size(coord,dim=1) == n)

    status = nf90_inq_dimid(ncid, 'num_nodes', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    INSIST(size(coord,dim=2) == n)

    status = nf90_inq_dimid(ncid, 'num_faces', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    INSIST(size(fsize) == n)

    status = nf90_inq_dimid(ncid, 'num_fnode', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    INSIST(size(fnode) == n)

    !! Read the fsize data.
    allocate(xface(size(fsize)+1))
    status = nf90_inq_varid(ncid, 'xface', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, xface)
    call handle_netcdf_error (status)
    fsize = xface(2:) - xface(:size(fsize))
    deallocate(xface)

    !! Read the fnode data.
    status = nf90_inq_varid(ncid, 'fnode', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, fnode)
    call handle_netcdf_error (status)

    !! Read the coord data.
    status = nf90_inq_varid(ncid, 'coord', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, coord)
    call handle_netcdf_error (status)

  end subroutine ERF_get_surface

  subroutine ERF_put_group_info (ncid, gnum, group_ids)

    integer, intent(in) :: ncid
    integer, intent(in) :: gnum(:), group_ids(:)

    integer :: status, dimid, num_faces, gnum_id, group_id

    !! Get the num_faces dimension from the database.
    status = nf90_inq_dimid(ncid, 'num_faces', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=num_faces)
    call handle_netcdf_error (status)

    INSIST(size(gnum) == num_faces)
    INSIST(minval(gnum) > 0)
    INSIST(maxval(gnum) <= size(group_ids))

    !! Put the dataset into define mode.
    status = nf90_redef (ncid)
    call handle_netcdf_error (status)

    !! Define the gnum variable.
    status = nf90_def_var (ncid, 'gnum', NF90_INT, (/dimid/), gnum_id)
    call handle_netcdf_error (status)

    !! Define the group_ids variable.
    status = nf90_def_dim(ncid, 'num_group', size(group_ids), dimid)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'group_ids', NF90_INT, (/dimid/), group_id)
    call handle_netcdf_error (status)

    !! Take the dataset out of define mode.
    status = nf90_enddef (ncid)
    call handle_netcdf_error (status)

    !! Write the gnum data.
    status = nf90_put_var (ncid, gnum_id, gnum)
    call handle_netcdf_error (status)

    !! Write the group_ids data.
    status = nf90_put_var (ncid, group_id, group_ids)
    call handle_netcdf_error (status)

  end subroutine ERF_put_group_info

  subroutine ERF_get_group_info (ncid, gnum, group_ids)

    integer, intent(in) :: ncid
    integer, intent(out) :: gnum(:), group_ids(:)

    integer :: status, varid

    status = nf90_inq_varid(ncid, 'gnum', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, gnum)
    call handle_netcdf_error (status)

    status = nf90_inq_varid(ncid, 'group_ids', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, group_ids)
    call handle_netcdf_error (status)

  end subroutine ERF_get_group_info

  subroutine ERF_put_symmetry_info (ncid, mirror, rot_axis, num_rot)

    integer, intent(in) :: ncid
    logical, intent(in) :: mirror(:)
    integer, intent(in) :: rot_axis
    integer, intent(in) :: num_rot

    integer :: status, dim1, varid
    integer :: symmetry(5)

    INSIST(size(mirror) == 3)
    INSIST(rot_axis >= 0 .and. rot_axis <=3)

    !! Put the dataset into define mode.
    status = nf90_redef (ncid)
    call handle_netcdf_error (status)

    !! Define the symmetry variable.
    status = nf90_def_dim(ncid, 'num_symmetry', size(symmetry), dim1)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'symmetry', NF90_INT, (/dim1/), varid)
    call handle_netcdf_error (status)

    !! Take the dataset out of define mode.
    status = nf90_enddef (ncid)
    call handle_netcdf_error (status)

    !! Write the serialized symmetry info.
    symmetry = 0
    where (mirror) symmetry(1:3) = 1
    symmetry(4) = rot_axis
    if (rot_axis > 0) then
      INSIST(num_rot > 1)
      symmetry(5) = num_rot
    end if
    status = nf90_put_var (ncid, varid, symmetry)
    call handle_netcdf_error (status)

  end subroutine ERF_put_symmetry_info

  subroutine ERF_get_symmetry_info (ncid, mirror, rot_axis, num_rot)

    integer, intent(in)  :: ncid
    logical, intent(out) :: mirror(:)
    integer, intent(out) :: rot_axis
    integer, intent(out) :: num_rot

    integer :: status, varid
    integer :: symmetry(5)

    INSIST(size(mirror) == 3)

    !! Read the serialized symmetry info.
    status = nf90_inq_varid(ncid, 'symmetry', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, symmetry)
    call handle_netcdf_error (status)

    !! Unpack the serialized symmetry info.
    mirror  = (symmetry(1:3) /= 0)
    rot_axis = symmetry(4)
    num_rot  = symmetry(5)

  end subroutine ERF_get_symmetry_info

  subroutine ERF_put_source_info (ncid, src_elem, src_side)

    integer, intent(in) :: ncid
    integer, intent(in) :: src_elem(:), src_side(:)

    integer :: status, dimid, num_faces, src_elem_id, src_side_id

    !! Get num_faces from the database.
    status = nf90_inq_dimid(ncid, 'num_faces', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=num_faces)
    call handle_netcdf_error (status)

    INSIST(size(src_elem) == num_faces)
    INSIST(size(src_side) == num_faces)

    !! Put the dataset into define mode.
    status = nf90_redef (ncid)
    call handle_netcdf_error (status)

    !! Define the src_elem and src_side variables.
    status = nf90_def_var (ncid, 'src_elem', NF90_INT, (/dimid/), src_elem_id)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'src_side', NF90_INT, (/dimid/), src_side_id)
    call handle_netcdf_error (status)

    !! Take the dataset out of define mode.
    status = nf90_enddef (ncid)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, src_elem_id, src_elem)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, src_side_id, src_side)
    call handle_netcdf_error (status)

  end subroutine ERF_put_source_info

  subroutine ERF_get_source_info (ncid, src_elem, src_side)

    integer, intent(in)  :: ncid
    integer, intent(out) :: src_elem(:), src_side(:)

    integer :: status, varid

    status = nf90_inq_varid(ncid, 'src_elem', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, src_elem)
    call handle_netcdf_error (status)

    status = nf90_inq_varid(ncid, 'src_side', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, src_side)
    call handle_netcdf_error (status)

  end subroutine ERF_get_source_info

  subroutine ERF_init_vf (ncid, nnonz)

    integer, intent(in) :: ncid
    integer(c_size_t), intent(in) :: nnonz

    integer :: status, dim1, dim2, varid, nface

    status = nf90_redef(ncid)
    call handle_netcdf_error (status)

    !! Get the number of faces and its dimension ID.
    status = nf90_inq_dimid(ncid, 'num_faces', dim1)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dim1, len=nface)
    call handle_netcdf_error (status)

    !! Define the VAL and JA arrays.
    !status = nf90_def_dim(ncid, 'num_nonzero', nnonz, dim2)
    status = nc_def_dim(ncid, 'num_nonzero'//c_null_char, nnonz, dim2)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'val', NF90_FLOAT, (/dim2+1/), varid)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'ja', NF90_INT, (/dim2+1/), varid)
    call handle_netcdf_error (status)

    ! Replace IA with ICOUNT; see Note 2.
    ! !! Define the IA array.
    ! status = nf90_def_dim(ncid, 'num_ia', nface+1, dim2)
    ! call handle_netcdf_error (status)
    ! status = nf90_def_var(ncid, 'ia', NF90_INT, (/dim2/), varid)
    ! call handle_netcdf_error (status)

    !! Define the ICOUNT array.
    status = nf90_def_var(ncid, 'icount', NF90_INT, (/dim1/), varid)
    call handle_netcdf_error (status)

    !! Define the ambient view factor array.
    status = nf90_def_var(ncid, 'ambient', NF90_FLOAT, (/dim1/), varid)
    call handle_netcdf_error (status)

    status = nf90_enddef(ncid)
    call handle_netcdf_error (status)

  end subroutine ERF_init_vf

  subroutine ERF_get_vf_dims (ncid, nface, nnonz)

    integer, intent(in)  :: ncid
    integer, intent(out) :: nface
    integer(c_size_t), intent(out) :: nnonz

    integer :: status, dimid

    !! Get the number of faces.
    status = nf90_inq_dimid(ncid, 'num_faces', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=nface)
    call handle_netcdf_error (status)

    !! Get the number of nonzeros in the VF matrix.
    status = nf90_inq_dimid(ncid, 'num_nonzero', dimid)
    call handle_netcdf_error (status)
    !status = nf90_inquire_dimension(ncid, dimid, len=nnonz)
    status = nc_inq_dimlen(ncid, dimid-1, nnonz)  ! see Note 3
    call handle_netcdf_error (status)

  end subroutine ERF_get_vf_dims

  !! Keep for backward compatibility with client code.
  subroutine ERF_put_ia (ncid, ia)
    integer, intent(in) :: ncid, ia(:)
    call ERF_put_icount(ncid, ia(2:)-ia(:size(ia)-1))
  end subroutine ERF_put_ia

  subroutine ERF_put_icount (ncid, icount)
    integer, intent(in) :: ncid, icount(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'icount', varid)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, varid, icount)
    call handle_netcdf_error (status)
  end subroutine ERF_put_icount

  !! Keep for backward compatibility with client code.
  subroutine ERF_get_ia (ncid, ia)
    integer, intent(in)  :: ncid
    integer, intent(out) :: ia(:)
    integer :: j
    call ERF_get_icount (ncid, ia(2:))
    ia(1) = 1
    do j = 2, size(ia)
      ia(j) = ia(j) + ia(j-1)
    end do
  end subroutine ERF_get_ia

  subroutine ERF_get_icount (ncid, icount)
    integer, intent(in)  :: ncid
    integer, intent(out) :: icount(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'icount', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, icount)
    call handle_netcdf_error (status)
  end subroutine ERF_get_icount

  subroutine ERF_get_vf_rowcount (ncid, rowcount)
    integer, intent(in)  :: ncid
    integer, intent(out) :: rowcount(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'icount', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, rowcount)
    call handle_netcdf_error (status)
  end subroutine ERF_get_vf_rowcount

  subroutine ERF_put_vf_rows (ncid, val, ja, start)
    integer, intent(in) :: ncid, ja(:)
    integer(i8), intent(in) :: start
    real, intent(in) :: val(:)
    integer :: status, varid
    ASSERT( size(ja) == size(val) )
    if (size(ja) == 0) return
    !! Should check that we're not trying to write too much data
    status = nf90_inq_varid(ncid, 'val', varid)
    call handle_netcdf_error (status)
    !status = nf90_put_var(ncid, varid, val, (/start/))
    status = nc_put_vara_float(ncid, varid-1, [start-1], [size(val,kind=i8)], val) ! see Note 3
    call handle_netcdf_error (status)
    status = nf90_inq_varid(ncid, 'ja', varid)
    call handle_netcdf_error (status)
    !status = nf90_put_var(ncid, varid, ja, (/start/))
    status = nc_put_vara_int(ncid, varid-1, [start-1], [size(ja,kind=i8)], ja) ! see Note 3
    call handle_netcdf_error (status)
  end subroutine

  subroutine ERF_get_vf_rows (ncid, val, ja, start)
    integer, intent(in)  :: ncid
    integer(i8), intent(in) :: start
    integer, intent(out) :: ja(:)
    real,    intent(out) :: val(:)
    integer :: status, varid
    ASSERT( size(ja) == size(val) )
    if (size(ja) == 0) return
    !! Should check that we're not trying to read too much data
    status = nf90_inq_varid(ncid, 'val', varid)
    call handle_netcdf_error (status)
    !status = nf90_get_var(ncid, varid, val, (/start/))
    status = nc_get_vara_float(ncid, varid-1, [start-1], [size(val,kind=i8)], val) ! see Note 3
    call handle_netcdf_error (status)
    status = nf90_inq_varid(ncid, 'ja', varid)
    call handle_netcdf_error (status)
    !status = nf90_get_var(ncid, varid, ja, (/start/))
    status = nc_get_vara_int(ncid, varid-1, [start-1], [size(ja,kind=i8)], ja) ! see Note 3
    call handle_netcdf_error (status)
  end subroutine

  subroutine ERF_put_ambient (ncid, ambient)
    integer, intent(in) :: ncid
    real,    intent(in) :: ambient(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'ambient', varid)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, varid, ambient)
    call handle_netcdf_error (status)
  end subroutine ERF_put_ambient

  subroutine ERF_get_ambient (ncid, ambient)
    integer, intent(in)  :: ncid
    real,    intent(out) :: ambient(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'ambient', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, ambient)
    call handle_netcdf_error (status)
  end subroutine ERF_get_ambient

  subroutine handle_netcdf_error (status)
    integer, intent(in) :: status
    if (status == NF90_NOERR) return
    print *, trim(nf90_strerror(status))
    stop
  end subroutine handle_netcdf_error

end module ER_file
