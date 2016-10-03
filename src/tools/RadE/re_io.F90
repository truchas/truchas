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

module re_io

  use netcdf
  use,intrinsic :: iso_fortran_env, only: i8 => int64
  use,intrinsic :: iso_c_binding, only: c_int, c_size_t, c_char, c_null_char, c_float
  implicit none
  private
  
  public :: re_open_rw, re_close, re_init_vf, re_put_ia, re_put_vf_rows, re_put_ambient
  public :: re_open_ro, re_get_vf_dims, re_get_ia, re_get_vf_rows, re_get_ambient
  public :: re_put_icount, re_get_icount

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
    function nc_inq_dimid(ncid, name, dimidp) result(status) bind(c,name='nc_inq_dimid')
      import c_int, c_char
      integer(c_int), value :: ncid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: dimidp
      integer(c_int) :: status
    end function
    function nc_inq_dimlen(ncid, dimid, lenp) result(status) bind(c,name='nc_inq_dimlen')
      import c_int, c_size_t
      integer(c_int), value :: ncid, dimid
      integer(c_size_t) :: lenp
      integer(c_int) :: status
    end function
    function nc_def_var(ncid, name, xtype, ndims, dimids, varidp) result(status) &
        bind(c,name='nc_def_var')
      import c_int, c_char
      integer(c_int), value :: ncid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), value :: xtype, ndims
      integer(c_int), intent(in) :: dimids(*)
      integer(c_int), intent(out) :: varidp
      integer(c_int) :: status
    end function
    function nc_inq_varid(ncid, name, varidp) result(status) bind(c,name='nc_inq_varid')
      import c_int, c_char
      integer(c_int), value :: ncid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: varidp
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

  !! Values from netcdf.h (we count on them not changing)
  integer, parameter :: NC_INT = 4
  integer, parameter :: NC_FLOAT = 5

contains
  
  subroutine re_open_rw (path, ncid)
    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid
    integer :: status
    status = nf90_open(path, NF90_WRITE, ncid)
    call handle_netcdf_error (status)
  end subroutine re_open_rw
  
  subroutine re_open_ro (path, ncid)
    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid
    integer :: status
    status = nf90_open(path, NF90_NOWRITE, ncid)
    call handle_netcdf_error (status)
  end subroutine re_open_ro
  
  subroutine re_close (ncid)
    integer, intent(in) :: ncid
    integer :: status
    status = nf90_close (ncid)
    call handle_netcdf_error (status)
  end subroutine re_close
  
  subroutine re_init_vf (ncid, nnonz)
  
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
    status = nc_def_dim(ncid, 'num_nonzero'//c_null_char, nnonz, dim2)
    call handle_netcdf_error (status)
    status = nc_def_var(ncid, 'val'//c_null_char, NC_FLOAT, 1, [dim2], varid)
    call handle_netcdf_error (status)
    status = nc_def_var(ncid, 'ja'//c_null_char, NC_INT, 1, [dim2], varid)
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
    
  end subroutine re_init_vf
  
  subroutine re_get_vf_dims (ncid, nface, nnonz)
  
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
    status = nc_inq_dimid(ncid, 'num_nonzero'//c_null_char, dimid)
    call handle_netcdf_error (status)
    status = nc_inq_dimlen(ncid, dimid, nnonz)
    call handle_netcdf_error (status)

  end subroutine re_get_vf_dims

  !! Keep for backward compatibility with client code.
  subroutine re_put_ia (ncid, ia)
    integer, intent(in) :: ncid, ia(:)
    call re_put_icount(ncid, ia(2:)-ia(:size(ia)-1))
  end subroutine re_put_ia

  subroutine re_put_icount (ncid, icount)
    integer, intent(in) :: ncid, icount(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'icount', varid)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, varid, icount)
    call handle_netcdf_error (status, 'writing matrix row nonzero counts')
  end subroutine re_put_icount

  !! Keep for backward compatibility with client code.
  subroutine re_get_ia (ncid, ia)
    integer, intent(in)  :: ncid
    integer, intent(out) :: ia(:)
    integer :: j
    call re_get_icount (ncid, ia(2:))
    ia(1) = 1
    do j = 2, size(ia)
      ia(j) = ia(j) + ia(j-1)
    end do
  end subroutine re_get_ia

  subroutine re_get_icount (ncid, icount)
    integer, intent(in)  :: ncid
    integer, intent(out) :: icount(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'icount', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, icount)
    call handle_netcdf_error (status, 'reading matrix row nonzero counts')
  end subroutine re_get_icount

  subroutine re_put_vf_rows (ncid, val, ja, start)
    integer, intent(in) :: ncid, ja(:)
    integer(i8), intent(in) :: start
    real, intent(in) :: val(:)
    integer :: status, varid
    ASSERT( size(ja) == size(val) )
    if (size(ja) == 0) return
    !! Should check that we're not trying to write too much data
    status = nc_inq_varid(ncid, 'val'//c_null_char, varid)
    call handle_netcdf_error (status)
    status = nc_put_vara_float(ncid, varid, [start-1], [size(val,kind=i8)], val)
    call handle_netcdf_error (status, 'writing matrix row values')
    status = nc_inq_varid(ncid, 'ja'//c_null_char, varid)
    call handle_netcdf_error (status)
    status = nc_put_vara_int(ncid, varid, [start-1], [size(ja,kind=i8)], ja)
    call handle_netcdf_error (status, 'writing matrix row indices')
  end subroutine
  
  subroutine re_get_vf_rows (ncid, val, ja, start)
    integer, intent(in)  :: ncid
    integer(i8), intent(in) :: start
    integer, intent(out) :: ja(:)
    real,    intent(out) :: val(:)
    integer :: status, varid
    ASSERT( size(ja) == size(val) )
    if (size(ja) == 0) return
    !! Should check that we're not trying to read too much data
    status = nc_inq_varid(ncid, 'val'//c_null_char, varid)
    call handle_netcdf_error (status)
    status = nc_get_vara_float(ncid, varid, [start-1], [size(val,kind=i8)], val)
    call handle_netcdf_error (status, 'reading matrix row values')
    status = nc_inq_varid(ncid, 'ja'//c_null_char, varid)
    call handle_netcdf_error (status)
    status = nc_get_vara_int(ncid, varid, [start-1], [size(ja,kind=i8)], ja)
    call handle_netcdf_error (status, 'reading matrix row indices')
  end subroutine
  
  subroutine re_put_ambient (ncid, ambient)
    integer, intent(in) :: ncid
    real,    intent(in) :: ambient(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'ambient', varid)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, varid, ambient)
    call handle_netcdf_error (status, 'writing ambient view factors')
  end subroutine re_put_ambient
  
  subroutine re_get_ambient (ncid, ambient)
    integer, intent(in)  :: ncid
    real,    intent(out) :: ambient(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'ambient', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, ambient)
    call handle_netcdf_error (status, 'reading ambient view factors')
  end subroutine re_get_ambient
  
  subroutine handle_netcdf_error (status, errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    integer, intent(in) :: status
    character(*), intent(in), optional :: errmsg
    if (status == NF90_NOERR) return
    if (present(errmsg)) then
      write(error_unit,'(4a)') 'ERROR: ', errmsg, ': ', trim(nf90_strerror(status))
    else
      write(error_unit,'(2a)') 'ERROR: ', trim(nf90_strerror(status))
    end if
    stop
  end subroutine handle_netcdf_error
  
end module re_io
