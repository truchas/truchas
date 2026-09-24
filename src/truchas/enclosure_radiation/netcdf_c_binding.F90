!!
!! NETCDF_C_BINDING
!!
!! Raw bindings to a subset of the netCDF C library interface.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! October 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!!  These are just the raw bindings. See NETCDF_FILE_TYPE for a Fortran-centric
!!  interface that handles null termination of string arguments, reversing the
!!  order of array indices, etc.
!!

module netcdf_c_binding

  use,intrinsic :: iso_c_binding, only: c_int, c_long_long, c_ptr, c_char, c_size_t
  use,intrinsic :: iso_c_binding, only: c_float, c_double
  implicit none
  private

  public :: nc_create, nc_open, nc_close, nc_redef, nc_enddef, nc_strerror
  public :: nc_def_dim, nc_inq_dimid, nc_inq_dimlen, nc_def_var, nc_inq_varid
  public :: nc_put_var_int, nc_put_var_longlong, nc_put_var_float, nc_put_var_double
  public :: nc_get_var_int, nc_get_var_longlong, nc_get_var_float, nc_get_var_double
  public :: nc_put_vara_int, nc_put_vara_longlong, nc_put_vara_float, nc_put_vara_double
  public :: nc_get_vara_int, nc_get_vara_longlong, nc_get_vara_float, nc_get_vara_double

  !! File creation modes from netcdf.h
  integer(c_int), parameter, public :: NC_CLOBBER   = int(z'0000')
  integer(c_int), parameter, public :: NC_NOCLOBBER = int(z'0004')
  integer(c_int), parameter, public :: NC_NETCDF4   = int(z'1000')
  integer(c_int), parameter, public :: NC_64BIT_OFFSET = int(z'0200')

  !! File open modes from netcdf.h
  integer(c_int), parameter, public :: NC_NOWRITE = 0, NC_WRITE = 1

  !! NetCDF-4 types from netcdf.h
  integer(c_int), parameter, public :: NC_INT    = 4  ! signed 4-byte integer
  integer(c_int), parameter, public :: NC_FLOAT  = 5  ! single precision fp number
  integer(c_int), parameter, public :: NC_DOUBLE = 6  ! double precision fp number
  integer(c_int), parameter, public :: NC_INT64  = 10 ! signed 8-byte integer

  !! Select error codes from netcdf.h
  integer(c_int), parameter, public :: NC_NOERR = 0
  integer(c_int), parameter, public :: NC_EBADDIM = -46

  interface
    function nc_create(path, cmode, ncidp) result(ncerr) bind(c,name='nc_create')
      import c_int, c_char
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: cmode
      integer(c_int) :: ncidp, ncerr
    end function
    function nc_open(path, omode, ncidp) result(ncerr) bind(c,name='nc_open')
      import c_int, c_char
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: omode
      integer(c_int) :: ncidp, ncerr
    end function
    function nc_close(ncid) result(ncerr) bind(c,name='nc_close')
      import c_int
      integer(c_int), value :: ncid
      integer(c_int) :: ncerr
    end function
    function nc_redef(ncid) result(ncerr) bind(c,name='nc_redef')
      import c_int
      integer(c_int), value :: ncid
      integer(c_int) :: ncerr
    end function
    function nc_enddef(ncid) result(ncerr) bind(c,name='nc_enddef')
      import c_int
      integer(c_int), value :: ncid
      integer(c_int) :: ncerr
    end function
    function nc_def_dim(ncid, name, len, idp) result(ncerr) bind(c,name='nc_def_dim')
      import c_char, c_int, c_size_t
      integer(c_int), value :: ncid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_size_t), value :: len
      integer(c_int) :: idp
      integer(c_int) :: ncerr
    end function
    function nc_inq_dimid(ncid, name, dimidp) result(ncerr) bind(c,name='nc_inq_dimid')
      import c_int, c_char
      integer(c_int), value :: ncid
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: dimidp
      integer(c_int) :: ncerr
    end function
    function nc_inq_dimlen(ncid, dimid, lenp) result(ncerr) bind(c,name='nc_inq_dimlen')
      import c_int, c_size_t
      integer(c_int), value :: ncid, dimid
      integer(c_size_t) :: lenp
      integer(c_int) :: ncerr
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
    function nc_put_var_int(ncid, varid, op) result(status) &
        bind(c,name='nc_put_var_int')
      import c_int
      integer(c_int), value :: ncid, varid
      integer(c_int), intent(in) :: op(*)
      integer(c_int) :: status
    end function
    function nc_put_var_longlong(ncid, varid, op) result(status) &
        bind(c,name='nc_put_var_longlong')
      import c_int, c_long_long
      integer(c_int), value :: ncid, varid
      integer(c_long_long), intent(in) :: op(*)
      integer(c_int) :: status
    end function
    function nc_put_var_float(ncid, varid, op) result(status) &
        bind(c,name='nc_put_var_float')
      import c_int, c_float
      integer(c_int), value :: ncid, varid
      real(c_float), intent(in) :: op(*)
      integer(c_int) :: status
    end function
    function nc_put_var_double(ncid, varid, op) result(status) &
        bind(c,name='nc_put_var_double')
      import c_int, c_double
      integer(c_int), value :: ncid, varid
      real(c_double), intent(in) :: op(*)
      integer(c_int) :: status
    end function
    function nc_get_var_int(ncid, varid, ip) result(status) &
        bind(c,name='nc_get_var_int')
      import c_int
      integer(c_int), value :: ncid, varid
      integer(c_int), intent(out) :: ip(*)
      integer(c_int) :: status
    end function
    function nc_get_var_longlong(ncid, varid, ip) result(status) &
        bind(c,name='nc_get_var_longlong')
      import c_int, c_long_long
      integer(c_int), value :: ncid, varid
      integer(c_long_long), intent(out) :: ip(*)
      integer(c_int) :: status
    end function
    function nc_get_var_float(ncid, varid, ip) result(status) &
        bind(c,name='nc_get_var_float')
      import c_int, c_float
      integer(c_int), value :: ncid, varid
      real(c_float), intent(out) :: ip(*)
      integer(c_int) :: status
    end function
    function nc_get_var_double(ncid, varid, ip) result(status) &
        bind(c,name='nc_get_var_double')
      import c_int, c_double
      integer(c_int), value :: ncid, varid
      real(c_double), intent(out) :: ip(*)
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
    function nc_put_vara_longlong(ncid, varid, startp, countp, op) result(status) &
        bind(c,name='nc_put_vara_longlong')
      import c_int, c_size_t, c_long_long
      integer(c_int), value :: ncid, varid
      integer(c_size_t), intent(in) :: startp(*), countp(*)
      integer(c_long_long), intent(in) :: op(*)
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
    function nc_put_vara_double(ncid, varid, startp, countp, op) result(status) &
        bind(c,name='nc_put_vara_double')
      import c_int, c_size_t, c_double
      integer(c_int), value :: ncid, varid
      integer(c_size_t), intent(in) :: startp(*), countp(*)
      real(c_double), intent(in) :: op(*)
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
    function nc_get_vara_longlong(ncid, varid, startp, countp, ip) result(status) &
        bind(c,name='nc_get_vara_longlong')
      import c_int, c_size_t, c_long_long
      integer(c_int), value :: ncid, varid
      integer(c_size_t), intent(in) :: startp(*), countp(*)
      integer(c_long_long), intent(out) :: ip(*)
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
    function nc_get_vara_double(ncid, varid, startp, countp, ip) result(status) &
        bind(c,name='nc_get_vara_double')
      import c_int, c_size_t, c_double
      integer(c_int), value :: ncid, varid
      integer(c_size_t), intent(in) :: startp(*), countp(*)
      real(c_double), intent(out) :: ip(*)
      integer(c_int) :: status
    end function
    function nc_strerror(ncerr) result(strp) bind(c,name='nc_strerror')
      import c_int, c_ptr
      integer(c_int), value :: ncerr
      type(c_ptr) :: strp
    end function
  end interface

end module netcdf_c_binding
