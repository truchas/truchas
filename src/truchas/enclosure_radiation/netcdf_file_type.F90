!!
!! NETCDF_FILE_TYPE
!!
!! This module defines an object-oriented interface for reading and writing
!! data to and from a netCDF file.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! October 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module netcdf_file_type

  use netcdf_c_binding
  use,intrinsic :: iso_c_binding, only: c_int, c_size_t, c_char, c_null_char, c_ptr, c_f_pointer
  use,intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  integer, parameter, public :: NF_INT32  = NC_INT
  integer, parameter, public :: NF_INT64  = NC_INT64
  integer, parameter, public :: NF_REAL32 = NC_FLOAT
  integer, parameter, public :: NF_REAL64 = NC_DOUBLE

  integer, parameter, public :: NF_EBADDIM = NC_EBADDIM

  !! File creation mode flags
  integer(c_int), parameter, public :: NF_CLOBBER   = NC_CLOBBER
  integer(c_int), parameter, public :: NF_NOCLOBBER = NC_NOCLOBBER
  integer(c_int), parameter, public :: NF_NETCDF4   = NC_NETCDF4

  !! File open mode flags
  integer(c_int), parameter, public :: NF_NOWRITE = NC_NOWRITE
  integer(c_int), parameter, public :: NF_WRITE   = NC_WRITE

  type, public :: netcdf_file
    private
    integer(c_int) :: ncid = -1
  contains
    procedure :: create
    procedure :: open
    procedure :: close
    procedure :: redef
    procedure :: enddef
    procedure :: def_dim
    procedure :: inq_dimid
    procedure :: inq_dimlen
    generic   :: def_var => def_var_1, def_var_n
    procedure :: inq_varid
    generic   :: put_var => put_var_int32, put_var_int64, put_var_real32, put_var_real64, &
                            put_vara_int32, put_vara_int64, put_vara_real32, put_vara_real64, &
                            put_var_real64_2
    generic   :: get_var => get_var_int32, get_var_int64, get_var_real32, get_var_real64, &
                            get_vara_int32, get_vara_int64, get_vara_real32, get_vara_real64, &
                            get_var_real64_2
    procedure, private :: def_var_1
    procedure, private :: def_var_n
    procedure, private :: put_var_int32
    procedure, private :: put_var_int64
    procedure, private :: put_var_real32
    procedure, private :: put_var_real64
    procedure, private :: put_var_real64_2
    procedure, private :: put_vara_int32
    procedure, private :: put_vara_int64
    procedure, private :: put_vara_real32
    procedure, private :: put_vara_real64
    procedure, private :: get_var_int32
    procedure, private :: get_var_int64
    procedure, private :: get_var_real32
    procedure, private :: get_var_real64
    procedure, private :: get_var_real64_2
    procedure, private :: get_vara_int32
    procedure, private :: get_vara_int64
    procedure, private :: get_vara_real32
    procedure, private :: get_vara_real64
  end type netcdf_file

contains

  subroutine create(this, path, cmode, stat, errmsg)
    class(netcdf_file), intent(out) :: this
    character(*,kind=c_char), intent(in) :: path
    integer(c_int), intent(in) :: cmode
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_create(path//c_null_char, cmode, this%ncid)
    call handle_status(status, 'NETCDF_FILE%CREATE', stat, errmsg)
    status = nc_enddef(this%ncid)
    call handle_status(status, 'NETCDF_FILE%CREATE', stat, errmsg)
  end subroutine create

  subroutine open(this, path, omode, stat, errmsg)
    class(netcdf_file), intent(out) :: this
    character(*,kind=c_char), intent(in) :: path
    integer(c_int), intent(in) :: omode
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_open(path//c_null_char, omode, this%ncid)
    call handle_status(status, 'NETCDF_FILE%OPEN', stat, errmsg)
  end subroutine open

  subroutine close(this, stat, errmsg)
    class(netcdf_file), intent(inout) :: this
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_close(this%ncid)
    call handle_status(status, 'NETCDF_FILE%CLOSE', stat, errmsg)
    this%ncid = -1
  end subroutine close

  subroutine redef(this, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_redef(this%ncid)
    call handle_status(status, 'NETCDF_FILE%REDEF', stat, errmsg)
  end subroutine redef

  subroutine enddef(this, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_enddef(this%ncid)
    call handle_status(status, 'NETCDF_FILE%ENDDEF', stat, errmsg)
  end subroutine enddef

  subroutine def_dim(this, name, len, dimid, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    character(*,kind=c_char), intent(in) :: name
    integer(c_size_t), intent(in) :: len
    integer(c_int), intent(out) :: dimid
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_def_dim(this%ncid, name//c_null_char, len, dimid)
    call handle_status(status, 'NETCDF_FILE%DEF_DIM', stat, errmsg)
  end subroutine def_dim

  subroutine inq_dimid(this, name, dimid, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    character(*,kind=c_char), intent(in) :: name
    integer(c_int), intent(out) :: dimid
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_inq_dimid(this%ncid, name//c_null_char, dimid)
    call handle_status(status, 'NETCDF_FILE%INQ_DIMID', stat, errmsg)
  end subroutine inq_dimid

  subroutine inq_dimlen(this, dimid, len, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: dimid
    integer(c_size_t), intent(out) :: len
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_inq_dimlen(this%ncid, dimid, len)
    call handle_status(status, 'NETCDF_FILE%INQ_DIMLEN', stat, errmsg)
  end subroutine inq_dimlen

  subroutine def_var_1(this, name, xtype, dimids, varid, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    character(*,kind=c_char), intent(in) :: name
    integer(c_int), intent(in) :: xtype
    integer(c_int), intent(in) :: dimids
    integer(c_int), intent(out) :: varid
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_def_var(this%ncid, name//c_null_char, xtype, 1_c_int, [dimids], varid)
    call handle_status(status, 'NETCDF_FILE%DEF_VAR', stat, errmsg)
  end subroutine def_var_1

  subroutine def_var_n(this, name, xtype, dimids, varid, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    character(*,kind=c_char), intent(in) :: name
    integer(c_int), intent(in) :: xtype
    integer(c_int), intent(in) :: dimids(:)
    integer(c_int), intent(out) :: varid
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_def_var(this%ncid, name//c_null_char, xtype, size(dimids,kind=c_int), &
        dimids(size(dimids):1:-1), varid)
    call handle_status(status, 'NETCDF_FILE%DEF_VAR', stat, errmsg)
  end subroutine def_var_n

  subroutine def_inq_varid(this, name, varid, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    character(*,kind=c_char), intent(in) :: name
    integer(c_int), intent(out) :: varid
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_inq_varid(this%ncid, name//c_null_char, varid)
    call handle_status(status, 'NETCDF_FILE%INQ_VARID', stat, errmsg)
  end subroutine def_inq_varid

  subroutine inq_varid(this, name, varid, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    character(*,kind=c_char), intent(in) :: name
    integer(c_int), intent(out) :: varid
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_inq_varid(this%ncid, name//c_null_char, varid)
    call handle_status(status, 'NETCDF_FILE%INQ_VARID', stat, errmsg)
  end subroutine inq_varid

  subroutine put_var_int32(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(int32), intent(in) :: var(*)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_put_var_int(this%ncid, varid, var)
    call handle_status(status, 'NETCDF_FILE%PUT_VAR_INT32', stat, errmsg)
  end subroutine put_var_int32

  subroutine put_var_int64(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(int64), intent(in) :: var(*)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_put_var_longlong(this%ncid, varid, var)
    call handle_status(status, 'NETCDF_FILE%PUT_VAR_INT64', stat, errmsg)
  end subroutine put_var_int64

  subroutine put_var_real32(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    real(real32), intent(in) :: var(*)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_put_var_float(this%ncid, varid, var)
    call handle_status(status, 'NETCDF_FILE%PUT_VAR_REAL32', stat, errmsg)
  end subroutine put_var_real32

  subroutine put_var_real64(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    real(real64), intent(in) :: var(*)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_put_var_double(this%ncid, varid, var)
    call handle_status(status, 'NETCDF_FILE%PUT_VAR_REAL64', stat, errmsg)
  end subroutine put_var_real64

  subroutine put_var_real64_2(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    real(real64), intent(in) :: var(:,:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    call put_var_real64(this, varid, var, stat, errmsg)
  end subroutine put_var_real64_2

  subroutine put_vara_int32(this, varid, start, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(c_size_t), intent(in) :: start
    integer(int32), intent(in) :: var(:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_put_vara_int(this%ncid, varid, &
        [start-1_c_size_t], [size(var,kind=c_size_t)], var)
    call handle_status(status, 'NETCDF_FILE%PUT_VARA_INT32', stat, errmsg)
  end subroutine put_vara_int32

  subroutine put_vara_int64(this, varid, start, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(c_size_t), intent(in) :: start
    integer(int64), intent(in) :: var(:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_put_vara_longlong(this%ncid, varid, &
        [start-1_c_size_t], [size(var,kind=c_size_t)], var)
    call handle_status(status, 'NETCDF_FILE%PUT_VARA_INT64', stat, errmsg)
  end subroutine put_vara_int64

  subroutine put_vara_real32(this, varid, start, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(c_size_t), intent(in) :: start
    real(real32), intent(in) :: var(:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_put_vara_float(this%ncid, varid, &
        [start-1_c_size_t], [size(var,kind=c_size_t)], var)
    call handle_status(status, 'NETCDF_FILE%PUT_VARA_REAL32', stat, errmsg)
  end subroutine put_vara_real32

  subroutine put_vara_real64(this, varid, start, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(c_size_t), intent(in) :: start
    real(real64), intent(in) :: var(:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_put_vara_double(this%ncid, varid, &
        [start-1_c_size_t], [size(var,kind=c_size_t)], var)
    call handle_status(status, 'NETCDF_FILE%PUT_VARA_REAL64', stat, errmsg)
  end subroutine put_vara_real64

  subroutine get_var_int32(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(int32), intent(out) :: var(*)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_get_var_int(this%ncid, varid, var)
    call handle_status(status, 'NETCDF_FILE%GET_VAR_INT32', stat, errmsg)
  end subroutine get_var_int32

  subroutine get_var_int64(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(int64), intent(out) :: var(*)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_get_var_longlong(this%ncid, varid, var)
    call handle_status(status, 'NETCDF_FILE%GET_VAR_INT64', stat, errmsg)
  end subroutine get_var_int64

  subroutine get_var_real32(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    real(real32), intent(out) :: var(*)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_get_var_float(this%ncid, varid, var)
    call handle_status(status, 'NETCDF_FILE%GET_VAR_REAL32', stat, errmsg)
  end subroutine get_var_real32

  subroutine get_var_real64(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    real(real64), intent(out) :: var(*)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_get_var_double(this%ncid, varid, var)
    call handle_status(status, 'NETCDF_FILE%GET_VAR_REAL64', stat, errmsg)
  end subroutine get_var_real64

  subroutine get_var_real64_2(this, varid, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    real(real64), intent(out) :: var(:,:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    call get_var_real64(this, varid, var, stat, errmsg)
  end subroutine get_var_real64_2

  subroutine get_vara_int32(this, varid, start, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(c_size_t), intent(in) :: start
    integer(int32), intent(out) :: var(:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_get_vara_int(this%ncid, varid, &
        [start-1_c_size_t], [size(var,kind=c_size_t)], var)
    call handle_status(status, 'NETCDF_FILE%GET_VARA_INT32', stat, errmsg)
  end subroutine get_vara_int32

  subroutine get_vara_int64(this, varid, start, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(c_size_t), intent(in) :: start
    integer(int64), intent(out) :: var(:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_get_vara_longlong(this%ncid, varid, &
        [start-1_c_size_t], [size(var,kind=c_size_t)], var)
    call handle_status(status, 'NETCDF_FILE%GET_VARA_INT64', stat, errmsg)
  end subroutine get_vara_int64

  subroutine get_vara_real32(this, varid, start, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(c_size_t), intent(in) :: start
    real(real32), intent(out) :: var(:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_get_vara_float(this%ncid, varid, &
        [start-1_c_size_t], [size(var,kind=c_size_t)], var)
    call handle_status(status, 'NETCDF_FILE%GET_VARA_REAL32', stat, errmsg)
  end subroutine get_vara_real32

  subroutine get_vara_real64(this, varid, start, var, stat, errmsg)
    class(netcdf_file), intent(in) :: this
    integer(c_int), intent(in) :: varid
    integer(c_size_t), intent(in) :: start
    real(real64), intent(out) :: var(:)
    integer(c_int), intent(out), optional :: stat
    character(:,kind=c_char), allocatable, intent(out), optional :: errmsg
    integer(c_int) :: status
    status = nc_get_vara_double(this%ncid, varid, &
        [start-1_c_size_t], [size(var,kind=c_size_t)], var)
    call handle_status(status, 'NETCDF_FILE%GET_VARA_REAL64', stat, errmsg)
  end subroutine get_vara_real64

  subroutine handle_status(status, prefix, stat, errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
#if defined(NAGFOR) | defined(NAG_COMPILER)
    use,intrinsic :: f90_unix, only: exit
#endif
    integer(c_int), intent(in) :: status
    character(*), intent(in) :: prefix
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    character(:,kind=c_char), pointer :: string
    if (present(stat)) stat = status
    if (status == NC_NOERR) return
    string => f_string_pointer(nc_strerror(status))
    if (present(stat)) then
      if (present(errmsg)) errmsg = string
    else
      write(error_unit,'(3a)') prefix, ': ', string
      call exit(1)
    end if
  end subroutine handle_status

  !! Auxiliary function that converts a C string pointer to a Fortran
  !! character pointer.  Taken from the MIT-licensed Petaca library,
  !! https://github.com/nncarlson/petaca

  function f_string_pointer(cptr) result(fptr)
    type(c_ptr), intent(in) :: cptr
    character(:,kind=c_char), pointer :: fptr
    interface ! to strlen from the standard C library
      function strlen(s) result(len) bind(c)
        import c_ptr, c_size_t
        type(c_ptr), value :: s
        integer(c_size_t) :: len
      end function
    end interface
    fptr => f_string_pointer_aux(cptr, strlen(cptr))
  contains
    function f_string_pointer_aux(cptr, len) result(fptr)
      type(c_ptr), intent(in) :: cptr
      integer(c_size_t), intent(in) :: len
      character(len,kind=c_char), pointer :: fptr
      call c_f_pointer(cptr, fptr)
    end function
  end function f_string_pointer

end module netcdf_file_type
