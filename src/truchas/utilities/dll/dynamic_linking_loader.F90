!!
!! DYNAMIC_LINKING_LOADER
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 10 Apr 2006; last revised 18 Apr 2006
!!
!! This experimental module provides access to the dynamic linking loader.
!! The POSIX C functions dlopen, dlclose, dlsym, and dlerror from the system
!! DLL library (libdl on linux) are not directly callable from Fortran 95,
!! so the complementary piece to this module is the Fortran-callable, C wrapper
!! functions defined in dlwrap.c.
!!
!! NOTE: When compiled on a platform with 64-bit addresses, the pre-processor
!! macro M64 must be defined (e.g. -DM64 on the compile commandline).
!!
!! PROGRAMMING INTERFACE
!!
!! Each subroutine has the optional intent-out arguments STAT and ERRMSG.
!! If STAT is present, it is assigned the value 0 if the subroutine completes
!! successfully, and a nonzero value if an error occurs.  In the latter case,
!! the character string ERRMSG, if present, is assigned the error string
!! returned by the underlying DL library.  If STAT is not present and an
!! error occurs, the error string is written to stderr and the program exits
!! with a nonzero status.
!!
!!  CALL DLL_OPEN (PATH, MODE, HANDLE [,STAT] [,ERRMSG])
!!    CHARACTER(LEN=*), INTENT(IN) :: PATH
!!    INTEGER, INTENT(IN) :: MODE
!!    INTEGER(C_PTR_KIND) :: HANDLE
!!    INTEGER, INTENT(OUT), OPTIONAL :: STAT
!!    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: ERRMSG
!!
!!    This call opens the shared library file specified by PATH and returns
!!    a handle for the library.  This handle is intended only to be passed
!!    back to DLL_SYMBOL and DLL_CLOSE.  This can be called multiple times
!!    for the same library; the same handle is returned.  The value supplied
!!    for MODE must be one of the module parameters RTLD_LAZY or RTLD_NOW,
!!    optionally or'ed with the module parameter RTLD_GLOBAL.  See the man
!!    page for dlopen(3) for a detailed description of the behavior.
!!    
!!  CALL DLL_CLOSE (HANDLE [,STAT] [,ERRMSG])
!!    INTEGER(C_PTR_KIND) :: HANDLE
!!    INTEGER, INTENT(OUT), OPTIONAL :: STAT
!!    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: ERRMSG
!!
!!    This call decrements the reference count on the specified shared library
!!    handle.  When the reference count reaches zero, the shared library is
!!    unloaded.  See the man page for dlclose(3) for more details.
!!
!!  CALL DLL_SYMBOL (HANDLE, SYMBOL, FUNPTR [,STAT] [,ERRMSG])
!!    INTEGER(C_PTR_KIND) :: HANDLE
!!    CHARACTER(LEN=*), INTENT(IN) :: SYMBOL
!!    INTEGER(C_PTR_KIND), INTENT(OUT) :: FUNPTR
!!    INTEGER, INTENT(OUT), OPTIONAL :: STAT
!!    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: ERRMSG
!!
!!    Given the handle of a shared library returned by DLL_OPEN and the name of
!!    a symbol in that library, this routine returns the address of the symbol
!!    in FUNPTR.  The returned address is suitable only for passing back to a
!!    C function (not included in this module) which can invoke the procedure
!!    at that address with appropriate arguments.  For more details see the man
!!    page for dlsym(3).
!!

#ifdef M64
# define _C_PTR_KIND selected_int_kind(18)
#else
# define _C_PTR_KIND selected_int_kind(9)
#endif

module dynamic_linking_loader

  implicit none
  private

  public :: dll_open, dll_close, dll_symbol

  integer, parameter, public :: c_ptr_kind = _C_PTR_KIND  ! integer kind to hold C pointer

  !! Linking mode flags; values from bits/dlfcn.h on Linux.
  integer, parameter, public :: RTLD_LAZY   = 1
  integer, parameter, public :: RTLD_NOW    = 2
  integer, parameter, public :: RTLD_LOCAL  = 0
  integer, parameter, public :: RTLD_GLOBAL = 256

  !! Special pseudo-handles; values from dlfcn.h on Linux.
  integer(c_ptr_kind), parameter, public :: RTLD_DEFAULT = 0_c_ptr_kind
  integer(c_ptr_kind), parameter, public :: RTLD_NEXT = -1_c_ptr_kind

  !! Interfaces to the C wrapper functions in dlwrap.c.
  interface
    function f_dlopen(filename, flag) result(handle)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: flag
      integer(_C_PTR_KIND) :: handle
      end function
    function f_dlsym(handle, symbol) result(addr)
      integer(_C_PTR_KIND) :: handle
      character(len=*), intent(in) :: symbol
      integer(_C_PTR_KIND) :: addr
      end function
    integer function f_dlclose(handle)
      integer(_C_PTR_KIND) :: handle
      end function
    function f_dlerror() result(error)
      integer(_C_PTR_KIND) :: error
      end function
  end interface

contains

  subroutine dll_open (path, mode, handle, stat, errmsg)

    character(*), intent(in) :: path
    integer, intent(in) :: mode
    integer(c_ptr_kind), intent(out) :: handle
    integer, intent(out), optional :: stat
    character(*), intent(out), optional :: errmsg

    if (present(stat)) stat = 0    
    handle = f_dlopen(trim(path)//char(0), mode)
    if (handle == 0) call error_handler ('DLL_OPEN', f_dlerror(), stat, errmsg)

  end subroutine dll_open


  subroutine dll_symbol (handle, symbol, funptr, stat, errmsg)

    integer(c_ptr_kind) :: handle
    character(*), intent(in) :: symbol
    integer(c_ptr_kind), intent(out) :: funptr
    integer, optional, intent(out) :: stat
    character(*), optional, intent(out) :: errmsg

    integer(c_ptr_kind) :: cptr

    if (present(stat)) stat = 0
    cptr = f_dlerror()  ! clear any existing error
    funptr = f_dlsym(handle, trim(symbol)//char(0))
    cptr = f_dlerror()
    if (cptr /= 0) call error_handler ('DLL_SYM', cptr, stat, errmsg)

  end subroutine dll_symbol


  subroutine dll_close (handle, stat, errmsg)

    integer(c_ptr_kind) :: handle
    integer, optional, intent(out) :: stat
    character(*), optional, intent(out) :: errmsg

    if (present(stat)) stat = 0
    if (f_dlclose(handle) /= 0) call error_handler ('DLL_CLOSE', f_dlerror(), stat, errmsg)

  end subroutine dll_close


  subroutine error_handler (proc, cptr, stat, errmsg)

#ifdef NAG_COMPILER
    use f90_unix, only: exit
#endif

    character(*), intent(in) :: proc
    integer(c_ptr_kind), intent(in) :: cptr
    integer, optional, intent(out) :: stat
    character(*), optional, intent(out) :: errmsg

    integer, parameter :: MAX_LEN = 1023
    character(MAX_LEN) :: error

    !! Function from dlwrap.c for getting at the C string.    
    interface
      subroutine c_string_to_f (cptr, c)
        integer(_C_PTR_KIND) :: cptr
        character(len=*) :: c
      end subroutine
    end interface

    call c_string_to_f (cptr, error)

    if (present(stat)) then
      stat = 1
      if (present(errmsg)) errmsg = error
    else
      write(0,'(3a)') proc, ': ', trim(error)
      call exit(1)
    end if

  end subroutine error_handler

end module dynamic_linking_loader
