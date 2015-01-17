MODULE PGSLib_Init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Provide the routines needed for initialization and finalization.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use PGSLib_Instrument
  use PGSLib_Stats
  USE PGSLib_Type_MODULE
  USE PGSLib_Misc_Utility
  use pgslib_c_binding
  
  IMPLICIT NONE
  SAVE
  PUBLIC:: PGSLib_Initialize, PGSLib_Finalize

  ! maximum length of a token on the command line
  ! *MUST* agree with what''s in 
  integer, parameter, public :: PGSLIB_CL_MAX_TOKEN_LENGTH = 1023

  ! $Id: pgslib_init.F,v 1.2 2002/09/12 20:52:34 lally Exp $
  ! Variables for use by source which USEs this module
  integer (PGSLib_Int_TYPE), PRIVATE, PARAMETER :: FILE_NAME_LENGTH = 2048

  ! Interfaces for routines in this module.
  ! The reason for these generic interfaces is so that the user needn''t use the _F
  ! post-fix for the calls.

  INTERFACE PGSLib_Initialize
     MODULE PROCEDURE PGSLib_Initialize_F
  END INTERFACE!PGSLib_Initialize

  INTERFACE PGSLib_Finalize
     MODULE PROCEDURE PGSLib_Finalize_F
  END INTERFACE!PGSLib_Finalize

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE PGSLib_Initialize_F(IO_PE, FILE_PREFIX, FILE_PER_PE, argv)
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE,   ONLY : PGSLib_PEInfo
    IMPLICIT NONE
    INTEGER (PGSLib_INT_Type), OPTIONAL, INTENT(IN) ::IO_PE
    character (LEN=*), intent(IN   ), OPTIONAL :: FILE_PREFIX
    logical (PGSLib_Log_Type), OPTIONAL, INTENT(IN) :: FILE_PER_PE
    character (LEN=*), dimension(:), pointer, optional :: argv

    ! Local variables
    INTEGER (PGSLib_Int_TYPE) :: nPE, thisPE, IO_ROOT_PE, lf
    integer (PGSLib_Int_Type) :: FILE_PER_PE_C
    character (LEN=FILE_NAME_LENGTH) :: FILE_PREFIX_C
    integer :: argc
    integer :: i

    ! If FILE_PREFIX provided, use that, otherwise use default
    if (PRESENT(FILE_PREFIX)) then
       if (LEN_TRIM(FILE_PREFIX) > FILE_NAME_LENGTH-5) then
          PRINT *, 'ERROR: FILE_PREFIX too long in PGSLib_Initialize'
       else
          lf = LEN_TRIM(FILE_PREFIX)
          FILE_PREFIX_C(1:lf) = FILE_PREFIX(1:lf)
          FILE_PREFIX_C(lf+1:) = ACHAR(0)
       end if
    else
       FILE_PREFIX_C = "PGSLib"
       lf = LEN_TRIM(FILE_PREFIX_C)
       FILE_PREFIX_C(lf+1:) = ACHAR(0)
    end if

    if (PRESENT(IO_PE)) then
       IO_ROOT_PE = IO_PE - 1
    else 
       IO_ROOT_PE = 0
    end if
       
    if (PRESENT(FILE_PER_PE)) then
       if (FILE_PER_PE) then
          FILE_PER_PE_C = PGSLib_TRUE
       else
          FILE_PER_PE_C = PGSLib_FALSE
       end if
    else
       ! Default is each PE gets its own file
       FILE_PER_PE_C = PGSLib_TRUE
    end if

    Call PGSLib_Initialize_C(nPE, thisPE, IO_ROOT_PE, FILE_PER_PE_C, FILE_PREFIX_C)

    ! The only way we can count on to get the argument vector back is one
    ! character string at a time.  That makes this kind of gruesome, but
    ! reasonably failsafe (if this fails, the mpich fortran initialization
    ! fails as well).  If we try and use an array of character
    ! variables, we end up with an array descriptor on the C side of
    ! things.  PGSLib_Initialize_C (above) created argc and argv, and cached
    ! the values for later.  If we got an argument to pass back the command
    ! line (argv), we now fill it up from the cached data.

    ! if we got passed a cptr, fill it up with the command line arguments
    if (PRESENT(argv)) then

       ! get argc from the cached data
       call pgslib_get_argc(argc)

       allocate(argv(0:argc-1))

       ! get argv values from the cached data
       do i = 0, argc - 1
          call pgslib_get_argv(i, LEN(argv), argv(i))
       end do
    end if

    ! tell C we''re finished with the command line (delete the cached data)
    call pgslib_cleanup_argv ()

    ! FORTRAN Uses 1 based numbering
    PGSLib_PEInfo%nPE = nPE
    PGSLib_PEInfo%thisPE = thisPE + 1
    PGSLib_PEInfo%IO_ROOT_PE = IO_ROOT_PE + 1

    ! Initialize the internal timers
    call Initialize_Instrument_Array()
    
    RETURN
  END SUBROUTINE PGSLib_Initialize_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE PGSLib_Finalize_F()
    USE PGSLib_Type_MODULE
    IMPLICIT NONE

    Call PGSLib_Finalize_C()

    RETURN
  END SUBROUTINE PGSLib_Finalize_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PGSLib_Init

#ifdef NAG_COMPILER
!! utility-c.c looks for Fortran iargc and getarg.
!! Some compilers provide these as intrinsics.  NAG provides them as
!! module procedures which are not directly callable from C.
integer function nag_iargc ()
  use f90_unix
  nag_iargc = iargc()
end function nag_iargc
subroutine nag_getarg (k, arg)
  use f90_unix
  integer, intent(in) :: k
  character(len=*), intent(out) :: arg
  call getarg (k, arg)
end subroutine nag_getarg
#endif
