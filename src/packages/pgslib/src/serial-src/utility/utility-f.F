! The subroutines in this directory are dummy stubs.
! They are replaced by C routines in the actual library

! $Id: utility-f.F,v 1.2 2002/09/12 20:52:35 lally Exp $

SUBROUTINE PGSLib_Initialize_C(nPE, thisPE, IO_ROOT_PE, FILE_PER_PE, FILE_PREFIX)
      IMPLICIT NONE
      INTEGER, intent(OUT)   :: nPE, thisPE
      INTEGER, intent(INOUT) :: IO_ROOT_PE
      INTEGER, intent(IN   ) :: FILE_PER_PE
      character (LEN=*), intent(IN   ) :: FILE_PREFIX

      ! File names are setup and used by output procedures.
      character (LEN=2048) :: F_OUT_NAME, F_ERR_NAME
      integer :: ERR_UNIT, OUT_UNIT
      COMMON /F_NAMES/F_OUT_NAME, F_ERR_NAME
      COMMON /F_UNITS/ERR_UNIT, OUT_UNIT

! In the C version, calls to MPI will determine values of nPE and thisPE.
! There may be some calls to open files on IO_ROOT_PE.
      ! Local variables
      logical, save :: been_initialized = .FALSE.
      integer       :: lf

      ! This routine may get called more than once, but can do init only once.
      IF (.NOT. been_initialized) THEN
         been_initialized = .TRUE.
         lf = LEN_TRIM(FILE_PREFIX)
         ! Drop the NULL at the end of the strings (only needed by C routines)
         F_OUT_NAME         = " "
         F_OUT_NAME(1:lf-1) = FILE_PREFIX(1:lf-1)
         F_OUT_NAME(lf:lf+3)= "-out"
         F_OUT_NAME         = TRIM(F_OUT_NAME)
         OUT_UNIT           = 29
         F_ERR_NAME         = ' '
         F_ERR_NAME(1:lf-1) = FILE_PREFIX(1:lf-1)
         F_ERR_NAME(lf:lf+3)= '-err'
         ERR_UNIT           = 30
      END IF

      ! NOTE: This routine must return valid values on all calls.
      nPE = 1
      thisPE = 0
      IO_ROOT_PE = 0

      RETURN
END SUBROUTINE PGSLib_Initialize_C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PGSLib_Finalize_C()
      IMPLICIT NONE

! In the C version, calls to MPI to finalize, and perhaps some calls to close output.

      RETURN
END SUBROUTINE PGSLib_Finalize_C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PGSLib_Error_C(Estring)
      IMPLICIT NONE
      CHARACTER (*):: Estring

      ! File names are setup and used by output procedures.
      character (LEN=2048) :: F_OUT_NAME, F_ERR_NAME
      integer :: ERR_UNIT, OUT_UNIT
      COMMON /F_NAMES/F_OUT_NAME, F_ERR_NAME
      COMMON /F_UNITS/ERR_UNIT, OUT_UNIT

      logical, save :: Err_Opened = .FALSE.
      integer       :: ioerror
! In the C version, the PE number is pre-pended to the string, then the whole
!          business is output to stderr (which may be a file).
! Since this is FORTRAN version, drop the NULL which PGSLib_Error_F
      ! tacked onto the end.

      if (.NOT. Err_Opened) then
         OPEN(UNIT=ERR_UNIT,FILE= TRIM(F_ERR_NAME), IOSTAT=ioerror)
         if (ioerror /= 0) then
            WRITE(*,*) 'Error opening error output file:',TRIM(F_ERR_NAME)
         end if
         Err_Opened = .TRUE.
      end if
      WRITE(ERR_UNIT,*) 'ERROR: ', EString(1:(len_trim(EString)-1))
      RETURN
END SUBROUTINE PGSLib_Error_C


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PGSLib_Output_C(Estring)
      IMPLICIT NONE
      CHARACTER (*):: Estring

      ! File names are setup and used by output procedures.
      character (LEN=2048) :: F_OUT_NAME, F_ERR_NAME
      integer :: ERR_UNIT, OUT_UNIT
      COMMON /F_NAMES/F_OUT_NAME, F_ERR_NAME
      COMMON /F_UNITS/ERR_UNIT, OUT_UNIT

      logical, save :: Out_Opened = .FALSE.
      integer       :: ioerror
! In the C version, the PE number is pre-pended to the string, then the whole
!          business is output to stdout (which may be a file).
! Since this is FORTRAN version, drop the NULL which PGSLib_Output_F
      ! tacked onto the end.

      if (.NOT. Out_Opened) then
         OPEN(UNIT=OUT_UNIT,FILE= TRIM(F_OUT_NAME), IOSTAT=ioerror)
         if (ioerror /= 0) then
            WRITE(*,*) 'Error opening standard output file:',TRIM(F_OUT_NAME)
         end if
         Out_Opened = .TRUE.
      end if
      WRITE(OUT_UNIT,*) 'Output: ', EString(1:(len_trim(EString)-1))
      RETURN
END SUBROUTINE PGSLib_Output_C


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PGSLib_Flush_Output_C()
      IMPLICIT NONE
      ! There is no Flush in fortran

      RETURN
    END SUBROUTINE PGSLib_Flush_Output_C


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PGSLib_Abort_C()

! In the C version, will try to bring all PEs down together.
!          For this version, just stop

      STOP
      END SUBROUTINE PGSLib_Abort_C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PGSLib_Barrier_C
        
! In the C version this calls MPI_Barrier.  Since this 
! is the serial version there is nothing to do.
        RETURN
      END SUBROUTINE PGSLib_Barrier_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PGSLib_Barrier_Time_C(BT)
      implicit none
      real :: BT
        
! In the C version this reports barrier time  Since this 
! is the serial version return 0
      BT = 0.0 

        RETURN
      END SUBROUTINE PGSLib_Barrier_Time_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PGSLib_SR_Time_C(T)
      implicit none
      real :: T
        
! In the C version this reports send/receive time  Since this 
! is the serial version return 0.
      T = 0.0 

        RETURN
      END SUBROUTINE PGSLib_SR_Time_C

! in the C version, these three routines are used to get the command
! line from C back to Fortran

      subroutine pgslib_get_argc (argc)
#ifdef NAG_COMPILER
         use f90_unix, only: iargc
#endif
         implicit none
         integer :: argc
#ifndef NAG_COMPILER
         integer :: iargc
#endif

         argc = iargc() + 1
         return
      end subroutine pgslib_get_argc

      subroutine pgslib_get_argv (i, string_l, string)
#ifdef NAG_COMPILER
         use f90_unix, only: getarg
#endif
         implicit none
         integer :: i
         integer :: string_l
         character(len=*) :: string

         call getarg (i, string)
         return
      end subroutine pgslib_get_argv

      subroutine pgslib_cleanup_argv
         implicit none
         return
      end subroutine pgslib_cleanup_argv
