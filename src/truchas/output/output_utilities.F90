!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE OUTPUT_UTILITIES
   !----------------------------------------------------------------------------
   ! Purpose:
   !
   !   common output utilities
   !
   !   public interface:
   !
   !   call ANNOUNCE_FILE_WRITE (description, filename)
   !----------------------------------------------------------------------------

   implicit none
   private

   ! public procedures
   public :: ANNOUNCE_FILE_WRITE, ANNOUNCE

CONTAINS

   !----------------------------------------------------------------------------
   Subroutine ANNOUNCE_FILE_WRITE (description, filename)
      !-------------------------------------------------------------------------
      ! Purpose:
      !
      !    announce the various dumps
      !-------------------------------------------------------------------------
      use time_step_module, only: t
      use utilities_module, only: TIMESTAMP
      use truchas_logging_services

      character(*), intent(in) :: description, filename

      character(32) :: run_date
      character(256) :: message

      if (TLS_verbosity >= TLS_VERB_NOISY) then
        call TIMESTAMP (date_time=run_date)
        write(message,'(2a,es12.5,4a)') trim(description), ' written at t = ', t, &
                                        ' on ', run_date(5:22), ' to ', trim(filename)
        call TLS_info ('')
        call TLS_info (message)
      end if

    End Subroutine ANNOUNCE_FILE_WRITE

    SUBROUTINE ANNOUNCE (string)
      !---------------------------------------------------------------------------
      ! Purpose:
      !
      !   announce a particular processing phase
      !---------------------------------------------------------------------------
      use truchas_logging_services

      character(LEN=*) :: string

      call TLS_info ('')
      call TLS_info (repeat('=',80))
      call TLS_info ('')
      call TLS_info (string)
      !call TLS_info ('')

    END SUBROUTINE ANNOUNCE

End MODULE OUTPUT_UTILITIES
