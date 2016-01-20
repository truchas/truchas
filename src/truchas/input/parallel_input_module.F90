!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module PARALLEL_INPUT_MODULE
   !=======================================================================
   ! purpose:
   !
   !    specify parameters for parallel computation
   !
   !    Public Interface:
   !
   !       call PARALLEL_PARAMETERS_INPUT ()
   !
   !       defaults, reads, checks, and broadcasts input variables
   !       in the PARALLEL_PARAMETERS namelist
   !
   ! contains: PARALLEL_PARAMETERS_INPUT
   !           PARALLEL_PARAMETERS_CHECK
   !           PARALLEL_PARAMETERS_DEFAULT
   !           PARALLEL_PARAMETERS_PARALLEL
   !
   ! Author(s): Robert C. Ferrell (ferrell@cpca.com)
   !            Bryan Lally (lally@lanl.gov)
   !=======================================================================
   use parameter_module, only: string_len
   use mesh_parameter_module, only: ndim
   use truchas_logging_services
   implicit none
   private

   ! public procedures
   public :: PARALLEL_PARAMETERS_INPUT

   character(string_len), save :: partitioner
   integer, dimension(ndim), save :: processor_array

Contains

   subroutine PARALLEL_PARAMETERS_INPUT (lun)
      !=======================================================================
      ! purpose:
      !
      !   read PARALLEL_PARAMETERS namelist
      !=======================================================================
      use input_utilities, only: seek_to_namelist
      use parallel_info_module, only: p_info
      
      integer, intent(in) :: lun

      ! define PARALLEL_PARAMETERS namelist
      namelist /PARALLEL_PARAMETERS/ partitioner, processor_array

      ! local variables
      logical :: fatal
      logical :: found
      integer :: ioerror

      ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      fatal = .FALSE.

      ! inform user that the PARALLEL_PARAMETERS namelist is being read
      call TLS_info (' Reading PARALLEL_PARAMETERS Namelist ...')

      ! set  variables in this namelist to their defaults
      call PARALLEL_PARAMETERS_DEFAULT ()

      ! read the namelist only on the I/O processor
      IO_PE_ONLY: if (p_info%IOP) then

         ! find namelist
         rewind lun
         call seek_to_namelist (lun, 'PARALLEL_PARAMETERS', found)

         ! if we found the namelist, read it, otherwise warn
         if (found) then
            read (lun, NML = PARALLEL_PARAMETERS, IOSTAT = ioerror)
            fatal = (ioerror /= 0)
         else
            call TLS_info ('PARALLEL_PARAMETERS namelist not found; using defaults.')
         end if

      end if IO_PE_ONLY
    
      ! abort if there was a fatal error reading namelist
      call TLS_fatal_if_any (fatal, 'PARALLEL_PARAMETERS_INPUT: error reading PARALLEL_PARAMETERS namelist')

      ! broadcast all varitables in the PARALLEL_PARAMETERS namelist
      call PARALLEL_INPUT_PARALLEL ()

      ! check for fatal input errors
      call PARALLEL_PARAMETERS_CHECK (fatal)
    
      ! abort if error found
      call TLS_fatal_if_any (fatal, 'PARALLEL_PARAMETERS_INPUT: PARALLEL_PARAMETERS namelist input error')

   end subroutine PARALLEL_PARAMETERS_INPUT

   !----------------------------------------------------------------------------

   subroutine PARALLEL_PARAMETERS_DEFAULT ()
      !=======================================================================
      ! purpose:
      !
      !   default PARALLEL_PARAMETERS namelist
      !=======================================================================
      use partitioner_data, only: GET_PARTITIONER_DEFAULT, GET_PROCESSOR_ARRAY_DEFAULT

      ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      partitioner = GET_PARTITIONER_DEFAULT ()
      processor_array = GET_PROCESSOR_ARRAY_DEFAULT ()

   end subroutine PARALLEL_PARAMETERS_DEFAULT

   !----------------------------------------------------------------------------

   subroutine PARALLEL_INPUT_PARALLEL ()
      !=======================================================================
      ! purpose:
      !
      !   broadcast PARALLEL_PARAMETERS namelist variables
      !=======================================================================
      use pgslib_module, only: PGSLib_BCAST

      ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      
      call PGSLib_BCAST (partitioner)
      call PGSLib_BCAST (processor_array)

   end subroutine PARALLEL_INPUT_PARALLEL

   !----------------------------------------------------------------------------

   subroutine PARALLEL_PARAMETERS_CHECK (fatal)
      !=======================================================================
      ! purpose:
      !
      !   check that PARALLEL_PARAMETERS are consistent
      !=======================================================================
      use partitioner_data, only: SET_PARTITIONER, SET_PROCESSOR_ARRAY

      ! arguments
      logical, intent(INOUT) :: fatal

      ! local variables
      integer :: status
      character(64) :: message

      ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      fatal = .false.
      status = SET_PARTITIONER (partitioner)

      if (status /= 0) then
         fatal = .true.
         status=0
         call TLS_error ('unknown partitioner requested: ' // trim(partitioner))
      end if

      status = SET_PROCESSOR_ARRAY(processor_array)
      if (status /= 0) then
         fatal = .true.
         status=0
         call TLS_error ('invalid processor_array requested')
      end if

   end subroutine PARALLEL_PARAMETERS_CHECK

   !----------------------------------------------------------------------------

End Module PARALLEL_INPUT_MODULE  
