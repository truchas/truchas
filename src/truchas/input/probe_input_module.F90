!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PROBE_INPUT_MODULE 
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures for the input of probes.
  !
  !   Public Interface:
  !
  !     * call PROBE_INPUT ()
  !
  !       Default, read, check, and broadcast variables associated
  !       with the PROBE namelist.
  !
  ! Contains: PROBE_CHECK
  !           PROBE_DEFAULT
  !           PROBE_INPUT
  !           PROBE_INPUT_PARALLEL
  !
  ! Author(s): Sharen Cummins (scummins@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  ! Public Subroutines
  public :: PROBE_INPUT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE PROBE_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check for obvious errors in the PROBE namelist input variables.
    !
    !=======================================================================

    use probe_data_module,      only: probe_name, probe_description, &
                                      probe_coords, probe_coords_scale
    use mesh_input_module,      only: coordinate_scale_factor
    use input_utilities,        only: NULL_R
    use parameter_module,       only: nprobes
 
    ! Arguments
    logical, intent(INOUT)  :: fatal

    ! Local Variables
    integer :: i
    character(256) :: myName 
    character(128) :: message

   ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    fatal = .false.

    write (message, 10) nprobes
10  format(9x,'Identified ',i0,' probe(s)')
    call TLS_info ('')
    call TLS_info (message)

    do i = 1, nprobes

      if (probe_name(i) == 'Unnamed') then
        write(myName,1) 'PROBE_',i 
1       FORMAT(a,i2.2) 
        probe_name(i) = myName
      end if

      if (probe_description(i) == 'None') then
        write(myName,2) 'DESCRIPTION FOR PROBE_',i 
2       FORMAT(a,i2.2) 
        probe_description(i) = myName
      end if


      if ((probe_coords(1, i)  == NULL_R) .or. &
           (probe_coords(2, i) == NULL_R) .or. &
           (probe_coords(3, i) == NULL_R)) then
         write (message, 30) i
30       format ('probe coordinates in PROBE namelist ',i2,' is invalid')
         call TLS_error (message)
         fatal = .true.
         exit 
      end if

      if (coordinate_scale_factor /= 1 .and. probe_coords_scale(i) == NULL_R) then
         write (message, 40) i
40       format ('Mesh coordinate_scale_factor is specified, so probe_coords_scale factor reqd in PROBE namelist ',i2)
         call TLS_error (message)
         fatal = .true.
         exit 
      end if

      if (coordinate_scale_factor == 1 .and. probe_coords_scale(i) == NULL_R) then
         probe_coords_scale(i) = 1
      end if

   end do

 END SUBROUTINE PROBE_CHECK
  
 SUBROUTINE PROBE_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default variables in and related to the PROBE Namelist.
    !
    !=======================================================================
     use probe_data_module,   only: probe_name, probe_description, &
                                    probe_coords, probe_coords_scale
     use input_utilities,     only: NULL_R

     ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

     ! Default probe variables
     probe_name         = 'Unnamed'
     probe_description  = 'None'
     probe_coords       = NULL_R 
     probe_coords_scale = NULL_R

  END SUBROUTINE PROBE_DEFAULT

  SUBROUTINE PROBE_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read PROBE namelist. The variables set in this routine are listed
    !   in the PROBE namelist. If there is no PROBE namelist in the input deck,
    !   then defaults are used. If there is an error in the PROBE namelist,
    !   then execution terminates.
    !
    !=======================================================================

    use probe_data_module,      only: probe_name, probe_description, &
                                      probe_coords, probe_coords_scale
    use input_utilities,        only: seek_to_namelist, NULL_I, NULL_R
    use parallel_info_module,   only: p_info
    use parameter_module,       only: string_dim, string_len, nprobes
    use pgslib_module,          only: PGSLIB_BCAST

    integer, intent(in) :: lun

    ! Local Variables
    character(string_len), dimension(string_dim) :: Fatal_Error_String
    logical :: fatal, read_done, found
    integer :: probe_number
    integer :: ioerror
    character(128) :: message

    ! Define PROBE  Namelist
    namelist /PROBE/  &
                       probe_name,             &
                       probe_description,      &
                       probe_coords_scale,     &
                       probe_coords
   
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! fatal flag will indicate error
    fatal     = .false.
    read_done = .false.

    ! Default namelist
    call PROBE_DEFAULT ()

    ! Read Notice
    call TLS_info ('')
    call TLS_info (' Reading PROBE namelists ...')

    ! Input is done only on the IO PE.  Results of read are broadcast to all PE's
    ! For initialization, need to rewind input file.  Only do this once
    REWIND_IO_PE_ONLY: if (p_info%IOP) then
       rewind lun
    end if REWIND_IO_PE_ONLY

    ! Read namelist (reads multiple namelist block)
    ! This loops until either end of namelist or a fatal error is encountered

    nprobes = 0

    READ_PROBE: do

       ! Re-initialize critical input values

       probe_number              = NULL_I

       probe_name(0)             = 'Unnamed'
       probe_description(0)      = 'None'
       probe_coords(:,0)         = NULL_R
       probe_coords_scale(0)     = NULL_R

       ! Read next PROBE
       READ_IO_PE_ONLY: if (p_info%IOP) then

          call seek_to_namelist (lun, 'PROBE', found)
          read_done = .not.found
          if (read_done) then
             write (message, 5) 
5            format (9x,'This is the last PROBE Namelist.')
             call TLS_info (message)
          else
             ioerror = 0
             ! Inform usr that a PROBE  namelist is being read
             write (message, 10) nprobes  + 1
10           format ('         Reading PROBE namelist #',i2,' ...')
             call TLS_info ('')
             call TLS_info (message)

             ! Read the namelist
             read (lun, PROBE, IOSTAT = ioerror)
          end if

       end if READ_IO_PE_ONLY

       call PGSLIB_BCAST (read_done)

       if (read_done) exit READ_PROBE

       call PGSLIB_BCAST (ioerror)

       if (ioerror /= 0) then
          fatal = .true.
          exit READ_PROBE
       end if


       ! Broadcast data just read to all PE's. Also broadcast the
       ! flags read_done and ioerror. If we didn't actually do a
       ! read, then this bcast is unnecessary (since default values
       ! are bcast). The code is cleaner, though, so do the extra work.
       ! Within PROBE_input_parallel, only the 0th element of each
       ! array is broadcast, since that is the only element just read.

    call PROBE_INPUT_PARALLEL 

       ! New abbreviated input w/o indexing...temp, do in
       ! separate routine with checking & error msgs, etc. later

       nprobes  = nprobes + 1
       PROBE_ASSIGNED: if (probe_number == NULL_I) then
          probe_number  = nprobes 
          write (message, 15) nprobes , nprobes 
15        format (9x,'INFO: Assigned probe_number ',i0,' to PROBE namelist ',i0)
          call TLS_info (message)
       end if PROBE_ASSIGNED

!       ! Put input variables in proper slot

       probe_name(nprobes)          = probe_name(0)
       probe_description(nprobes)   = probe_description(0)
       probe_coords(:,nprobes)      = probe_coords(:,0)
       probe_coords_scale(nprobes)  = probe_coords_scale(0)          
         
       if (fatal) exit

    end do READ_PROBE


    ! If we had a fatal error, abort now
    call TLS_fatal_if_any (fatal, 'PROBE_INPUT: PROBE_FUNCTION input error')

    ! Check PROBE namelist
    call PROBE_CHECK (fatal)
    call TLS_fatal_if_any (fatal, 'PROBE_INPUT: PROBE_FUNCTION input error')

  END SUBROUTINE PROBE_INPUT


  SUBROUTINE PROBE_INPUT_PARALLEL ()

    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast variables in the PROBE_FUNCTION namelist to all processors.
    !
    !======================================================================

    use probe_data_module,   only: probe_name, probe_description, &
                                   probe_coords, probe_coords_scale

    use pgslib_module,        only: PGSLIB_BCAST

    ! Argument List

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

       call PGSLIB_BCAST (probe_name(0))
       call PGSLIB_BCAST (probe_description(0))
       call PGSLIB_BCAST (probe_coords(:,0))
       call PGSLIB_BCAST (probe_coords_scale(0))

  END SUBROUTINE PROBE_INPUT_PARALLEL

END MODULE PROBE_INPUT_MODULE
