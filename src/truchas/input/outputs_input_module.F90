!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE OUTPUTS_INPUT_MODULE
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures for the input of outputs parameters
  !
  !   Public Interface:
  !
  !     * call OUTPUTS_INPUT ()
  !
  !       Defaults, reads, checks, and broadcasts input variables
  !       in the OUTPUTS namelist.
  !
  ! Contains: OUTPUTS_CHECK
  !           OUTPUTS_DEFAULT
  !           OUTPUTS_INPUT
  !           OUTPUTS_INPUT_PARALLEL
  !
  ! Author(s): Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !=======================================================================
  use truchas_logging_services
  Implicit None
  Private

  ! public procedures
  Public :: OUTPUTS_INPUT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE OUTPUTS_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read OUTPUTS namelist.
    !
    !=======================================================================
    use edit_module,             only: Short_Output_Dt_Multiplier
    use input_utilities,         only: seek_to_namelist, NULL_I, NULL_C
    use output_control,          only: Output_Dt, Output_T
    use parallel_info_module,    only: p_info
    use pgslib_module,           only: PGSLIB_BCAST
    use output_control,          only: Output_Dt_Multiplier
    use output_control,          only: part, part_path
    use toolpath_table,          only: toolpath_ptr
         
    integer, intent(in) :: lun

    ! Local Variables
    logical :: fatal, found
    integer :: ioerror
    character(128) :: fatal_error_string
    integer :: move_block_ids(32)
    character(32) :: move_toolpath_name

    ! Define OUTPUTS namelist
    namelist /OUTPUTS/ Output_Dt_Multiplier,       &
                       Output_Dt,                      &
                       Output_T,                       &
                       Short_Output_Dt_Multiplier,     &
                       move_block_ids, move_toolpath_name

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize fatal flag
    fatal = .false.

    ! Read notice
    call TLS_info ('')
    call TLS_info (' Reading OUTPUTS Namelist ...')

    ! Default namelist
    call OUTPUTS_DEFAULT ()

    ! Read on IO PE only
    IO_PE_ONLY: if (p_info%IOP) then
       ! Find namelist
       rewind lun
       call seek_to_namelist (lun, 'OUTPUTS', found)
       fatal = .not.found

       ! Namelist not found
       fatal_error_string = 'OUTPUTS_INPUT: OUTPUTS namelist not found'

       if (.not. fatal) then
          ! Read namelist
          move_block_ids = NULL_I
          move_toolpath_name = NULL_C
          read (lun, outputs, IOSTAT = ioerror)
          fatal_error_string = 'OUTPUTS_INPUT: Error reading OUTPUTS namelist'
          fatal = (ioerror/=0)
       end if
    end if IO_PE_ONLY

    ! Check read errors
    call TLS_fatal_if_any (fatal, fatal_error_string)

    ! Broadcast data if no errors so far
    call OUTPUTS_INPUT_PARALLEL ()
    call PGSLIB_BCAST (move_block_ids)
    call PGSLIB_BCAST (move_toolpath_name)

    ! Check namelist
    call OUTPUTS_CHECK (fatal)
    call TLS_fatal_if_any (fatal, 'OUTPUTS_INPUT: OUTPUTS Namelist input error')

    if (any(move_block_ids /= NULL_I)) then
      if (move_toolpath_name == NULL_C) call TLS_fatal('MOVE_TOOLPATH_NAME not specified')
    else if (move_toolpath_name /= NULL_C) then
      call TLS_fatal('MOVE_BLOCK_IDS not specified')
    end if

    allocate(part(count(move_block_ids/=NULL_I)))
    part = pack(move_block_ids,mask=move_block_ids/=NULL_I)
    if (move_toolpath_name /= NULL_C) then
      part_path => toolpath_ptr(move_toolpath_name)
      if (.not.associated(part_path)) &
          call TLS_fatal('no such toolpath: ' // trim(move_toolpath_name))
    end if

  END SUBROUTINE OUTPUTS_INPUT

  SUBROUTINE OUTPUTS_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default OUTPUTS namelist.
    !
    !=======================================================================
    use kinds, only: r8
    use edit_module,             only: short_edit, Short_Output_Dt_Multiplier
    use output_control,          only: Output_Dt, Output_T
    use output_control,          only: Output_Dt_Multiplier

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Output & delta times
    Output_Dt = 0
    Output_T  = 0

    ! Edit variables
    short_edit                 = .false.
    Short_Output_Dt_Multiplier = 0

    ! User output
    Output_Dt_Multiplier  = 1

  END SUBROUTINE OUTPUTS_DEFAULT

  SUBROUTINE OUTPUTS_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast outputs namelist data to all PE's.
    !
    !======================================================================
    use edit_module,             only: Short_Output_Dt_Multiplier
    use output_control,          only: Output_Dt, Output_T
    use parallel_info_module,    only: p_info
    use pgslib_module,           only: PGSLIB_BCAST
    use output_control,          only: Output_Dt_Multiplier

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    if (.NOT. p_info%UseGlobalServices) then
       call PGSLIB_BCAST (Output_Dt)
       call PGSLIB_BCAST (Output_T)
       call PGSLIB_BCAST (Short_Output_Dt_Multiplier)
       call PGSLIB_BCAST (Output_Dt_Multiplier)
    end if

  END SUBROUTINE OUTPUTS_INPUT_PARALLEL

  SUBROUTINE OUTPUTS_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check OUTPUTS namelist.
    !
    !=======================================================================
    use edit_module,             only: short_edit, Short_Output_Dt_Multiplier
    use output_control,          only: nops, Output_Dt, Output_T
    use parameter_module,        only: mops
    use time_step_module,        only: t
    use utilities_module,        only: STRING_COMPARE

    ! Arguments
    logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: n
    logical :: dt_okay, strings_match
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize fatal flag
    fatal = .false.

    ! Check number of output times and output delta times.
    dt_okay = .false.
    OUTPUT_DT_CHECK: do nops = mops,1,-1
       ! Found last positive dt, exit.
       if (Output_Dt(nops) > 0) then
          dt_okay = .true.
          exit OUTPUT_DT_CHECK
       end if
    end do OUTPUT_DT_CHECK

    ! Check output times
    if (.not.dt_okay) then

       ! Fatal: No output times specified.
       call TLS_error ('No output times specified')
       fatal   = .true.

    else

       OUTPUT_T_CHECK: do n = 1,nops

          if (Output_T(n) >= Output_T(n+1)) then
             write (message,3) n, Output_T(n), n+1, Output_T(n+1)
3            format('Output time Output_T(',i2,') =',1pe10.3, &
                    ' >= Output_T(',i2,') =',1pe10.3)
             call TLS_error (message)
             fatal = .true.
          end if

          if (Output_Dt(n) <= 0) then
             write (message,4) n, Output_Dt(n)
4            format('Output delta time Output_Dt(',i2,') =',1pe10.3,' <= 0.0')
             call TLS_error (message)
             fatal = .true.
          end if

       end do OUTPUT_T_CHECK

       ! Check last output time
       if (Output_T(nops+1) < t) then
          write (message,5) nops+1, Output_T(nops+1), t
5         format('The last output time Output_T(',i2,') =',es10.3,'is less than current problem time t =',es10.3)
          call TLS_error (message)
          fatal = .true.
       end if

    end if

    ! Initialize output flags
    short_edit       = .false.

    ! Short edit output flag.
    if (ANY(Short_Output_Dt_Multiplier    /= 0)) short_edit     = .true.

  END SUBROUTINE OUTPUTS_CHECK

END MODULE OUTPUTS_INPUT_MODULE
