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
    use parallel_communication,  only: is_IOP, broadcast
    use output_control,          only: Output_Dt_Multiplier
    use output_control,          only: part, part_path, write_mesh_partition
         
    integer, intent(in) :: lun

    ! Local Variables
    logical :: fatal, found
    integer :: ios
    character(128) :: iom
    character(128) :: fatal_error_string
    integer :: move_block_ids(32)
    character(32) :: move_toolpath_name

    ! Define OUTPUTS namelist
    namelist /OUTPUTS/ Output_Dt_Multiplier,       &
                       Output_Dt,                      &
                       Output_T,                       &
                       Short_Output_Dt_Multiplier,     &
                       move_block_ids, move_toolpath_name, write_mesh_partition

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize fatal flag
    fatal = .false.

    ! Read notice
    call TLS_info ('Reading OUTPUTS namelist ...')

    ! Default namelist
    call OUTPUTS_DEFAULT ()

    ! Read on IO PE only
    IO_PE_ONLY: if (is_IOP) then
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
          write_mesh_partition = .false.
          read(lun,nml=outputs,iostat=ios,iomsg=iom)
          fatal_error_string = 'error reading OUTPUTS namelist: ' // trim(iom)
          fatal = (ios/=0)
       end if
    end if IO_PE_ONLY

    ! Check read errors
    call TLS_fatal_if_any (fatal, fatal_error_string)

    ! Broadcast data if no errors so far
    call OUTPUTS_INPUT_PARALLEL ()
    call broadcast (move_block_ids)
    call broadcast (move_toolpath_name)
    call broadcast (write_mesh_partition)

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
      block
        use toolpath_driver, only: alloc_toolpath
        integer :: stat
        character(:), allocatable :: errmsg
        call alloc_toolpath(part_path, move_toolpath_name, stat, errmsg)
        if (stat /= 0) call TLS_fatal('unable to create toolpath: ' // errmsg)
      end block
    end if

  END SUBROUTINE OUTPUTS_INPUT

  SUBROUTINE OUTPUTS_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default OUTPUTS namelist.
    !
    !=======================================================================
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
    use parallel_communication,  only: broadcast
    use output_control,          only: Output_Dt_Multiplier

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    call broadcast (Output_Dt)
    call broadcast (Output_T)
    call broadcast (Short_Output_Dt_Multiplier)
    call broadcast (Output_Dt_Multiplier)

  END SUBROUTINE OUTPUTS_INPUT_PARALLEL

  SUBROUTINE OUTPUTS_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check OUTPUTS namelist.
    !
    !=======================================================================
    use edit_module,             only: short_edit, Short_Output_Dt_Multiplier
    use output_control,          only: mops, nops, Output_Dt, Output_T
    use time_step_module,        only: t

    ! Arguments
    logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: n
    logical :: dt_okay
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize fatal flag
    fatal = .false.

    ! Check number of output times and output delta times.
    dt_okay = .false.
    OUTPUT_DT_CHECK: do n = mops,1,-1
       ! Found last positive dt, exit.
       if (Output_Dt(n) > 0) then
          dt_okay = .true.
          exit OUTPUT_DT_CHECK
       end if
    end do OUTPUT_DT_CHECK
    nops = n

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
