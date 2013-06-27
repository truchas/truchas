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
    use edit_module,             only: Long_Edit_Bounding_Coords, &
                                       Long_Output_Dt_Multiplier, &
                                       Short_Output_Dt_Multiplier
    use input_utilities,         only: seek_to_namelist
    use interface_output_module, only: Int_Output_Dt_Multiplier
    use output_control,          only: Output_Dt, Output_T, precise_output
    use parallel_info_module,    only: p_info
    use timing_tree
    use output_control,          only: Output_Dt_Multiplier, retain_last_step
    use probe_output_module,     only: Probe_Output_Cycle_Multiplier
#ifdef USE_TBROOK
    use output_data_module,      only: XML_Data_Format
#endif
         
    integer, intent(in) :: lun

    ! Local Variables
    logical :: fatal, found
    integer :: ioerror
    character(128) :: fatal_error_string

    ! Define OUTPUTS namelist
    namelist /OUTPUTS/ Output_Dt_Multiplier,       &
#ifdef USE_TBROOK
                       XML_Data_Format,                &
#endif
                       Int_Output_Dt_Multiplier,       &
                       Long_Edit_Bounding_Coords,      &
                       Long_Output_Dt_Multiplier,      &
                       Output_Dt,                      &
                       Output_T,                       &
                       Short_Output_Dt_Multiplier,     &
                       Probe_Output_Cycle_Multiplier,  &
                       precise_output,                 &
                       retain_last_step

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
          read (lun, outputs, IOSTAT = ioerror)
          fatal_error_string = 'OUTPUTS_INPUT: Error reading OUTPUTS namelist'
          fatal = (ioerror/=0)
       end if
    end if IO_PE_ONLY

    ! Check read errors
    call TLS_fatal_if_any (fatal, fatal_error_string)

    ! Broadcast data if no errors so far
    call OUTPUTS_INPUT_PARALLEL ()

    ! Check namelist
    call OUTPUTS_CHECK (fatal)
    call TLS_fatal_if_any (fatal, 'OUTPUTS_INPUT: OUTPUTS Namelist input error')

  END SUBROUTINE OUTPUTS_INPUT

  SUBROUTINE OUTPUTS_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default OUTPUTS namelist.
    !
    !=======================================================================
    use kinds, only: r8
    use edit_module,             only: long_edit, Long_Edit_Bounding_Coords,  &
                                       Long_Edit_Cell_List, long_edit_cells,  &
                                       Long_Output_Dt_Multiplier, short_edit, &
                                       Short_Output_Dt_Multiplier
    use interface_output_module, only: Int_Output_Dt_Multiplier, interface_dump
    use output_control,          only: Output_Dt, Output_T, precise_output
    use timing_tree
    use output_control,          only: Output_Dt_Multiplier, retain_last_step
    use probe_output_module,     only: Probe_Output_Cycle_Multiplier
#ifdef USE_TBROOK
    use output_data_module,      only: XML_Data_Format
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Output & delta times
    Output_Dt = 0
    Output_T  = 0

    ! do output at exactly the specified time, deprecated
    precise_output = .false.

    ! Edit variables
    short_edit                 = .false.
    long_edit                  = .false.
    Long_Edit_Cell_List        = 0
    long_edit_cells            = 0
    Short_Output_Dt_Multiplier = 0
    Long_Output_Dt_Multiplier  = 0
    Long_Edit_Bounding_Coords(1,:,:) =  HUGE(1.0_r8)    ! Min coord is defaulted large positive
    Long_Edit_Bounding_Coords(2,:,:) = -HUGE(1.0_r8)    ! Max coord is defaulted large negative

    ! interface output
    interface_dump           = .false.
    Int_Output_Dt_Multiplier = 0

    ! User output
    retain_last_step          = .false.

    ! XML output, always produced
#ifdef USE_TBROOK
    XML_Data_Format           = 'binary'
#endif
    Output_Dt_Multiplier  = 1
    Probe_Output_Cycle_Multiplier = 1

  END SUBROUTINE OUTPUTS_DEFAULT

  SUBROUTINE OUTPUTS_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast outputs namelist data to all PE's.
    !
    !======================================================================
    use edit_module,             only: Long_Edit_Bounding_Coords,    &
                                       Long_Output_Dt_Multiplier,    &
                                       Short_Output_Dt_Multiplier
    use interface_output_module, only: Int_Output_Dt_Multiplier
    use output_control,          only: Output_Dt, Output_T, precise_output
    use parallel_info_module,    only: p_info
    use pgslib_module,           only: PGSLIB_BCAST
    use timing_tree
    use output_control,          only: Output_Dt_Multiplier, retain_last_step
    use probe_output_module,     only: Probe_Output_Cycle_Multiplier  
#ifdef USE_TBROOK
    use output_data_module,      only: XML_Data_Format
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    if (.NOT. p_info%UseGlobalServices) then
       call PGSLIB_BCAST (Int_Output_Dt_Multiplier)
       call PGSLIB_BCAST (Long_Edit_Bounding_Coords)
       call PGSLIB_BCAST (Long_Output_Dt_Multiplier)
       call PGSLIB_BCAST (Output_Dt)
       call PGSLIB_BCAST (Output_T)
       call PGSLIB_BCAST (precise_output)
       call PGSLIB_BCAST (Short_Output_Dt_Multiplier)
       call PGSLIB_BCAST (retain_last_step)
       call PGSLIB_BCAST (Output_Dt_Multiplier)
       call PGSLIB_BCAST (Probe_Output_Cycle_Multiplier)
#ifdef USE_TBROOK
       call PGSLIB_BCAST (XML_DATA_FORMAT)
#endif
    end if

  END SUBROUTINE OUTPUTS_INPUT_PARALLEL

  SUBROUTINE OUTPUTS_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check OUTPUTS namelist.
    !
    !=======================================================================
    use edit_module,             only: long_edit, Long_Output_Dt_Multiplier, &
                                       short_edit, Short_Output_Dt_Multiplier
    use interface_output_module, only: interface_dump, Int_Output_Dt_Multiplier
    use output_control,          only: nops, Output_Dt, Output_T
    use parameter_module,        only: mops
    use time_step_module,        only: t
    use utilities_module,        only: STRING_COMPARE
#ifdef USE_TBROOK
    use tbrook_module,           only: B_IFORM_BINARY, B_IFORM_XML
    use output_data_module,      only: XML_Data_Format, XML_iFormat, enable_tbrook_output
#endif

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

          if (Output_T(n) < 0) then
             write (message,2) n, Output_T(n)
2            format('Output time Output_T(',i2,') =',1pe10.3,' <= 0.0')
             call TLS_error (message)
             fatal = .true.
          end if

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
    interface_dump   = .false.
    long_edit        = .false.
    short_edit       = .false.

#ifdef USE_TBROOK
    if (enable_tbrook_output) then
    XML_IFORMAT=B_IFORM_BINARY
    if (XML_Data_Format /= 'binary') then
       XML_Data_Format = ADJUSTL(XML_Data_Format)
       call STRING_COMPARE (TRIM(XML_Data_Format), 'binary', strings_match)
       if (.not. strings_match) then
          call STRING_COMPARE (TRIM(XML_Data_Format), 'xml', strings_match)
          XML_IFORMAT=B_IFORM_XML
       end if
       if (.not. strings_match) then
          write (message,637) TRIM(XML_Data_Format)
637       format ('xml Data file format XML_Data_Format = ',a,' not valid! Valid formats are "binary", or "xml"')
          call TLS_error (message)
          fatal = .true.
          return
       end if
    end if
    end if
#endif

    ! Interface output flag.
    if (ANY(Int_Output_Dt_Multiplier      /= 0)) interface_dump = .true.

    ! Long edit output flag.
    if (ANY(Long_Output_Dt_Multiplier     /= 0)) long_edit      = .true.

    ! Short edit output flag.
    if (ANY(Short_Output_Dt_Multiplier    /= 0)) short_edit     = .true.

  END SUBROUTINE OUTPUTS_CHECK

END MODULE OUTPUTS_INPUT_MODULE
