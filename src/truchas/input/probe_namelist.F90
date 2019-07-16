!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The probe_namelist module contains the procedure to read the probe namelists.
! The format for these namelists is:
!
! &PROBE
!   Name  = 'string'
!   Coord = <3 reals>
!   Data  = 'Temperature', 'Pressure', 'Density', OR 'Velocity'
!   Description = 'string' (optional)
!   Coord_Scale_Factor = <1 real> (optional)
! /
! 
! The probe data will be output into a file in the output directory named with 
! the Name string, and an appended '.dat'.

module probe_namelist

  use parameter_list_type
  implicit none
  private

  ! Public function.

  public :: read_probe_namelists

  ! Output variable. This is set from the name lists by the read_probe_namelists
  ! call and accessed via use association wherever needed.

  type(parameter_list), public :: params

contains

  !-----------------------------------------------------------------------------
  ! Procedure read_probe_namelists:
  !
  ! Read any PROBE namelists that are found and set the global module variable
  ! params to include their information.
  !-----------------------------------------------------------------------------

  subroutine read_probe_namelists (lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use truchas_logging_services    ! Includes TLS_* procedures.

    ! Input variable.

    integer, intent(in) :: lun   ! Logical unit number for input file.

    ! Output variable.

    ! The global module variable params is updated with info from any
    ! PROBE namelists that are found.

    ! Internal variables.

    type(parameter_list), pointer :: sublist  ! params sublist pointer. 
    character(80) :: iom                      ! IO message.
    integer :: ios                            ! IO status number.
    logical :: found                          ! True if a given probe is found.
    ! Namelist variables:
    character(80) :: name                     ! Name for the probe.
    character(80) :: description              ! Description of the probe.
    character(80) :: data                     ! Field name to be output.
    real(r8), dimension(3) :: coord           ! Spatial coordinates of probe.
    real(r8) :: coord_scale_factor            ! Coordinate scale factor.
    integer :: nprobes                        ! Probe counter / total number.
    character(:), allocatable :: label        ! Error output label.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Define the PROBE Namelist.

    namelist /PROBE/ &
                     name,               &
                     description,        &
                     coord_scale_factor, &
                     coord,              &
                     data

    ! Rewind the input file on the IO PE.

    if (is_IOP) rewind (lun)

    ! Read PROBE namelists -- reads any and all PROBE namelist blocks.

    call TLS_info ('')
    call TLS_info ('Reading PROBE namelists ...')
    nprobes = 0
    do

      ! Search for a PROBE namelist location.

      if (is_IOP) call seek_to_namelist (lun, 'PROBE', found, iostat=ios)
      call broadcast (ios)
      if (ios /= 0) &
        call TLS_fatal ('Error reading input file: iostat=' // i_to_c(ios))
      call broadcast (found)
      if (.not.found) exit

      ! PROBE namelist found.

      nprobes = nprobes + 1
      label = 'PROBE[' // i_to_c(nprobes) // ']'

      ! Set default probe variables.

      name               = NULL_C
      description        = NULL_C
      coord              = NULL_R 
      coord_scale_factor = NULL_R
      data               = NULL_C

      ! Read the PROBE namelist.

      if (is_IOP) read (lun, nml=PROBE, iostat=ios, iomsg=iom)

      ! Error output and stop for namelist read problem.

      call broadcast (ios)
      if (ios /= 0) then
        call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))
      end if
      call TLS_info ('  Reading PROBE namelist: ' // trim(name) // ' ...')

      ! Broadcast the PROBE namelist variables.

      call broadcast (name)
      call broadcast (description)
      call broadcast (coord(:))
      call broadcast (coord_scale_factor)
      call broadcast (data)

      ! Error check for lack of field name.

      if (data == NULL_C) then
        data = 'Temperature'
        call TLS_info ('  WARNING: ' // trim(name) // &
                       ' probe does not have a field name. Setting equal')
        call TLS_info ('  to ' // trim(data) // ' and continuing.')
        call TLS_info ('')
      end if

      ! Error check for unique probe name. 

      if (name == NULL_C) then
        call TLS_fatal (label // ': NAME not specified')
      else if (params%is_sublist(name)) then
        call TLS_fatal (label // ': another PROBE namelist has this NAME: ' &
                        // trim(name))
      end if

      ! Put namelist input variables into params.

      sublist => params%sublist(trim(name))
      ! Required:
      call sublist%set ('name',                 trim(name))
      call sublist%set ('coord',                coord)
      call sublist%set ('data',                 trim(data))
      ! Optional:
      if (description /= NULL_C) &
        call sublist%set ('description',        trim(description))
      if (coord_scale_factor /= NULL_R) &
        call sublist%set ('coord-scale-factor', coord_scale_factor)

    end do

    ! Output total number of probes.

    call TLS_info ('  Total number of probes defined: ' // i_to_c(nprobes))

  end subroutine read_probe_namelists

end module probe_namelist
