!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module induction_coil_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use string_utilities, only: i_to_c
  use truchas_logging_services
  use parameter_list_type
  implicit none
  private

  public :: read_induction_coil_namelists

  integer, parameter :: MAXSV = 32

contains

  subroutine read_induction_coil_namelists(lun, params)

    use input_utilities, only: seek_to_namelist, NULL_I, NULL_R
    use parallel_communication, only: is_IOP, broadcast

    integer, intent(in) :: lun
    type(parameter_list), pointer :: params

    integer :: n, ios
    logical :: found
    character(128) :: iom
    character(:), allocatable :: label
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer :: num_loops
    real(r8) :: center(3), radius, length, current(MAXSV)
    namelist /induction_coil/ center, radius, length, num_loops, current

    call TLS_info('Reading INDUCTION_COIL namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all INDUCTION_COIL namelists have been read

      if (is_IOP) call seek_to_namelist(lun, 'induction_coil', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'INDUCTION_COIL[' // i_to_c(n) // ']'

      center = NULL_R
      radius = NULL_R
      length = NULL_R
      num_loops = NULL_I
      current = NULL_R

      if (is_IOP) read(lun,nml=induction_coil,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(center)
      call broadcast(radius)
      call broadcast(length)
      call broadcast(num_loops)
      call broadcast(current)

      plist => params%sublist('COIL' // i_to_c(n))

      if (num_loops == NULL_I) then
        call TLS_fatal(label // ': NUM_LOOPS must be assigned a value')
      else if (num_loops <= 0) then
        call TLS_fatal(label // ': NUM_LOOPS must be > 0')
      end if
      call plist%set('num-loops', num_loops)

      if (radius == NULL_R) then
        call TLS_fatal(label // ': RADIUS must be assigned a value')
      else if (radius <= 0.0_r8) then
        call TLS_fatal(label // ': RADIUS must be > 0.0')
      end if
      call plist%set('radius', radius)

      if (num_loops > 1) then
        if (length == NULL_R) then
          call TLS_fatal(label // ': LENGTH must be assigned a value')
        else if (length < 0.0_r8) then
          call TLS_fatal(label // ': LENGTH must be >= 0.0')
        end if
        call plist%set('length', length)
      end if

      if (any(center /= NULL_R)) then
        if (any(center == NULL_R)) call TLS_fatal(label // ': CENTER requires 3 values')
        call plist%set('center', center)
      end if

      call plist%set('current', pack(current, mask=(current/=NULL_R)))
    end do

    select case (n)
    case (0)
      call TLS_info('  none found')
    case (1)
      call TLS_info('  read 1 INDUCTION_COIL namelist')
    case default
      call TLS_info('  read ' // i_to_c(n) // ' INDUCTION_COIL namelists')
    end select

  end subroutine read_induction_coil_namelists

end module induction_coil_namelist
