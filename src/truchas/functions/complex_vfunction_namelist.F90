!!
!! COMPLEX_VFUNCTION_NAMELIST
!!
!! This module provides procedures for reading the COMPLEX_VFUNCTION namelists.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_COMPLEX_VFUNCTION_NAMELISTS(LUN) reads all the so-named namelists
!!    and creates the specified functions as COMPLEX_VECTOR_FUNC class objects,
!!    and adds them to the truchas function table.
!!

#include "f90_assert.fpp"

module complex_vfunction_namelist

  implicit none
  private

  public :: read_complex_vfunction_namelists

  integer, parameter :: MAX_PARAM = 16

contains

  subroutine read_complex_vfunction_namelists(lun)
  
    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c, lower_case
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use complex_vector_func_factories
    use func_table, only: known_func, insert_func
    use truchas_logging_services

    integer, intent(in) :: lun

    logical :: found
    character(80) :: iom
    character(:), allocatable :: label
    integer :: n, ios, npar
    class(complex_vector_func), allocatable :: f

    !! Namelist variables
    character(31)  :: name, type
    character(128) :: library_path, library_symbol
    real(r8) :: parameters(MAX_PARAM)
    integer  :: dim

    namelist /complex_vfunction/ name, type, library_path, library_symbol, dim, parameters

    call TLS_info('Reading COMPLEX_VFUNCTION namelists ...')

    if (is_IOP) rewind(lun)
    n = 0

    do  ! until all COMPLEX_VFUNCTION namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'complex_vfunction', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit ! no further COMPLEX_VFUNCTION namelists found

      n = n + 1
      label = 'COMPLEX_VFUNCTION[' // i_to_c(n) // ']'

      !! Default values
      name = NULL_C
      type = NULL_C
      library_path = NULL_C
      library_symbol = NULL_C
      dim = NULL_I
      parameters = NULL_R

      !! Read the namelist
      if (is_IOP) read(lun,nml=complex_vfunction,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(type)
      call broadcast(library_path)
      call broadcast(library_symbol)
      call broadcast(dim)
      call broadcast(parameters)

      !! Check the user-supplied name
      if (name == NULL_C .or. name == '') call TLS_fatal(label // ': NAME must be assigned a nonempty value')
      if (known_func(name)) then
        call TLS_fatal(label // ': another function namelist has this NAME: ' // trim(name))
      end if

      !! Check the type
      select case (lower_case(type))
      case ('library')
#ifndef ENABLE_DYNAMIC_LOADING
        call TLS_fatal(label // ': the configuration of this executable does not support TYPE="LIBRARY"')
#endif
      case (NULL_C)
        call TLS_fatal(label // ': TYPE must be assigned a value')
      case default
        call TLS_fatal(label // ': unknown value for TYPE: ' // trim(type))
      end select

      !! Create the specified function and add it to the function table
      select case (lower_case(type))
#ifdef ENABLE_DYNAMIC_LOADING
      case ('library')
        !! Identify the number of parameters
        do npar = size(parameters), 1, -1
          if (parameters(npar) /= NULL_R) exit
        end do
        if (dim == NULL_I) call TLS_fatal(label // ': DIM not specified')
        if (dim <= 0) call TLS_fatal(label // ': DIM must be >= 1')
        call alloc_dl_complex_vector_func(f, library_path, library_symbol, dim, parameters(:npar))
        call insert_func(name, f)
#endif
      end select
      call TLS_info('  read namelist "' // trim(name) // '"')
    end do

    if (n == 0) call TLS_info('  none found')

  end subroutine read_complex_vfunction_namelists

end module complex_vfunction_namelist
