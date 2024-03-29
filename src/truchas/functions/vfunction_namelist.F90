!!
!! VFUNCTION_NAMELIST
!!
!! This module provides procedures for reading the VFUNCTION namelists and for
!! accessing the functions they define.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_VFUNCTION_NAMELISTS (LUN) reads all the so-named namelists and
!!    creates the specified functions as VECTOR_FUNC class objects, which are
!!    held as private module data.
!!
!!  CALL LOOKUP_FUNC (NAME, F) returns the named function.  F is an allocatable
!!    class VECTOR_FUNC variable.  It is allocated by the subroutine and returns
!!    a copy (via sourced allocation) of the named function.  NAME is the value
!!    assigned to the NAME variable in the VFUNCTION namelist.  F is returned
!!    unallocated if no function with the given name exists.
!!

#include "f90_assert.fpp"

module vfunction_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use string_utilities, only: raise_case, i_to_c
  use input_utilities, only: seek_to_namelist
  use vector_func_factories
  use vector_func_table, only: known_func, insert_func
  use toolhead_driver, only: table
  use truchas_logging_services
  implicit none
  private

  public :: read_vfunction_namelists

  integer, parameter :: MAX_VAR = 10
  integer, parameter :: MAX_COEF = 64
  integer, parameter :: MAX_PARAM = 16

  integer,   parameter :: NULL_I = HUGE(1)
  real(r8),  parameter :: NULL_R = HUGE(1.0_r8)
  character, parameter :: NULL_C = char(0)

contains

  subroutine read_vfunction_namelists(lun)

    integer, intent(in) :: lun

    !! local variables
    logical :: found
    character(:), allocatable :: label
    integer :: n, ios, npts, rank, npar
    class(vector_func), allocatable :: f
    character(80) :: iom

    !! namelist variables
    character(len=31)  :: name, type, toolhead
    character(len=128) :: library_path, library_symbol
    real(r8) :: parameters(MAX_PARAM)
    real(r8) :: tabular_data(10,100), axis(3)
    integer  :: tabular_dim, dim

    namelist /vfunction/ name, type, tabular_data, tabular_dim, axis, toolhead, &
        library_path, library_symbol, dim, parameters

    call TLS_info('Reading VFUNCTION namelists ...')

    if (is_IOP) rewind(lun)
    n = 0

    do  ! until all VFUNCTION namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'VFUNCTION', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))

      call broadcast(found)
      if (.not.found) exit  ! no further VFUNCTION namelists found

      n = n + 1
      label = 'VFUNCTION[' // i_to_c(n) // ']'

      !! Default values
      name = NULL_C
      type = NULL_C
      tabular_data = NULL_R
      tabular_dim = NULL_I
      axis = NULL_R
      toolhead = NULL_C
      library_path = NULL_C
      library_symbol = NULL_C
      dim = NULL_I
      parameters = NULL_R

      !! Read the VFUNCTION namelist
      if (is_IOP) read(lun,nml=vfunction,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      !! Broadcast the namelist variables
      call broadcast(name)
      call broadcast(type)
      call broadcast(tabular_data)
      call broadcast(tabular_dim)
      call broadcast(axis)
      call broadcast(toolhead)
      call broadcast(library_path)
      call broadcast(library_symbol)
      call broadcast(dim)
      call broadcast(parameters)

      !! Check the user-supplied name.
      if (name == NULL_C .or. name == '') call TLS_fatal(label // ': NAME must be assigned a nonempty value')
      if (known_func(name)) then
        call TLS_fatal(label // ': another VFUNCTION namelist has this NAME: ' // trim(name))
      end if

      !! Identify the number of parameters.
      do npar = size(parameters), 1, -1
        if (parameters(npar) /= NULL_R) exit
      end do

      !! Check the type
      select case (raise_case(type))
      case ('TABULAR')
      case ('DIV-RADIAL-CYL-FLOW')
      case ('TOOLHEAD-LASER')
      case ('LIBRARY')
#ifndef ENABLE_DYNAMIC_LOADING
        call TLS_fatal (label // ': the configuration of this executable does not support TYPE="LIBRARY"')
#endif
      case default
        call TLS_fatal(label // ': unknown value for TYPE: ' // trim(type))
      end select

#ifdef ENABLE_DYNAMIC_LOADING
      if (raise_case(type) /= 'LIBRARY') then
        if (library_path /= NULL_C) &
            call TLS_warn (label // ': LIBRARY_PATH is ignored when TYPE="' // trim(type) // '"')
        if (library_symbol /= NULL_C) &
            call TLS_warn (label // ': LIBRARY_SYMBOL is ignored when TYPE="' // trim(type) // '"')
        if (dim /= NULL_I) &
            call TLS_warn (label // ': DIM is ignored when TYPE="' // trim(type) // '"')
        if (npar /= 0) &
            call TLS_warn (label // ': PARAMETERS is not used when TYPE="' // trim(type) // '"')
      end if
#endif

      !! Create the specified vfunction and add it to the vfunction table.

      select case (raise_case(type))
      case ('TABULAR')

        !! Identify and check the user-specified table.
        do npts = size(tabular_data,dim=2), 1, -1
          if (any(tabular_data(:,npts) /= NULL_R)) exit
        end do
        if (npts == 0) call TLS_fatal(label // ': TABULAR_DATA not specified')
        if (npts < 2) call TLS_fatal(label // ': TABULAR_DATA requires at least two data points')
        do rank = size(tabular_data,dim=1), 1, -1
          if (any(tabular_data(rank,:) /= NULL_R)) exit
        end do
        if (rank < 2) call TLS_fatal(label // ': TABULAR_DATA is incomplete')
        if (any(tabular_data(:rank,:npts) == NULL_R)) call TLS_fatal(label // ': TABULAR_DATA is incomplete')
        associate (xleft => tabular_data(1,1:npts-1), xright => tabular_data(1,2:npts))
          if (any(xleft >= xright)) call TLS_fatal(label // ': TABULAR_DATA is not ordered')
        end associate

        if (tabular_dim == NULL_I) then
          tabular_dim = 1
        else if (tabular_dim <= 0) then
          call TLS_fatal('TABULAR_DIM must be >0')
        end if

        call alloc_tabular_vector_func(f, tabular_data(1,:npts), tabular_data(2:rank,:npts), tabular_dim)
        call insert_func(name, f)

      case ('DIV-RADIAL-CYL-FLOW')

        if (any(axis == NULL_R)) call TLS_fatal(label // ': AXIS requires a 3-vector value')
        call alloc_div_radial_cyl_flow_func(f, axis)
        call insert_func(name, f)

      case ('TOOLHEAD-LASER')

        ! hackery to workaround conflicting input variable with type
        associate (thname => toolhead)
          block
            use toolhead_type
            type(toolhead), pointer :: th
            if (thname == NULL_C) call TLS_fatal(label // ': TOOLHEAD namelist name not specified')
            th => table%toolhead_ptr(thname)
            if (associated(th)) then
              call th%alloc_laser_func(f)
              call insert_func(name, f)
            else
              call TLS_fatal(label // ': unknown TOOLHEAD namelist: ' // trim(thname))
            end if
          end block
        end associate

#ifdef ENABLE_DYNAMIC_LOADING
      case ('LIBRARY')

        call alloc_dl_vector_func(f, library_path, library_symbol, dim, parameters(:npar))
        call insert_func(name, f)

#endif

      end select
      call TLS_info('  read namelist "' // trim(name) // '"')
    end do

    if (n == 0) call TLS_info('  none found')

  end subroutine read_vfunction_namelists

end module vfunction_namelist
