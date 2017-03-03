!!
!! FUNCTION_NAMELIST
!!
!! This module provides procedures for reading the FUNCTION namelists and for
!! accessing the functions they define.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_FUNCTION_NAMELISTS (LUN) reads all the so-named namelists and
!!    creates the specified functions as SCALAR_FUNC class objects, which are
!!    held as private module data.
!!
!!  CALL LOOKUP_FUNC (NAME, F) returns the named function.  F is an allocatable
!!    class SCALAR_FUNC variable.  It is allocated by the subroutine and returns
!!    a copy (via sourced allocation) of the named function.  NAME is the value
!!    assigned to the NAME variable in the FUNCTION namelist.  F is returned
!!    unallocated if no function with the given name exists.
!!

#include "f90_assert.fpp"

module function_namelist

  use kinds
  use parallel_communication
  use string_utilities, only: lower_case, raise_case, i_to_c
  use input_utilities, only: seek_to_namelist
  use scalar_func_factories
  use scalar_func_table, only: known_func, insert_func
  use ded_head_driver, only: alloc_ded_head_laser_func
  use truchas_logging_services
  implicit none
  private

  public :: read_function_namelists

  integer, parameter :: MAX_VAR = 10
  integer, parameter :: MAX_COEF = 64
  integer, parameter :: MAX_PARAM = 16

  integer,   parameter :: NULL_I = HUGE(1)
  real(r8),  parameter :: NULL_R = HUGE(1.0_r8)
  character, parameter :: NULL_C = char(0)

contains

  subroutine read_function_namelists (lun)

    integer, intent(in) :: lun

    !! local variables
    logical :: found, tabular_smooth
    integer :: n, stat, ncoef, nvar, npar, npts
    class(scalar_func), allocatable :: f

    !! namelist variables
    character(len=31)  :: name, type
    character(len=128) :: library_path, library_symbol
    character(len=8)   :: tabular_interp, tabular_extrap
    real(r8) :: poly_coefficients(MAX_COEF), poly_refvars(MAX_VAR), tabular_data(2,100)
    integer  :: poly_exponents(MAX_VAR,MAX_COEF), tabular_dim
    real(r8) :: parameters(MAX_PARAM)
    real(r8) :: smooth_step_x0, smooth_step_y0, smooth_step_x1, smooth_step_y1

    namelist/function/ name, type, parameters, library_path, library_symbol, &
        tabular_data, tabular_dim, tabular_interp, tabular_extrap, &
        poly_coefficients, poly_exponents, poly_refvars, &
        smooth_step_x0, smooth_step_y0, smooth_step_x1, smooth_step_y1

    call TLS_info ('')
    call TLS_info ('Reading FUNCTION namelists ...')

    if (is_IOP) rewind(lun)
    n = 0

    do  ! until all FUNCTION namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'FUNCTION', found, iostat=stat)
      call broadcast(stat)
      if (stat /= 0) call TLS_fatal ('error reading input file')

      call broadcast (found)
      if (.not.found) return  ! no further FUNCTION namelists found

      n = n + 1
      call TLS_info ('  Reading FUNCTION namelist #' // i_to_c(n))

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        name              = NULL_C
        type              = NULL_C
        poly_coefficients = NULL_R
        poly_exponents    = NULL_I
        poly_refvars      = NULL_R
        library_path      = NULL_C
        library_symbol    = NULL_C
        parameters        = NULL_R
        tabular_data      = NULL_R
        tabular_dim       = NULL_I
        tabular_interp    = NULL_C
        tabular_extrap    = NULL_C
        smooth_step_x0    = NULL_R
        smooth_step_y0    = NULL_R
        smooth_step_x1    = NULL_R
        smooth_step_y1    = NULL_R
        read(lun, nml=function, iostat=stat)
      end if

      call broadcast(stat)
      if (stat /= 0) call TLS_fatal ('error reading FUNCTION namelist')

      !! Broadcast the namelist variables.
      call broadcast (name)
      call broadcast (type)
      call broadcast (poly_coefficients)
      call broadcast (poly_exponents)
      call broadcast (poly_refvars)
      call broadcast (library_path)
      call broadcast (library_symbol)
      call broadcast (parameters)
      call broadcast (tabular_data)
      call broadcast (tabular_dim)
      call broadcast (tabular_interp)
      call broadcast (tabular_extrap)
      call broadcast (smooth_step_x0)
      call broadcast (smooth_step_y0)
      call broadcast (smooth_step_x1)
      call broadcast (smooth_step_y1)

      !! Check the user-supplied name.
      if (name == NULL_C .or. name == '') call TLS_fatal ('NAME must be assigned a nonempty value')
      if (known_func(name)) then
        call TLS_fatal ('already read a FUNCTION namelist with this name: ' // trim(name))
      end if

      !! Identify the number of parameters.
      do npar = size(parameters), 1, -1
        if (parameters(npar) /= NULL_R) exit
      end do

      !! Check the type.
      select case (raise_case(type))
      case ('POLYNOMIAL')
      case ('LIBRARY')
#ifndef ENABLE_DYNAMIC_LOADING
        call TLS_fatal ('the configuration of this executable does not support TYPE="LIBRARY"')
#endif
      case ('TABULAR')
      case ('SMOOTH STEP')
      case ('DED HEAD LASER')
      case (NULL_C)
        call TLS_fatal ('TYPE must be assigned a value')
      case default
        call TLS_fatal ('unknown value for TYPE: ' // trim(type))
      end select

      if (raise_case(type) /= 'POLYNOMIAL') then
        if (any(poly_coefficients /= NULL_R)) &
            call TLS_warn ('POLY_COEFFICIENTS is ignored when TYPE="' // trim(type) // '"')
        if (any(poly_exponents /= NULL_I)) &
            call TLS_warn ('POLY_EXPONENTS is ignored when TYPE="' // trim(type) // '"')
        if (any(poly_refvars /= NULL_R)) &
            call TLS_warn ('POLY_REFVARS is ignored when TYPE="' // trim(type) // '"')
      end if

#ifdef ENABLE_DYNAMIC_LOADING
      if (raise_case(type) /= 'LIBRARY') then
        if (library_path /= NULL_C) &
            call TLS_warn ('LIBRARY_PATH is ignored when TYPE="' // trim(type) // '"')
        if (library_symbol /= NULL_C) &
            call TLS_warn ('LIBRARY_SYMBOL is ignored when TYPE="' // trim(type) // '"')
        if (npar /= 0) &
            call TLS_warn ('PARAMETERS is not used when TYPE="' // trim(type) // '"')
      end if
#endif

      if (raise_case(type) /= 'TABULAR') then
        if (any(tabular_data /= NULL_R)) &
            call TLS_warn ('TABULAR_DATA is ignored when TYPE="' // trim(type) // '"')
        if (tabular_dim /= NULL_I) &
            call TLS_warn ('TABULAR_DIM is ignored when TYPE="' // trim(type) // '"')
        if (tabular_interp /= NULL_C) &
            call TLS_warn ('TABULAR_INTERP is ignored when TYPE="' // trim(type) // '"')
        if (tabular_extrap /= NULL_C) &
            call TLS_warn ('TABULAR_EXTRAP is ignored when TYPE="' // trim(type) // '"')
      end if

      if (raise_case(type) /= 'SMOOTH STEP') then
        if (smooth_step_x0 /= NULL_R) &
            call TLS_warn ('SMOOTH_STEP_X0 is ignored when TYPE="' // trim(type) // '"')
        if (smooth_step_y0 /= NULL_R) &
            call TLS_warn ('SMOOTH_STEP_Y0 is ignored when TYPE="' // trim(type) // '"')
        if (smooth_step_x1 /= NULL_R) &
            call TLS_warn ('SMOOTH_STEP_X1 is ignored when TYPE="' // trim(type) // '"')
        if (smooth_step_y1 /= NULL_R) &
            call TLS_warn ('SMOOTH_STEP_Y1 is ignored when TYPE="' // trim(type) // '"')
      end if

      !! Create the specified function and add it to the function table.

      select case (raise_case(type))
#ifdef ENABLE_DYNAMIC_LOADING
      case ('LIBRARY')

        call alloc_dl_scalar_func (f, library_path, library_symbol, parameters(:npar))
        call insert_func (name, f)

#endif
      case ('POLYNOMIAL')

        !! Identify the user-specified coefficients.  We require at least one,
        !! and they must fill the first part of the array -- no holes.
        do ncoef = size(poly_coefficients), 1, -1
          if (poly_coefficients(ncoef) /= NULL_R) exit
        end do
        if (ncoef == 0) then
          call TLS_fatal ('POLY_COEFFICIENTS must be assigned at least one value; none found')
        end if
        if (any(poly_coefficients(:ncoef) == NULL_R)) then
          call TLS_fatal ('values assigned to POLY_COEFFICIENTS are not packed')
        end if

        !! Identify the user-specified number of variables and exponents.
        !! We require a complete set of exponents for each coefficient
        do nvar = size(poly_exponents,dim=1), 1, -1
          if (any(poly_exponents(nvar,:ncoef) /= NULL_I)) exit
        end do
        if (nvar == 0) then
          call TLS_fatal ('POLY_EXPONENTS must be assigned values; none found')
        end if
        if (any(poly_exponents(:nvar,:ncoef) == NULL_I)) then
          call TLS_fatal ('not all required POLY_EXPONENTS values have been assigned')
        end if
        if (any(poly_exponents(:,ncoef+1:) /= NULL_I)) then
          call TLS_warn ('some values assigned to POLY_EXPONENTS are unused')
        end if

        !! Set the default for the reference variables.
        if (any(poly_refvars(nvar+1:) /= NULL_R)) then
          call TLS_warn ('some values assigned to POLY_REFVARS are unused')
        else
          where (poly_refvars == NULL_R) poly_refvars = 0.0_r8
        end if

        if (nvar == 1) then
          call alloc_poly_scalar_func (f, c = poly_coefficients(1:ncoef), &
                                       e = poly_exponents(1, 1:ncoef), &
                                       x0 = poly_refvars(1))
          call insert_func (name, f)
        else
          call alloc_mpoly_scalar_func (f, c = poly_coefficients(1:ncoef), &
                                        e = poly_exponents(1:nvar, 1:ncoef), &
                                        x0 = poly_refvars(1:nvar) )
          call insert_func (name, f)
        endif

      case ('TABULAR')

        !! Identify and check the user-specified table.
        do npts = size(tabular_data,dim=2), 1, -1
          if (any(tabular_data(:,npts) /= NULL_R)) exit
        end do
        if (npts == 0) then
          call TLS_fatal ('no values assigned to TABULAR_DATA')
        end if
        if (any(tabular_data(:,:npts) == NULL_R)) then
          call TLS_fatal ('values assigned to TABULAR_DATA are not packed')
        end if
        if (npts < 2) then
          call TLS_fatal ('at least two data points are required for TABULAR_DATA')
        end if
        associate (xleft => tabular_data(1,1:npts-1), xright => tabular_data(1,2:npts))
          if (any(xleft >= xright)) call TLS_fatal ('TABULAR_DATA is not ordered')
        end associate

        if (tabular_dim == NULL_I) then
          tabular_dim = 1
        else if (tabular_dim <= 0) then
          call TLS_fatal ('TABULAR_DIM must be >0')
        end if

        if (tabular_interp == NULL_C) tabular_interp = 'linear'
        select case (raise_case(tabular_interp))
        case ('LINEAR')
          tabular_smooth = .false.
        case ('AKIMA')
          tabular_smooth = .true.
        case default
          call TLS_fatal ('unknown value for TABULAR_INTERP: ' // trim(tabular_interp))
        end select

        if (tabular_extrap == NULL_C) tabular_extrap = 'nearest'
        tabular_extrap = lower_case(tabular_extrap)
        select case (tabular_extrap)
        case ('nearest')
        case ('linear')
        case default
          call TLS_fatal ('unknown value for TABULAR_EXTRAP: ' // trim(tabular_extrap))
        end select

        call alloc_tabular_scalar_func (f, tabular_data(1,:npts), tabular_data(2,:npts), &
            tabular_dim, tabular_smooth, tabular_extrap)
        call insert_func (name, f)

      case ('SMOOTH STEP')

        if (smooth_step_x0 == NULL_R) call TLS_fatal ('SMOOTH_STEP_X0 must be assigned a value')
        if (smooth_step_y0 == NULL_R) call TLS_fatal ('SMOOTH_STEP_Y0 must be assigned a value')
        if (smooth_step_x1 == NULL_R) call TLS_fatal ('SMOOTH_STEP_X1 must be assigned a value')
        if (smooth_step_y1 == NULL_R) call TLS_fatal ('SMOOTH_STEP_Y1 must be assigned a value')
        if (smooth_step_x0 >= smooth_step_x1) call TLS_fatal ('require SMOOTH_STEP_X0 < SMOOTH_STEP_X1')

        call alloc_sstep_scalar_func (f, smooth_step_x0, smooth_step_y0, smooth_step_x1, smooth_step_y1)
        call insert_func (name, f)

      case ('DED HEAD LASER')

        call alloc_ded_head_laser_func (f)
        call insert_func (name, f)

      end select
    end do

  end subroutine read_function_namelists

end module function_namelist
