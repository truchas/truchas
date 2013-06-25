
#include "f90_assert.fpp"

module function_namelist

  use kinds
  use parallel_communication
  use string_utilities, only: raise_case, i_to_c
  use input_utilities, only: seek_to_namelist
  use scalar_functions
  use function_table
  use ds_utilities

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
    logical :: found
    integer :: n, stat, ncoef, nvar, npar
    type(scafun) :: f

    !! namelist variables
    character(len=FT_MAX_NAME_LEN) :: name, type
    character(len=128) :: library_path, library_symbol
    real(r8) :: poly_coefficients(MAX_COEF), poly_refvars(MAX_VAR)
    integer  :: poly_exponents(MAX_VAR,MAX_COEF)
    real(r8) :: parameters(MAX_PARAM)
    real(r8) :: smooth_step_x0, smooth_step_y0, smooth_step_x1, smooth_step_y1

    namelist/function/ name, type, parameters, library_path, library_symbol, &
        poly_coefficients, poly_exponents, poly_refvars, &
        smooth_step_x0, smooth_step_y0, smooth_step_x1, smooth_step_y1

    call ds_info ('')
    call ds_info ('Reading FUNCTION namelists ...')

    if (is_IOP) rewind(lun)
    n = 0

    do  ! until all FUNCTION namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'FUNCTION', found, iostat=stat)
      call broadcast(stat)
      if (stat /= 0) call ds_halt ('error reading input file')

      call broadcast (found)
      if (.not.found) return  ! no further FUNCTION namelists found

      n = n + 1
      call ds_info ('  Reading FUNCTION namelist #' // i_to_c(n))

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
        smooth_step_x0    = NULL_R
        smooth_step_y0    = NULL_R
        smooth_step_x1    = NULL_R
        smooth_step_y1    = NULL_R
        read(lun, nml=function, iostat=stat)
      end if

      call broadcast(stat)
      if (stat /= 0) call ds_halt ('error reading FUNCTION namelist')

      !! Broadcast the namelist variables.
      call broadcast (name)
      call broadcast (type)
      call broadcast (poly_coefficients)
      call broadcast (poly_exponents)
      call broadcast (poly_refvars)
      call broadcast (library_path)
      call broadcast (library_symbol)
      call broadcast (parameters)
      call broadcast (smooth_step_x0)
      call broadcast (smooth_step_y0)
      call broadcast (smooth_step_x1)
      call broadcast (smooth_step_y1)

      !! Check the user-supplied name.
      if (name == NULL_C .or. name == '') call ds_halt ('NAME must be assigned a nonempty value')
      if (ft_has_function(name)) then
        call ds_halt ('already read a FUNCTION namelist with this name: ' // trim(name))
      end if
      
      !! Identify the number of parameters.
      do npar = size(parameters), 1, -1
        if (parameters(npar) /= NULL_R) exit
      end do

      !! Check the type.
      select case (raise_case(type))
      case ('POLYNOMIAL')
      case ('LIBRARY')
      case ('SMOOTH STEP')
      case (NULL_C)
        call ds_halt('TYPE must be assigned a value')
      case default
        call ds_halt('unknown value for TYPE: ' // trim(type))
      end select

      if (raise_case(type) /= 'POLYNOMIAL') then
        if (any(poly_coefficients /= NULL_R)) &
            call ds_warn ('POLY_COEFFICIENTS is ignored when TYPE="' // trim(type) // '"')
        if (any(poly_exponents /= NULL_I)) &
            call ds_warn ('POLY_EXPONENTS is ignored when TYPE="' // trim(type) // '"')
        if (any(poly_refvars /= NULL_R)) &
            call ds_warn ('POLY_REFVARS is ignored when TYPE="' // trim(type) // '"')
      end if

      if (raise_case(type) /= 'LIBRARY') then
        if (library_path /= NULL_C) &
            call ds_warn ('LIBRARY_PATH is ignored when TYPE="' // trim(type) // '"')
        if (library_symbol /= NULL_C) &
            call ds_warn ('LIBRARY_SYMBOL is ignored when TYPE="' // trim(type) // '"')
        if (npar /= 0) &
            call ds_warn ('PARAMETERS is not used when TYPE="' // trim(type) // '"')
      end if

      if (raise_case(type) /= 'SMOOTH STEP') then
        if (smooth_step_x0 /= NULL_R) &
            call ds_warn ('SMOOTH_STEP_X0 is ignored when TYPE="' // trim(type) // '"')
        if (smooth_step_y0 /= NULL_R) &
            call ds_warn ('SMOOTH_STEP_Y0 is ignored when TYPE="' // trim(type) // '"')
        if (smooth_step_x1 /= NULL_R) &
            call ds_warn ('SMOOTH_STEP_X1 is ignored when TYPE="' // trim(type) // '"')
        if (smooth_step_y1 /= NULL_R) &
            call ds_warn ('SMOOTH_STEP_Y1 is ignored when TYPE="' // trim(type) // '"')
      end if

      !! Create the specified function and add it to the function table.

      select case (raise_case(type))
      case ('LIBRARY')

        call create_scafun_dll (f, lib=library_path, sym=library_symbol, p=parameters(:npar))
        call ft_add_function (name, f)
        call destroy (f)

      case ('POLYNOMIAL')

        !! Identify the user-specified coefficients.  We require at least one,
        !! and they must fill the first part of the array -- no holes.
        do ncoef = size(poly_coefficients), 1, -1
          if (poly_coefficients(ncoef) /= NULL_R) exit
        end do
        if (ncoef == 0) then
          call ds_halt ('POLY_COEFFICIENTS must be assigned at least one value; none found')
        end if
        if (any(poly_coefficients(:ncoef) == NULL_R)) then
          call ds_halt ('values assigned to POLY_COEFFICIENTS are not packed')
        end if

        !! Identify the user-specified number of variables and exponents.
        !! We require a complete set of exponents for each coefficient
        do nvar = size(poly_exponents,dim=1), 1, -1
          if (any(poly_exponents(nvar,:ncoef) /= NULL_I)) exit
        end do
        if (nvar == 0) then
          call ds_halt ('POLY_EXPONENTS must be assigned values; none found')
        end if
        if (any(poly_exponents(:nvar,:ncoef) == NULL_I)) then
          call ds_halt ('not all required POLY_EXPONENTS values have been assigned')
        end if
        if (any(poly_exponents(:,ncoef+1:) /= NULL_I)) then
          call ds_warn ('some values assigned to POLY_EXPONENTS are unused')
        end if

        !! Set the default for the reference variables.
        if (any(poly_refvars(nvar+1:) /= NULL_R)) then
          call ds_warn ('some values assigned to POLY_REFVARS are unused')
        else
          where (poly_refvars == NULL_R) poly_refvars = 0.0_r8
        end if

        if (nvar == 1) then
          call create_scafun_poly (f, c = poly_coefficients(1:ncoef), &
                                      e = poly_exponents(1, 1:ncoef), &
                                      x0 = poly_refvars(1))
          call ft_add_function (name, f)
          call destroy (f)
        else
          call create_scafun_mpoly (f, c = poly_coefficients(1:ncoef), &
                                       e = poly_exponents(1:nvar, 1:ncoef), &
                                       x0 = poly_refvars(1:nvar) )
          call ft_add_function (name, f)
          call destroy (f)
        endif

      case ('SMOOTH STEP')

        if (smooth_step_x0 == NULL_R) call ds_halt ('SMOOTH_STEP_X0 must be assigned a value')
        if (smooth_step_y0 == NULL_R) call ds_halt ('SMOOTH_STEP_Y0 must be assigned a value')
        if (smooth_step_x1 == NULL_R) call ds_halt ('SMOOTH_STEP_X1 must be assigned a value')
        if (smooth_step_y1 == NULL_R) call ds_halt ('SMOOTH_STEP_Y1 must be assigned a value')
        if (smooth_step_x0 >= smooth_step_x1) call ds_halt ('require SMOOTH_STEP_X0 < SMOOTH_STEP_X1')
        call create_scafun_sstep (f, smooth_step_x0, smooth_step_y0, smooth_step_x1, smooth_step_y1)
        call ft_add_function (name, f)
        call destroy (f)

      end select
    end do

  end subroutine read_function_namelists

end module function_namelist
