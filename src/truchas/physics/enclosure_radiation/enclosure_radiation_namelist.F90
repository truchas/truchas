!!
!! ENCLOSURE_RADIATION_NAMELIST
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Refactored, June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module enclosure_radiation_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use parameter_list_type
  use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
  use string_utilities, only: i_to_c, raise_case
  use truchas_env, only: input_dir
  use truchas_logging_services
  implicit none
  private

  public :: read_enclosure_radiation_namelists

  type(parameter_list), public :: params

  integer, parameter :: MAX_NAME_LEN = 31, MAX_FILE_LEN = 255, MAX_FACE_BLOCK_IDS = 32

contains

  subroutine read_enclosure_radiation_namelists(lun)

    use toolpath_table, only: known_toolpath

    integer, intent(in) :: lun

    integer :: n, ios
    logical :: found
    character(:), allocatable :: label
    character(128) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer  :: precon_iter
    logical  :: skip_geometry_check
    real(r8) :: coord_scale_factor, ambient_constant, error_tolerance
    character(MAX_FILE_LEN) :: enclosure_file
    character(MAX_NAME_LEN) :: name, ambient_function, precon_method, precon_coupling_method, toolpath
    namelist /enclosure_radiation/ name, enclosure_file, coord_scale_factor, skip_geometry_check, &
                                   ambient_constant, ambient_function, error_tolerance, &
                                   precon_method, precon_iter, precon_coupling_method, &
                                   toolpath

    call TLS_info('')
    call TLS_info('Reading ENCLOSURE_RADIATION namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all ENCLOSURE_RADIATION namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'enclosure_radiation', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'ENCLOSURE_RADIATION[' // i_to_c(n) // ']'

      !! Replicate the namelist variables on all processes before reading
      name = NULL_C
      toolpath = NULL_C
      enclosure_file = NULL_C
      coord_scale_factor = NULL_R
      skip_geometry_check = .false.
      ambient_constant = NULL_R
      ambient_function = NULL_C
      error_tolerance = NULL_R
      precon_method = NULL_C
      precon_iter = NULL_I
      precon_coupling_method = NULL_C

      if (is_IOP) read(lun,nml=enclosure_radiation,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      !! Replicate the namelist variables on all processes
      call broadcast(name)
      call broadcast(toolpath)
      call broadcast(enclosure_file)
      call broadcast(coord_scale_factor)
      call broadcast(skip_geometry_check)
      call broadcast(ambient_constant)
      call broadcast(ambient_function)
      call broadcast(error_tolerance)
      call broadcast(precon_method)
      call broadcast(precon_iter)
      call broadcast(precon_coupling_method)

      !! A unique NAME is required; becomes the enclosure sublist parameter name
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (params%is_parameter(trim(name))) then
        call TLS_fatal(label // ': another ENCLOSURE_RADIATION has this NAME: ' // trim(name))
      else
        plist => params%sublist(trim(name))
      end if

      !! Check optional TOOLPATH
      if (toolpath /= NULL_C) then
        if (.not.known_toolpath(toolpath)) &
            call TLS_fatal(label // ': unknown TOOLPATH: ' // trim(toolpath))
        call plist%set('toolpath', toolpath)
      end if

      !! Verify ENCLOSURE_FILE was assigned a value
      if (enclosure_file == NULL_C) then
        call TLS_fatal(label // ': ENCLOSURE_FILE not specified')
      end if

      !! Fix-up the ENCLOSURE_FILE path
      if (enclosure_file(1:1) /= '/') then  ! not an absolute path
        enclosure_file = trim(input_dir) // trim(enclosure_file)
      end if

      !! Check that ENCLOSURE_FILE can at least be found (when TOOLPATH not specified)
      if (.not.plist%is_parameter('toolpath')) then
        if (is_IOP) inquire(file=enclosure_file,exist=found)
        call broadcast(found)
        if (.not.found) call TLS_fatal(label // ': no such ENCLOSURE_FILE: ' // trim(enclosure_file))
      end if
      call plist%set('enclosure-filename', trim(enclosure_file))

      !! Check COORD_SCALE_FACTOR (optional)
      if (coord_scale_factor /= NULL_R) then
        if (coord_scale_factor <= 0.0_r8) call TLS_fatal(label // ': COORD_SCALE_FACTOR must be > 0')
        call plist%set('coord-scale-factor', coord_scale_factor)
      end if

      !! Check AMBIENT_CONSTANT/AMBIENT_FUNCTION
      if (ambient_constant /= NULL_R .and. ambient_function /= NULL_C) then
        call TLS_fatal(label // ': both AMBIENT_CONSTANT and AMBIENT_FUNCTION specified')
      else if (ambient_constant /= NULL_R) then
        call plist%set('ambient-temp', ambient_constant)
      else if (ambient_function /= NULL_C) then
        call plist%set('ambient-temp', trim(ambient_function))
      else
        call TLS_fatal(label // ': neither AMBIENT_CONSTANT nor AMBIENT_FUNCTION specified')
      end if

      !! Check ERROR_TOLERANCE
      if (error_tolerance == NULL_R) then
        error_tolerance = 1.0e-3_r8
      else if (error_tolerance <= 0.0_r8) then
        call TLS_fatal(label // ': ERROR_TOLERANCE must be > 0')
      else
        call plist%set('error-tol', error_tolerance)
      end if

      !! Check PRECON_METHOD (optional, default 'JACOBI')
      if (precon_method /= NULL_C) then
        precon_method = raise_case(precon_method)
        select case (precon_method)
        case ('JACOBI')
        case ('CHEBYSHEV')
        case default
          call TLS_fatal(label // ': unknown PRECON_METHOD: ' // trim(precon_method))
        end select
        call plist%set('precon-method', precon_method)
      end if

      !! Check PRECON_ITER (optional, default 1)
      if (precon_iter /= NULL_I) then
        if (precon_iter < 1) call TLS_fatal(label // ': PRECON_ITER must be > 0')
        call plist%set('precon-iter', precon_iter)
      end if

      !! Check PRECON_COUPLING_METHOD (optional, default 'BACKWARD GS')
      if (precon_coupling_method /= NULL_C) then
        precon_coupling_method = raise_case(precon_coupling_method)
        select case (precon_coupling_method)
        case ('JACOBI')
        case ('FORWARD GS')
        case ('BACKWARD GS')
        case ('FACTORIZATION')
        case default
          call TLS_fatal(label // ': unknown PRECON_COUPLING_METHOD: ' // trim(precon_coupling_method))
        end select
        call plist%set('precon-coupling-method', precon_coupling_method)
      end if

    end do

    call read_enclosure_surface_namelists(lun)

  end subroutine read_enclosure_radiation_namelists


  subroutine read_enclosure_surface_namelists(lun)

    integer, intent(in) :: lun

    integer :: n, ios
    logical :: found
    character(:), allocatable :: label
    character(128) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer :: face_block_ids(MAX_FACE_BLOCK_IDS)
    real(r8) :: emissivity_constant
    character(MAX_NAME_LEN) :: name, enclosure_name, emissivity_function
    namelist /enclosure_surface/ name, enclosure_name, face_block_ids, &
                                 emissivity_constant, emissivity_function

    call TLS_info('')
    call TLS_info('Reading ENCLOSURE_SURFACE namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all ENCLOSURE_SURFACE namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'enclosure_surface', found, iostat=ios)
      call broadcast (ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'ENCLOSURE_SURFACE[' // i_to_c(n) // ']'

      !! Assign default values to namelist variables before reading
      name = NULL_C
      enclosure_name = NULL_C
      face_block_ids = NULL_I
      emissivity_constant = NULL_R
      emissivity_function = NULL_C

      if (is_IOP) read(lun,nml=enclosure_surface,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      !! Replicate the namelist variables on all processes
      call broadcast(name)
      call broadcast(enclosure_name)
      call broadcast(face_block_ids)
      call broadcast(emissivity_constant)
      call broadcast(emissivity_function)

      !! Check ENCLOSURE_NAME (required)
      if (enclosure_name == NULL_C) then
        call TLS_fatal(label // ': ENCLOSURE_NAME not specified')
      else if (params%is_sublist(trim(enclosure_name))) then
        plist => params%sublist(trim(enclosure_name))
        plist => plist%sublist('surfaces')
      else
        call TLS_fatal(label // ': unknown ENCLOSURE_NAME: ' // trim(enclosure_name))
      end if

      !! A unique NAME is required; becomes the surface sublist name
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (plist%is_sublist(trim(name))) then
        call TLS_fatal(label // ': another ENCLOSURE_SURFACE has this NAME: ' // trim(name))
      else
        plist => plist%sublist(trim(name))
      end if

      !! Check for a non-empty FACE_BLOCK_IDS.
      if (count(face_block_ids /= NULL_I) == 0) then
        call TLS_fatal(label // ': FACE_BLOCK_IDS not specified')
      else
       call plist%set('face-block-ids', pack(face_block_ids, mask=(face_block_ids /= NULL_I)))
      end if

      !! Verify that only one of EMISSIVITY_CONSTANT and EMISSIVITY_FUNCTION was specified.
      if (emissivity_constant /= NULL_R .and. emissivity_function /= NULL_C) then
        call TLS_fatal(label // ': both EMISSIVITY_CONSTANT and EMISSIVITY_FUNCTION specified')
      else if (emissivity_constant /= NULL_R) then
        if (emissivity_constant < 0.0_r8 .or. emissivity_constant > 1.0_r8) then
          call TLS_fatal(label // ': EMISSIVITY_CONSTANT must be in [0,1]')
        else if (emissivity_constant == 0.0_r8) then
          call TLS_warn(label // ': EMISSIVITY_CONSTANT is 0')
        end if
        call plist%set('emissivity', emissivity_constant)
      else if (emissivity_function /= NULL_C) then
        call plist%set('emissivity', trim(emissivity_function))
      else
        call TLS_fatal(label // ': neither EMISSIVITY_CONSTANT nor EMISSIVITY_FUNCTION specified')
      end if

      call TLS_info (trim(label) // ' read ENCLOSURE_SURFACE namelist "' // trim(name) // '"')

    end do

  end subroutine read_enclosure_surface_namelists

end module enclosure_radiation_namelist
