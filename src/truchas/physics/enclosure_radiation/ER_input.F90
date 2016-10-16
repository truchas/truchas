!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ER_input

  use kinds
  use parallel_communication
  use function_namelist, only: lookup_func
  use scalar_func_class
  use truchas_logging_services
  use string_utilities, only: i_to_c, raise_case
  implicit none
  private

  public :: ERI_read_enclosure_radiation
  public :: ERI_num_enclosures, ERI_get_names
  public :: ERI_get_file, ERI_get_coord_scale_factor, ERI_get_ambient
  public :: ERI_get_error_tolerance, ERI_get_precon_method, ERI_get_precon_coupling_method
  public :: ERI_check_geometry

  public :: ERI_read_enclosure_surface
  public :: ERI_get_emissivity

  integer, parameter :: MAX_NAME_LEN = 31, MAX_FILE_LEN = 255, MAX_FACE_BLOCK_IDS = 32

  integer,   parameter :: NULL_I = HUGE(1)
  real(r8),  parameter :: NULL_R = HUGE(1.0_r8)
  character, parameter :: NULL_C = char(0)

  !! Database for the enclosure_radiation namelists
  type :: er_list_node
    integer :: seq
    character(len=MAX_NAME_LEN) :: name
    character(len=MAX_FILE_LEN) :: file
    real(r8) :: csf
    class(scalar_func), allocatable :: tamb
    real(r8) :: tol
    character(len=MAX_NAME_LEN) :: pm, pcm
    integer :: numitr
    logical :: geom = .true.
    type(er_list_node), pointer :: next => null()
  end type er_list_node
  type(er_list_node), pointer, save :: er_list => null(), er_last => null()

  !! Database for the enclosure_surface namelists
  type :: es_list_node
    integer :: seq
    character(len=MAX_NAME_LEN) :: name, encl_name
    integer, pointer :: face_block_ids(:) => null()
    class(scalar_func), allocatable :: eps
    type(es_list_node), pointer :: next => null()
  end type es_list_node
  type(es_list_node), pointer, save :: es_list => null(), es_last => null()

contains

  subroutine ERI_read_enclosure_radiation (lun)

    use input_utilities, only: seek_to_namelist
    use truchas_env, only: input_dir
    use scalar_func_factories, only: alloc_const_scalar_func

    integer, intent(in) :: lun

    logical :: found
    integer :: n, stat
    class(scalar_func), allocatable :: f
    character(len=7) :: label
    character(len=127) :: errmsg

    !! Namelist variables
    integer  :: precon_iter
    logical  :: skip_geometry_check
    real(r8) :: coord_scale_factor, ambient_constant, error_tolerance
    character(len=MAX_FILE_LEN) :: enclosure_file
    character(len=MAX_NAME_LEN) :: name, ambient_function, precon_method, precon_coupling_method
    namelist /enclosure_radiation/ name, enclosure_file, coord_scale_factor, skip_geometry_check, &
                                   ambient_constant, ambient_function, error_tolerance, &
                                   precon_method, precon_iter, precon_coupling_method

    call TLS_info ('')
    call TLS_info ('Reading ENCLOSURE_RADIATION namelists ...')

    if (is_IOP) rewind(lun)
    n = 0 ! namelist counter
    stat = 0

    do  ! until all RADIATION_ENCLOSURE namelists have been read or an error occurs.

      n = n + 1
      write(label,'(a,i0,a)') '  [', n, ']'

      if (is_IOP) call seek_to_namelist (lun, 'ENCLOSURE_RADIATION', found, iostat=stat)

      call broadcast (stat)
      if (stat /= 0) then
        errmsg = trim(label) // ' seek error, iostat=' // i_to_c(stat)
        exit
      end if

      call broadcast (found)
      if (.not.found) exit

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        name = NULL_C
        enclosure_file = NULL_C
        coord_scale_factor = NULL_R
        skip_geometry_check = .false.
        ambient_constant = NULL_R
        ambient_function = NULL_C
        error_tolerance = NULL_R
        precon_method = NULL_C
        precon_iter = NULL_I
        precon_coupling_method = NULL_C
        read(lun,nml=enclosure_radiation,iostat=stat)
      end if

      call broadcast (stat)
      if (stat /= 0) then
        errmsg = trim(label) // ' read error, iostat=' // i_to_c(stat)
        exit
      end if

      !! Replicate the namelist variables on all processes.
      call broadcast (name)
      call broadcast (enclosure_file)
      call broadcast (coord_scale_factor)
      call broadcast (skip_geometry_check)
      call broadcast (ambient_constant)
      call broadcast (ambient_function)
      call broadcast (error_tolerance)
      call broadcast (precon_method)
      call broadcast (precon_iter)
      call broadcast (precon_coupling_method)

      !! Check the user-supplied NAME.
      if (name == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: no value assigned to NAME'
        exit
      else if (name_exists(name)) then
        stat = -1
        errmsg = trim(label) // ' error: another ENCLOSURE_RADIATION has this NAME: ' // trim(name)
        exit
      end if

      call TLS_info (trim(label) // ' read ENCLOSURE_RADIATION namelist "' // trim(name) // '"')

      !! Verify ENCLOSURE_FILE was assigned a value.
      if (enclosure_file == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: no value assigned to ENCLOSURE_FILE'
        exit
      end if

      !! Fix-up the ENCLOSURE_FILE path.
      if (enclosure_file(1:1) /= '/') then  ! not an absolute path
        enclosure_file = trim(input_dir) // trim(enclosure_file)
      end if

      !! Check that ENCLOSURE_FILE can at least be found.
      if (is_IOP) inquire(file=enclosure_file,exist=found)
      call broadcast (found)
      if (.not.found) then
        stat = -1
        errmsg = trim(label) // ' error: no such ENCLOSURE_FILE: ' // trim(enclosure_file)
        exit
      end if

      !! Check COORD_SCALE_FACTOR and assign default if necessary.
      if (coord_scale_factor == NULL_R) then
        coord_scale_factor = 1.0_r8
      else if (coord_scale_factor <= 0.0_r8) then
        stat = -1
        errmsg = trim(label) // ' error: COORD_SCALE_FACTOR must be > 0'
        exit
      end if

      !! Verify that only one of AMBIENT_CONSTANT and AMBIENT_FUNCTION was specified.
      if (ambient_constant == NULL_R .eqv. ambient_function == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: exactly one of AMBIENT_CONSTANT and &
                                &AMBIENT_FUNCTION must be assigned a value'
        exit
      end if

      !! Get or create the specified ambient temperature function.
      if (ambient_constant /= NULL_R) then
        call alloc_const_scalar_func (f, ambient_constant)
      else
        call lookup_func (ambient_function, f)
        if (.not.allocated(f)) then
          stat = -1
          errmsg = trim(label) // ' error: unknown function name: ' // trim(ambient_function)
          exit
        end if
      end if
      
      !! Check ERROR_TOLERANCE.
      if (error_tolerance == NULL_R) then
        error_tolerance = 1.0e-3_r8
      else if (error_tolerance <= 0.0_r8) then
        stat = -1
        errmsg = trim(label) // ' error: ERROR_TOLERANCE must be > 0'
        exit
      end if

      !! Check PRECON_METHOD.
      if (precon_method == NULL_C) then
        precon_method = 'JACOBI'
        call TLS_info (trim(label)//' using default PRECON_METHOD="'//trim(precon_method)//'"')
      end if
      precon_method = raise_case(precon_method)
      select case (precon_method)
      case ('JACOBI')
      case ('CHEBYSHEV')
      case default
        stat = -1
        errmsg = trim(label) // ' error: unknown PRECON_METHOD value: "' // trim(precon_method) // '"'
        exit
      end select

      !! Check PRECON_ITER.
      if (precon_iter == NULL_I) then
        precon_iter = 1
        write(errmsg,'(2a,i0)' )trim(label), ' using default PRECON_ITER=', precon_iter
        call TLS_info (errmsg)
      else if (precon_iter < 1) then
        stat = -1
        errmsg = trim(label) // ' error: PRECON_ITER must be > 0'
        exit
      end if

      !! Check PRECON_COUPLING_METHOD.
      if (precon_coupling_method == NULL_C) then
        precon_coupling_method = 'BACKWARD GS'
        call TLS_info (trim(label)//' using default PRECON_COUPLING_METHOD="'//trim(precon_coupling_method)//'"')
      end if
      precon_coupling_method = raise_case(precon_coupling_method)
      select case (precon_coupling_method)
      case ('JACOBI')
      case ('FORWARD GS')
      case ('BACKWARD GS')
      case ('FACTORIZATION')
      case default
        stat = -1
        errmsg = trim(label) // ' error: unknown PRECON_COUPLING_METHOD value: "' // trim(precon_coupling_method) // '"'
        exit
      end select

      !! Append the data for this namelist to the list of namelist data.
      if (associated(er_last)) then
        allocate(er_last%next)
        er_last => er_last%next
      else  ! list is empty
        allocate(er_list)
        er_last => er_list
      end if
      er_last%seq  = n
      er_last%name = name
      er_last%file = enclosure_file
      er_last%csf  = coord_scale_factor
      er_last%geom = .not.skip_geometry_check
      call move_alloc (f, er_last%tamb)
      er_last%tol  = error_tolerance
      er_last%pm   = precon_method
      er_last%pcm  = precon_coupling_method
      er_last%numitr = precon_iter

    end do

    if (stat /= 0) then
      call TLS_info (trim(errmsg))
      call TLS_fatal ('error processing ENCLOSURE_RADIATION namelists')
    end if

  contains

    logical function name_exists (name)
      character(len=*), intent(in) :: name
      type(er_list_node), pointer :: l
      l => er_list
      name_exists = .true.
      do while (associated(l))
        if (l%name == name) return
        l => l%next
      end do
      name_exists = .false.
    end function name_exists

  end subroutine ERI_read_enclosure_radiation

  function ERI_num_enclosures () result (n)
    type(er_list_node), pointer :: nml
    integer :: n
    n = 0
    nml => er_list
    do while (associated(nml))
      n = n + 1
      nml => nml%next
    end do
  end function ERI_num_enclosures

  subroutine ERI_get_names (names)
    character(len=*), intent(out) :: names(:)
    type(er_list_node), pointer :: nml
    integer :: n
    n = 0
    nml => er_list
    do while (associated(nml))
      INSIST(size(names) > n)
      n = n + 1
      names(n) = nml%name
      nml => nml%next
    end do
  end subroutine ERI_get_names

  subroutine ERI_get_coord_scale_factor (name, csf)
    character(len=*), intent(in)  :: name
    real(r8), intent(out) :: csf
    type(er_list_node), pointer :: nml
    nml => er_namelist_with_name(name)
    INSIST(associated(nml))
    csf = nml%csf
  end subroutine ERI_get_coord_scale_factor
  
  logical function ERI_check_geometry (name)
    character(len=*), intent(in) :: name
    type(er_list_node), pointer :: nml
    nml => er_namelist_with_name(name)
    INSIST(associated(nml))
    ERI_check_geometry = nml%geom
  end function ERI_check_geometry

  subroutine ERI_get_file (name, file)
    character(len=*), intent(in)  :: name
    character(:), allocatable, intent(out) :: file
    type(er_list_node), pointer :: nml
    nml => er_namelist_with_name(name)
    INSIST(associated(nml))
    file = trim(nml%file)
  end subroutine ERI_get_file

  subroutine ERI_get_ambient (name, tamb)
    character(len=*), intent(in)  :: name
    class(scalar_func), allocatable, intent(out) :: tamb
    type(er_list_node), pointer :: nml
    nml => er_namelist_with_name(name)
    INSIST(associated(nml))
    call move_alloc (nml%tamb, tamb)
  end subroutine ERI_get_ambient

  subroutine ERI_get_error_tolerance (name, tol)
    character(len=*), intent(in)  :: name
    real(r8), intent(out) :: tol
    type(er_list_node), pointer :: nml
    nml => er_namelist_with_name(name)
    INSIST(associated(nml))
    tol = nml%tol
  end subroutine ERI_get_error_tolerance

  subroutine ERI_get_precon_method (name, method, numitr)
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: method
    integer, intent(out) :: numitr
    type(er_list_node), pointer :: nml
    nml => er_namelist_with_name(name)
    INSIST(associated(nml))
    INSIST(len(method) >= len_trim(nml%pm))
    method = nml%pm
    numitr = nml%numitr
  end subroutine ERI_get_precon_method

  subroutine ERI_get_precon_coupling_method (name, method)
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: method
    type(er_list_node), pointer :: nml
    nml => er_namelist_with_name(name)
    INSIST(associated(nml))
    INSIST(len(method) >= len_trim(nml%pcm))
    method = nml%pcm
  end subroutine ERI_get_precon_coupling_method

  function er_namelist_with_name (name) result (nml)
    character(len=*), intent(in) :: name
    type(er_list_node), pointer :: nml
    nml => er_list
    do while (associated(nml))
      if (name == nml%name) return
      nml => nml%next
    end do
  end function er_namelist_with_name

!&ENCLOSURE_SURFACE
!  name = string
!  enclosure_name = string
!  face_block_ids = list-of-integers
!  emissivity_constant = real
!  emissivity_function = string
!/

  subroutine ERI_read_enclosure_surface (lun)

    use input_utilities, only: seek_to_namelist
    use scalar_func_factories, only: alloc_const_scalar_func

    integer, intent(in) :: lun

    logical :: found
    integer :: n, stat
    class(scalar_func), allocatable :: f
    character(len=7) :: label
    character(len=127) :: errmsg

    !! Namelist variables
    integer :: face_block_ids(MAX_FACE_BLOCK_IDS)
    real(r8) :: emissivity_constant
    character(len=MAX_NAME_LEN) :: name, enclosure_name, emissivity_function
    namelist /enclosure_surface/ name, enclosure_name, face_block_ids, &
                                 emissivity_constant, emissivity_function

    call TLS_info ('')
    call TLS_info ('Reading ENCLOSURE_SURFACE namelists ...')

    if (is_IOP) rewind(lun)
    n = 0 ! namelist counter
    stat = 0

    do  ! until all ENCLOSURE_SURFACE namelists have been read or an error occurs.

      n = n + 1
      write(label,'(a,i0,a)') '  [', n, ']'

      if (is_IOP) call seek_to_namelist (lun, 'ENCLOSURE_SURFACE', found, iostat=stat)

      call broadcast (stat)
      if (stat /= 0) then
        errmsg = trim(label) // ' seek error, iostat=' // i_to_c(stat)
        exit
      end if

      call broadcast (found)
      if (.not.found) exit

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        name = NULL_C
        enclosure_name = NULL_C
        face_block_ids = NULL_I
        emissivity_constant = NULL_R
        emissivity_function = NULL_C
        read(lun,nml=enclosure_surface,iostat=stat)
      end if

      call broadcast (stat)
      if (stat /= 0) then
        errmsg = trim(label) // ' read error, iostat=' // i_to_c(stat)
        exit
      end if

      !! Replicate the namelist variables on all processes.
      call broadcast (name)
      call broadcast (enclosure_name)
      call broadcast (face_block_ids)
      call broadcast (emissivity_constant)
      call broadcast (emissivity_function)

      !! Check the user-supplied NAME.
      if (name == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: no value assigned to NAME'
        exit
      else if (name_exists(name)) then
        stat = -1
        errmsg = trim(label) // ' error: another ENCLOSURE_SURFACE has this NAME: ' // trim(name)
        exit
      end if

      !! Verify ENCLOSURE_NAME was assigned a value.
      if (enclosure_name == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: no value assigned to ENCLOSURE_NAME'
        exit
      end if

      !! Check for a non-empty FACE_BLOCK_IDS.
      if (count(face_block_ids /= NULL_I) == 0) then
        stat = -1
        errmsg = trim(label) // ' error: no values assigned to FACE_BLOCK_IDS'
        exit
      end if

      !! Verify that only one of EMISSIVITY_CONSTANT and EMISSIVITY_FUNCTION was specified.
      if (emissivity_constant == NULL_R .eqv. emissivity_function == NULL_C) then
        stat = -1
        errmsg = trim(label) // ' error: exactly one of EMISSIVITY_CONSTANT and &
                                &EMISSIVITY_FUNCTION must be assigned a value'
        exit
      end if

      !! Get or create the specified function.
      if (emissivity_constant /= NULL_R) then
        if (emissivity_constant < 0.0_r8 .or. emissivity_constant > 1.0_r8) then
          stat = -1
          errmsg = trim(label) // ' error: EMISSIVITY_CONSTANT not in [0,1]'
          exit
        else if (emissivity_constant == 0.0_r8) then
          call TLS_info (trim(label) // ' warning: emissivity is 0')
        end if
        call alloc_const_scalar_func (f, emissivity_constant)
      else
        call lookup_func (emissivity_function, f)
        if (.not.allocated(f)) then
          stat = -1
          errmsg = trim(label) // 'error: unknown function name: ' // trim(emissivity_function)
          exit
        end if
      end if

      !! Append the data for this namelist to the list of namelist data.
      if (associated(es_last)) then
        allocate(es_last%next)
        es_last => es_last%next
      else  ! list is empty
        allocate(es_list)
        es_last => es_list
      end if
      es_last%seq = n
      es_last%name = name
      es_last%encl_name = enclosure_name
      call move_alloc (f, es_last%eps)
      allocate(es_last%face_block_ids(count(face_block_ids /= NULL_I)))
      es_last%face_block_ids = pack(face_block_ids, mask=(face_block_ids /= NULL_I))

      call TLS_info (trim(label) // ' read ENCLOSURE_SURFACE namelist "' // trim(name) // '"')

    end do

    if (stat /= 0) then
      call TLS_info (trim(errmsg))
      call TLS_fatal ('error reading ENCLOSURE_SURFACE namelists')
    end if

  contains

    logical function name_exists (name)
      character(len=*), intent(in) :: name
      type(es_list_node), pointer :: l
      l => es_list
      name_exists = .true.
      do while (associated(l))
        if (l%name == name) return
        l => l%next
      end do
      name_exists = .false.
    end function name_exists

  end subroutine ERI_read_enclosure_surface

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ERI_GET_EMISSIVITY
 !!

  subroutine ERI_get_emissivity (encl_name, encl, eps)

    use rad_encl_type
    use rad_encl_func_type

    character(len=*), intent(in)  :: encl_name
    type(rad_encl), pointer :: encl
    type(rad_encl_func), intent(out) :: eps

    integer :: stat
    character(len=127) :: errmsg
    type(es_list_node), pointer :: l

    call TLS_info ('    Defining emissivity for enclosure "' // trim(encl_name) // '" ...')

    call eps%prep (encl)

    l => es_list
    do while (associated(l))
      if (l%encl_name == encl_name) then
        write(errmsg,'(6x,a,i0,2a)') 'using ENCLOSURE_SURFACE[', l%seq, ']: ', trim(l%name)
        call TLS_info (trim(errmsg)) ! not an error message!
        call eps%add_function (l%eps, l%face_block_ids, stat, errmsg)
        call TLS_fatal_if_any (stat /=0 , 'Error defining emissivity: ' // trim(errmsg))
      end if
      l => l%next
    end do

    call eps%done (stat, errmsg)
    call TLS_fatal_if_any (stat /= 0, 'Error defining emissivity: ' // trim(errmsg))

  end subroutine ERI_get_emissivity

end module ER_input
