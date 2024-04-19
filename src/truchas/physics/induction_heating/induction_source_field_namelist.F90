module induction_source_field_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_induction_source_field_namelist

contains

  subroutine read_induction_source_field_namelist(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_I, NULL_C, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios, j
    logical :: found
    character(128) :: iom
    character(:), allocatable :: pname, vname
    real(r8), allocatable :: array(:)
    type(parameter_list), pointer :: coils_plist, plist

    !! Namelist variables
    character :: orientation
    integer, parameter :: M = 32
    real(r8) :: times(M), frequency(M+1), uniform_strength(M+1)
    type :: induction_coil
      integer :: num_loops
      real(r8) :: center(3), radius, length, current(M+1)
    end type
    type(induction_coil) :: coils(32)
    namelist /induction_source_field/ orientation, times, frequency, uniform_strength, coils

    call TLS_info('Reading INDUCTION_SOURCE_FIELD namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'induction_source_field', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('induction_source_field namelist not found')

    orientation  = NULL_C
    times = NULL_R
    frequency = NULL_R
    uniform_strength = NULL_R
    call set_coils_default

    if (is_IOP) read(lun,nml=induction_source_field,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading INDUCTION_SOURCE_FIELD namelist: ' // trim(iom))

    call broadcast(orientation)
    call broadcast(times)
    call broadcast(frequency)
    call broadcast(uniform_strength)
    call broadcast_coils

    select case (orientation)
    case ('x','X','y','Y','z','Z')
    case (NULL_C)
      call TLS_info('  using default value "Z" for ORIENTATION')
      orientation = 'Z'
    case default
      call TLS_fatal('invalid ORIENTATION: ' // trim(orientation))
    end select
    call params%set('orientation', orientation)

    array = pack(times, mask=(times /= NULL_R))
    if (size(array) > 1) then
      n = size(array)
      if (any(array(2:n) <= array(:n-1))) &
          call TLS_fatal('TIMES values must be strictly increasing')
      call params%set('times', array)
    else
      n = 0
    end if

    array = pack(frequency, mask=(frequency /= NULL_R))
    if (size(array) == 0) then
      call TLS_fatal('FREQUENCY must be assigned a value')
    else if (any(array <= 0.0_r8)) then
      call TLS_fatal('FREQUENCY values must be > 0.0')
    else if (size(array) /= n+1) then
      call TLS_fatal(i_to_c(n+1) // ' values required for FREQUENCY')
    end if
    call params%set('frequency', array) !TODO: separate scalar and array cases?

    array = pack(uniform_strength, mask=(uniform_strength /= NULL_R))
    if (size(array) > 0) then ! this is optional
      if (size(array) /= n+1) then
        call TLS_fatal(i_to_c(n+1) // ' values required for UNIFORM_STRENGTH')
      end if
      call params%set('uniform-strength', array) !TODO: separate scalar and array cases?
    end if

    coils_plist => params%sublist('coils')
    do j = 1, size(coils)
      if (default_coil(coils(j))) cycle ! unused element
      plist => coils_plist%sublist('COIL' // i_to_c(j))

      pname = 'COILS(' // i_to_c(j) // ')%'
      vname = pname // 'NUM_LOOPS'
      if (coils(j)%num_loops == NULL_I) then
        call TLS_fatal(vname // ' must be assigned a value')
      else if (coils(j)%num_loops <= 0) then
        call TLS_fatal(vname // ' must be > 0')
      end if
      call plist%set('num-loops', coils(j)%num_loops)

      vname = pname // 'RADIUS'
      if (coils(j)%radius == NULL_R) then
        call TLS_fatal(vname // ' must be assigned a value')
      else if (coils(j)%radius <= 0.0_r8) then
        call TLS_fatal(vname // ' must be > 0.0')
      end if
      call plist%set('radius', coils(j)%radius)

      vname = pname // 'LENGTH'
      if (coils(j)%num_loops > 1) then
        if (coils(j)%length == NULL_R) then
          call TLS_fatal(vname // ' must be assigned a value')
        else if (coils(j)%length < 0.0_r8) then
          call TLS_fatal(vname // ' must be >= 0.0')
        end if
        call plist%set('length', coils(j)%length)
      end if

      vname = pname // 'CENTER'
      if (any(coils(j)%center /= NULL_R)) then
        if (any(coils(j)%center == NULL_R)) call TLS_fatal(vname // ' requires 3 values')
        call plist%set('center', coils(j)%center)
      end if

      vname = pname // 'CURRENT'
      array = pack(coils(j)%current, mask=(coils(j)%current/=NULL_R))
      if (size(array) == 0) then
        call TLS_fatal(vname // ' must be assigned a value')
      else if (size(array) /= n+1) then
        call TLS_fatal(i_to_c(n+1) // ' values required for ' // vname)
      end if
      call plist%set('current', array)
    end do

  contains

    subroutine set_coils_default
      integer :: j
      do j = 1, size(coils)
        coils(j)%num_loops = NULL_I
        coils(j)%center = NULL_R
        coils(j)%radius = NULL_R
        coils(j)%length = NULL_R
        coils(j)%current = NULL_R
      end do
    end subroutine

    subroutine broadcast_coils
      integer :: j
      do j = 1, size(coils)
        call broadcast(coils(j)%num_loops)
        call broadcast(coils(j)%center)
        call broadcast(coils(j)%radius)
        call broadcast(coils(j)%length)
        call broadcast(coils(j)%current)
      end do
    end subroutine

    logical function default_coil(coil)
      type(induction_coil), intent(in) :: coil
      default_coil = .false.
      if (coil%num_loops /= NULL_I) return
      if (any(coil%center /= NULL_R)) return
      if (coil%radius /= NULL_R) return
      if (coil%length /= NULL_R) return
      if (any(coil%current /= NULL_R)) return
      default_coil = .true.
    end function

  end subroutine read_induction_source_field_namelist

end module induction_source_field_namelist
