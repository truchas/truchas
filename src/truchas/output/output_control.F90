!!
!! OUTPUT_CONTROL
!!
!! This module was created to hold the variables from OUTPUT_DATA_MODULE and
!! OUTPUT_MODULE that control when the solution was output.  Everything else
!! in those modules is tied to TBrook and they are slated for deletion when
!! TBrook is ultimately deleted from the code.  I expect this module to exist
!! only temporarily until the control of output is refactored/reimplemented.
!!
!! Neil Carlson <nnc@lanl.gov>
!! October 2012
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module output_control

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use toolpath_type
  implicit none
  private

  public :: output_init, add_output_events, read_outputs_namelist

  integer, parameter :: MAX_TIMES = 21
  integer, parameter :: MAX_FIELDS = 128
  integer, parameter :: MAX_FIELD_NAME_LEN = 64
  real(r8), save, public :: output_t(MAX_TIMES)
  integer,  save, public :: output_n(MAX_TIMES-1)
  integer,  save, public :: ntimes

  integer, allocatable, public :: part(:)
  type(toolpath), allocatable, target, public :: part_path

  type(parameter_list), public :: params

contains

  subroutine output_init(t)
    real(r8), intent(in) :: t
    if (allocated(part_path)) call part_path%set_segment(t)
  end subroutine

  subroutine add_output_events(eventq)
    use sim_event_queue_type
    use toolpath_event_type
    type(sim_event_queue), intent(inout) :: eventq
    integer :: j
    if (allocated(part_path)) then
      block
        real(r8), allocatable :: times(:)
        call part_path%get_segment_starts(times)
        do j = 1, size(times)
          call eventq%add_event(times(j), toolpath_event(part_path))
        end do
      end block
    end if
  end subroutine


  subroutine read_outputs_namelist(lun)

    use input_utilities, only: seek_to_namelist, NULL_I, NULL_C, NULL_R
    use string_utilities, only: i_to_c
    use parallel_communication,  only: is_IOP, broadcast
    use truchas_logging_services
    use truchas_env, only: output_file_name, overwrite_output

    integer, intent(in) :: lun

    integer :: ios
    logical :: found, exists
    character(128) :: iom
    character(:), allocatable :: filename
    type(parameter_list), pointer :: stream, plist

    !! Namelist variables
    integer :: move_block_ids(32)
    character(32) :: move_toolpath_name
    character(MAX_FIELD_NAME_LEN) :: fields(MAX_FIELDS)
    integer :: subintervals(MAX_TIMES-1)
    real(r8) :: times(MAX_TIMES)
    namelist /outputs/ times, subintervals, move_block_ids, move_toolpath_name, &
        fields

    call TLS_info('Reading OUTPUTS namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'outputs', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('OUTPUTS namelist not found')

    times = NULL_R
    subintervals = NULL_I
    move_block_ids = NULL_I
    move_toolpath_name = NULL_C
    fields = NULL_C

    if (is_IOP) read(lun,nml=outputs,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading OUTPUTS namelist: ' // trim(iom))

    call broadcast(times)
    call broadcast(subintervals)
    call broadcast(move_block_ids)
    call broadcast(move_toolpath_name)
    call broadcast(fields)

    if (all(times == NULL_R)) call TLS_fatal('TIMES not specified')

    !! Determine the number of specified TIME values
    ntimes = findloc(times, NULL_R, dim=1)
    ntimes = modulo(ntimes-1, 1+size(times))
    if (any(times(ntimes+1:) /= NULL_R)) call TLS_fatal('invalid TIMES specification')

    !! Require at least 2 values and the values must be strictly increasing.
    if (ntimes < 2) call TLS_fatal('invalid TIMES specification')
    if (any(times(2:ntimes)-times(:ntimes-1) <= 0)) call TLS_fatal('invalid TIMES specification')

    stream => params%sublist('streams')
    ! The OUTPUTS namelist currently configures a single VTKHDF stream named "main".
    stream => stream%sublist('main')

    plist => stream%sublist('schedule')
    call plist%set('times', times(:ntimes))
    output_t = 0
    output_t(:ntimes) = times(:ntimes)

    output_n = 1
    if (any(subintervals /= NULL_I)) then
      if (any(subintervals(:ntimes-1) == NULL_I) .or. any(subintervals(:ntimes-1) < 1) .or. &
          any(subintervals(ntimes:) /= NULL_I)) &
          call TLS_fatal('invalid SUBINTERVALS specification')
      call plist%set('subintervals', subintervals(:ntimes-1))
      output_n(:ntimes-1) = subintervals(:ntimes-1)
    end if

    if (any(move_block_ids /= NULL_I)) then
      if (move_toolpath_name == NULL_C) call TLS_fatal('MOVE_TOOLPATH_NAME not specified')
      call stream%set('move-block-ids', pack(move_block_ids, mask=(move_block_ids/=NULL_I)))
      call stream%set('move-toolpath-name', trim(move_toolpath_name))
    else if (move_toolpath_name /= NULL_C) then
      call TLS_fatal('MOVE_BLOCK_IDS not specified')
    end if

    allocate(part(count(move_block_ids/=NULL_I)))
    part = pack(move_block_ids,mask=move_block_ids/=NULL_I)
    if (move_toolpath_name /= NULL_C) then
      block
        use toolpath_driver, only: alloc_toolpath
        integer :: stat
        character(:), allocatable :: errmsg
        call alloc_toolpath(part_path, move_toolpath_name, stat, errmsg)
        if (stat /= 0) call TLS_fatal('unable to create toolpath: ' // errmsg)
      end block
    end if

    filename = output_file_name('vtkhdf')
    if (is_IOP .and. .not. overwrite_output) then
      inquire(file=filename, exist=exists)
      if (exists) call TLS_panic('must specify "-f" flag to overwrite "' // filename // '"')
    endif
    call stream%set('filename', filename)

    if (any(fields /= NULL_C)) &
        call stream%set('fields', pack(fields, mask=(fields /= NULL_C)))

  end subroutine read_outputs_namelist

end module output_control
