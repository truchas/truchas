program main

  use,intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use parameter_list_json
  use timer_tree_type
  use th_electrostatics_sim_type
  implicit none

  integer :: n, num_arg, inlun, stat
  character(255) :: arg
  character(:), allocatable :: prog, infile, errmsg
  type(parameter_list), pointer :: params
  type(th_electrostatics_sim) :: sim

  !! Initialize MPI
  call init_parallel_communication
  call fhypre_initialize

  !! Get the program name from the command line.
  call get_command_argument(0, arg)
  n = scan(arg, '/', back=.true.)
  prog = trim(arg(n+1:))  ! remove the leading path component, if any

  !! Get the input file path from the command line.
  num_arg = command_argument_count()
  if (num_arg == 1) then
    call get_command_argument(1, arg)
    infile = trim(arg)
  else
    if (is_IOP) write(error_unit,'(3a)') 'usage: ', prog, ' INFILE'
    call halt_parallel_communication
    error stop
  end if

  !! Read the parameter list from the input file
  open(newunit=inlun,file=infile,action='read',access='stream')
  call parameter_list_from_json_stream(inlun, params, errmsg)
  if (.not.associated(params)) then
    if (is_IOP) write(error_unit,'(a)') 'error reading input file:', errmsg
    call halt_parallel_communication
    error stop
  end if
  close(inlun)

  !! Initialize the message logging system. TLS will write to <prefix>.log
  n = scan(infile, '/', back=.true.)
  prefix = trim(infile(n+1:)) ! remove leading path component, if any
  n = scan(prefix, '.', back=.true.)
  if (n > 0) prefix = prefix(:n-1) ! remove trailing suffix, if any (i.e., '.inp')
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

  !! Create the simulation and run it.
  call start_timer('simulation')
  call sim%init(params, stat, errmsg)
  if (stat == 0) call sim%run(stat, errmsg)
  call stop_timer('simulation')

  !! Write some timing info.
  call TLS_info('')
  call TLS_info('Timing Summary:')
  call TLS_info('')
  if (is_IOP) call write_timer_tree(output_unit, indent=3)

  !! And quit.
  call TLS_info('')
  if (stat == 0) then
    call TLS_exit
  else
    call TLS_fatal(errmsg)
  end if

end program main
