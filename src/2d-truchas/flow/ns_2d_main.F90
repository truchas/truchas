!!
!! NS_2D_MAIN
!!
!! A basic JSON-input driver for the two-dimensional isothermal incompressible
!! Navier--Stokes model.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

program ns_2d_main

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use,intrinsic :: iso_fortran_env, only: error_unit
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use parameter_list_json
  use ns_2d_sim_type
  implicit none

  integer :: n, num_arg, inlun, stat
  character(255) :: arg
  character(:), allocatable :: prog, infile, errmsg
  type(parameter_list), pointer :: params
  type(ns_2d_sim) :: sim

  call init_parallel_communication
  call fhypre_initialize
  call get_command_argument(0, arg)
  n = scan(arg, '/', back=.true.)
  prog = trim(arg(n+1:))
  num_arg = command_argument_count()
  if (num_arg /= 1) then
    if (is_IOP) write(error_unit,'(3a)') 'usage: ', prog, ' INFILE'
    call halt_parallel_communication
    stop 1
  end if
  call get_command_argument(1, arg)
  infile = trim(arg)
  open(newunit=inlun, file=infile, action='read', access='stream')
  call parameter_list_from_json_stream(inlun, params, errmsg)
  close(inlun)
  if (.not.associated(params)) then
    if (is_IOP) write(error_unit,'(a)') 'error reading input file: ' // errmsg
    call halt_parallel_communication
    call exit(1)
  end if

  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)
  call sim%init(params, stat, errmsg)
  if (stat == 0) call sim%run(stat, errmsg)
  if (stat /= 0) then
    if (is_IOP) write(error_unit,'(a)') 'flow simulation error: ' // errmsg
    call TLS_exit
    call halt_parallel_communication
    call exit(1)
  end if
  call TLS_exit
  call halt_parallel_communication

end program ns_2d_main
