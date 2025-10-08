module fdme_solver_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_fdme_solver_namelist

contains

  subroutine read_fdme_solver_namelist(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: ios
    logical :: found
    character(128) :: iom
    type(parameter_list), pointer :: plist, sublist

    real(r8) :: pc_damping
    logical :: use_mixed_form, graphics_output
    character(32) :: solver_type, precon_type
    namelist /fdme_solver/ use_mixed_form, graphics_output
    namelist /fdme_solver/ solver_type, precon_type, pc_damping

    !! Linear solver variables
    integer :: max_iter, print_level
    real(r8) :: rel_tol
    namelist /fdme_solver/ rel_tol, max_iter, print_level

    !! Built-in SSOR preconditioner
    integer  :: ssor_num_cycles
    real(r8) :: ssor_omega
    namelist /fdme_solver/ ssor_num_cycles, ssor_omega

    !! HYPRE BoomerAMG preconditioner variables
    integer  :: boomer_num_cycles
    real(r8) :: boomer_strong_threshold
    integer  :: boomer_coarsen_type, boomer_interp_type
    integer  :: boomer_relax_down_type, boomer_relax_up_type
    integer  :: boomer_print_level, boomer_debug_level, boomer_logging_level
    namelist /fdme_solver/ boomer_num_cycles, boomer_strong_threshold, boomer_coarsen_type, &
        boomer_interp_type, boomer_relax_down_type, boomer_relax_up_type, boomer_print_level, &
        boomer_debug_level, boomer_logging_level

    call TLS_info('Reading FDME_SOLVER namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'fdme_solver', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('FDME_SOLVER namelist not found')

    use_mixed_form = .false.
    graphics_output = .false.

    solver_type = NULL_C
    precon_type = NULL_C
    pc_damping = NULL_R

    rel_tol = NULL_R
    max_iter = NULL_I
    print_level = NULL_I

    boomer_num_cycles = NULL_I
    boomer_strong_threshold = NULL_R
    boomer_coarsen_type = NULL_I
    boomer_interp_type = NULL_I
    boomer_relax_down_type = NULL_I
    boomer_relax_up_type = NULL_I
    boomer_print_level = NULL_I
    boomer_debug_level = NULL_I
    boomer_logging_level = NULL_I
    ssor_num_cycles = NULL_I
    ssor_omega = NULL_R

    if (is_IOP) read(lun,nml=fdme_solver,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading FDME_SOLVER namelist: ' // trim(iom))

    call broadcast(use_mixed_form)
    call broadcast(graphics_output)

    call broadcast(solver_type)
    call broadcast(precon_type)
    call broadcast(pc_damping)

    call broadcast(rel_tol)
    call broadcast(max_iter)
    call broadcast(print_level)

    call broadcast(boomer_num_cycles)
    call broadcast(boomer_strong_threshold)
    call broadcast(boomer_coarsen_type)
    call broadcast(boomer_interp_type)
    call broadcast(boomer_relax_down_type)
    call broadcast(boomer_relax_up_type)
    call broadcast(boomer_print_level)
    call broadcast(boomer_debug_level)
    call broadcast(boomer_logging_level)
    call broadcast(ssor_num_cycles)
    call broadcast(ssor_omega)

    call params%set('use-mixed-form', use_mixed_form)
    call params%set('graphics-output', graphics_output)

    select case (solver_type)
    case ('minres')
    case ('mumps')
#ifndef USE_MUMPS
      call TLS_fatal('SOLVER_TYPE = "mumps" is not supported by this Truchas build')
#endif
    case (NULL_C)
      solver_type = 'minres'
    case default
      call TLS_fatal('invalid SOLVER_TYPE: ' // solver_type)
    end select
    call params%set('solver-type', solver_type)

    select case (solver_type)
    case ('minres')
      if (rel_tol /= NULL_R) call params%set('rel-tol', rel_tol)
      if (max_iter /= NULL_I) call params%set('max-iter', max_iter)
      if (print_level /= NULL_I) call params%set('print-level', print_level)
      if (pc_damping /= NULL_R) call plist%set('beta', pc_damping)
      plist => params%sublist('precon')
      select case (precon_type)
      case ('none')
      case ('boomer')
        sublist => plist%sublist('params')
        if (boomer_num_cycles /= NULL_I)  call sublist%set('num-cycles', boomer_num_cycles)
        if (boomer_strong_threshold /= NULL_R)  call sublist%set('strong-threshold', boomer_strong_threshold)
        if (boomer_coarsen_type /= NULL_I)  call sublist%set('coarsen-type', boomer_coarsen_type)
        if (boomer_interp_type /= NULL_I)  call sublist%set('interp-type', boomer_interp_type)
        if (boomer_relax_down_type /= NULL_I)  call sublist%set('relax-down-type', boomer_relax_down_type)
        if (boomer_relax_up_type /= NULL_I)  call sublist%set('relax-up-type', boomer_relax_up_type)
        if (boomer_print_level /= NULL_I)  call sublist%set('print-level', boomer_print_level)
        if (boomer_debug_level /= NULL_I)  call sublist%set('debug-level', boomer_debug_level)
        if (boomer_logging_level /= NULL_I)  call sublist%set('logging-level', boomer_logging_level)
      case ('ssor')
        sublist => plist%sublist('params')
        if (ssor_omega /= NULL_R)  call sublist%set('omega', ssor_omega)
        if (ssor_num_cycles /= NULL_I) call sublist%set('num-cycles', ssor_num_cycles)
      case (NULL_C)
        call TLS_fatal('PRECON_TYPE not specified')
      case default
        call TLS_fatal('invalid PRECON_TYPE: ' // precon_type)
      end select
      call plist%set('type', precon_type)
      call plist%set('method', precon_type) ! for pcsr_precon_factory (FIXME)
      select case (solver_type)
      case ('minres')
      end select
    end select

  end subroutine read_fdme_solver_namelist

end module fdme_solver_namelist
