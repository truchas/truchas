!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module diffusion_solver_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use parameter_list_type
  use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
  use string_utilities, only: i_to_c, raise_case
  use truchas_logging_services
  implicit none
  private

  public :: read_diffusion_solver_namelist

contains

  subroutine read_diffusion_solver_namelist(lun, ds_sys_type, num_species, params)

    integer, intent(in) :: lun, ds_sys_type, num_species
    type(parameter_list), intent(inout) :: params

    integer :: ios
    logical :: found
    character(:), allocatable :: label
    character(128) :: iom
    character(len=8)  :: string
    type(parameter_list), pointer :: plist

    !! Namelist variables

    !! Time stepping control parameters.
    integer :: pc_freq, max_step_tries
    logical :: verbose_stepping
    integer :: integrator
    integer, parameter :: DS_ADAPTIVE_BDF2 = 1
    integer, parameter :: DS_NONADAPTIVE_BDF1 = 2

    !! NLK nonlinear solver parameters for BDF2 step.
    integer  :: max_nlk_itr, max_nlk_vec
    real(r8) :: nlk_tol, nlk_vec_tol

    !! Preconditioner choices for the NLK nonlinear solver.
    integer :: ds_nlk_pc = 0
    integer, parameter :: DS_NLK_PC_SSOR = 1
    integer, parameter :: DS_NLK_PC_HYPRE_AMG = 2

    !! Preconditioner parameters
    integer  :: pc_ssor_sweeps, pc_amg_cycles
    real(r8) :: pc_ssor_relax, hypre_amg_strong_threshold
    integer  :: hypre_amg_print_level, hypre_amg_debug_level, hypre_amg_logging_level
    integer  :: hypre_amg_coarsen_type, hypre_amg_interp_type
    integer  :: hypre_amg_relax_down_type, hypre_amg_relax_up_type

    !! Error tolerances
    real(r8) :: abs_conc_tol, abs_temp_tol, abs_enthalpy_tol
    real(r8) :: rel_conc_tol, rel_temp_tol, rel_enthalpy_tol
    real(r8) :: residual_atol, residual_rtol

    real(r8) :: cond_vfrac_threshold
    character(len=32) :: nlk_preconditioner, stepping_method

    logical :: use_new_mfd  ! use the new MFD mass matrix (temporary)
    logical :: hypre_amg_debug
    real(r8) :: void_temperature

    namelist /diffusion_solver/ max_step_tries, abs_conc_tol, rel_conc_tol, &
        abs_temp_tol, rel_temp_tol, abs_enthalpy_tol, rel_enthalpy_tol, &
        max_nlk_itr, nlk_tol, max_nlk_vec, nlk_vec_tol, pc_freq, &
        pc_ssor_sweeps, pc_ssor_relax, verbose_stepping, stepping_method, &
        nlk_preconditioner, pc_amg_cycles, hypre_amg_print_level, &
        hypre_amg_debug, hypre_amg_logging_level, &
        hypre_amg_coarsen_type, hypre_amg_interp_type, hypre_amg_strong_threshold, &
        hypre_amg_relax_down_type, hypre_amg_relax_up_type, &
        cond_vfrac_threshold, residual_atol, residual_rtol, &
        use_new_mfd, void_temperature

    integer, parameter :: DS_SPEC_SYS = 1
    integer, parameter :: DS_TEMP_SYS = 2
    integer, parameter :: DS_TEMP_SPEC_SYS = 3

    call TLS_info('Reading DIFFUSION_SOLVER namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'diffusion_solver', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('DIFFUSION_SOLVER namelist not found')

    !! Assign default values to the namelist variables before reading
    max_nlk_itr    = NULL_I
    nlk_tol        = NULL_R
    max_nlk_vec    = NULL_I
    nlk_vec_tol    = NULL_R
    max_step_tries = NULL_I
    abs_conc_tol   = NULL_R
    rel_conc_tol   = NULL_R
    abs_temp_tol   = NULL_R
    rel_temp_tol   = NULL_R
    abs_enthalpy_tol = NULL_R
    rel_enthalpy_tol = NULL_R
    nlk_preconditioner = NULL_C
    pc_freq = NULL_I
    pc_ssor_sweeps = NULL_I
    pc_ssor_relax  = NULL_R
    pc_amg_cycles  = NULL_I
    hypre_amg_print_level = NULL_I
    hypre_amg_logging_level = NULL_I
    hypre_amg_debug = .false.
    hypre_amg_coarsen_type = NULL_I
    hypre_amg_interp_type = NULL_I
    hypre_amg_strong_threshold = NULL_R
    hypre_amg_relax_down_type = NULL_I
    hypre_amg_relax_up_type = NULL_I
    verbose_stepping = .false.
    stepping_method = NULL_C
    cond_vfrac_threshold = NULL_R
    residual_atol = NULL_R
    residual_rtol = NULL_R
    use_new_mfd = .true.
    void_temperature = 0.0_r8

    if (is_IOP) read(lun,nml=diffusion_solver,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading DIFFUSION_SOLVER namelist: ' // trim(iom))

    !! Replicate the namelist variables on all processes
    call broadcast(max_nlk_itr)
    call broadcast(nlk_tol)
    call broadcast(max_nlk_vec)
    call broadcast(nlk_vec_tol)
    call broadcast(max_step_tries)
    call broadcast(abs_conc_tol)
    call broadcast(rel_conc_tol)
    call broadcast(abs_temp_tol)
    call broadcast(rel_temp_tol)
    call broadcast(abs_enthalpy_tol)
    call broadcast(rel_enthalpy_tol)
    call broadcast(nlk_preconditioner)
    call broadcast(pc_freq)
    call broadcast(pc_ssor_sweeps)
    call broadcast(pc_ssor_relax)
    call broadcast(pc_amg_cycles)
    call broadcast(hypre_amg_print_level)
    call broadcast(hypre_amg_logging_level)
    call broadcast(hypre_amg_debug)
    call broadcast(hypre_amg_coarsen_type)
    call broadcast(hypre_amg_interp_type)
    call broadcast(hypre_amg_strong_threshold)
    call broadcast(hypre_amg_relax_down_type)
    call broadcast(hypre_amg_relax_up_type)
    call broadcast(verbose_stepping)
    call broadcast(stepping_method)
    call broadcast(cond_vfrac_threshold)
    call broadcast(residual_atol)
    call broadcast(residual_rtol)
    call broadcast(use_new_mfd)
    call broadcast(void_temperature)

    call params%set('void-temperature', void_temperature)

    if (stepping_method == NULL_C) then
      stepping_method = 'Adaptive BDF2'
      call TLS_info ('  using default STEPPING_METHOD value: "' // trim(stepping_method) // '"')
    end if
    select case (stepping_method)
    case ('Adaptive BDF2')
      integrator = DS_ADAPTIVE_BDF2
      call params%set('integrator', 'adaptive-bdf2')
    case ('Non-adaptive BDF1')
      integrator = DS_NONADAPTIVE_BDF1
      call params%set('integrator', 'nonadaptive-bdf1')
    case default
      call TLS_fatal ('unknown value for STEPPING_METHOD: "' // trim(stepping_method) // '"')
    end select

    if (max_nlk_itr == NULL_I) then
      max_nlk_itr = 5
      call TLS_info ('  using default MAX_NLK_ITR value: ' // i_to_c(max_nlk_itr))
    else if (max_nlk_itr < 2) then
      call TLS_fatal ('MAX_NLK_ITR must be > 1')
    end if
    call params%set('nlk-max-iter', max_nlk_itr)

    call params%set('verbose-stepping', verbose_stepping)

    select case (integrator)
    case (DS_ADAPTIVE_BDF2)
      if (nlk_tol == NULL_R) then
        nlk_tol = 0.1_r8
        write(string,'(es8.2)') nlk_tol
        call TLS_info ('  using default NLK_TOL value: ' // string)
      else if (nlk_tol <= 0.0_r8 .or. nlk_tol >= 1.0_r8) then
        call TLS_fatal ('NLK_TOL must be positive and < 1')
      end if
      call params%set('nlk-tol', nlk_tol)
      if (max_step_tries == NULL_I) then
        max_step_tries = 10
        call TLS_info ('  using default MAX_STEP_TRIES value: ' // i_to_c(max_step_tries))
      else if (max_step_tries < 1) then
        call TLS_fatal ('MAX_STEP_TRIES must be > 0')
      end if
      call params%set('max-step-tries', max_step_tries)
      if (pc_freq == NULL_I) then
        pc_freq = huge(pc_freq)
      else if (pc_freq < 1) then
        call TLS_fatal ('PC_FREQ must be > 0')
      end if
      call params%set('pc-freq', pc_freq)
      call params%set('hmin', tiny(1.0_r8)) ! HMIN enforced by cycle driver
    case default
      if (nlk_tol /= NULL_R) then
        call TLS_info ('  ignoring NLK_TOL value; not relevant to STEPPING_METHOD choice.')
      end if
      if (max_step_tries /= NULL_I) then
        call TLS_info ('  ignoring MAX_STEP_TRIES value; not relevant to STEPPING_METHOD choice.')
      end if
    end select

    !FIXME: max_nlk_vec already handled above
    if (max_nlk_vec == NULL_I) then
      max_nlk_vec = max_nlk_itr-1
      call TLS_info ('  using default MAX_NLK_VEC value: MAX_NLK_ITR - 1')
    else if (max_nlk_vec < 1) then
      call TLS_fatal ('MAX_NLK_VEC must be > 0')
    else if (max_nlk_vec >= max_nlk_itr) then
      max_nlk_vec = max_nlk_itr - 1
      call TLS_info ('  reducing MAX_NLK_VEC to MAX_NLK_ITR - 1')
    end if
    call params%set('nlk-max-vec', max_nlk_vec)

    if (nlk_vec_tol == NULL_R) then
      nlk_vec_tol = 1.0d-3
      write(string,'(es8.2)') nlk_vec_tol
      call TLS_info ('  using default NLK_VEC_TOL value: ' // string)
    else if (nlk_vec_tol <= 0.0_r8) then
      call TLS_fatal ('NLK_VEC_TOL must be > 0')
    end if
    call params%set('nlk-vec-tol', nlk_vec_tol)

    select case (integrator)
    case (DS_ADAPTIVE_BDF2)

      plist => params%sublist('norm')
      call plist%set('verbose', verbose_stepping)

      select case (ds_sys_type)
      case (DS_SPEC_SYS, DS_TEMP_SPEC_SYS)

        if (abs_conc_tol == NULL_R) then
          call TLS_fatal ('ABS_CONC_TOL must be assigned a value')
        else if (abs_conc_tol < 0.0_r8) then
          call TLS_fatal ('ABS_CONC_TOL must be >= 0')
        end if

        if (rel_conc_tol == NULL_R) then
          rel_conc_tol = 0.0_r8
          write(string,'(es8.2)') rel_conc_tol
          call TLS_info ('  using default REL_CONC_TOL value: ' // string)
        else if (rel_conc_tol < 0.0_r8 .or. rel_conc_tol >= 1.0_r8) then
          call TLS_fatal ('REL_CONC_TOL must be in [0,1)')
        end if

        if (abs_conc_tol == 0.0_r8) then
          if (rel_conc_tol == 0.0_r8) then
            call TLS_fatal ('ABS_CONC_TOL and REL_CONC_TOL are both 0')
          else
            call TLS_info ('  WARNING: using a pure relative error norm; conc must be bounded away from 0')
          end if
        end if

        call plist%set('abs-c-tol', spread(abs_conc_tol, dim=1, ncopies=num_species))
        call plist%set('rel-c-tol', spread(rel_conc_tol, dim=1, ncopies=num_species))

      case default  ! system doesn't involve species concentration

        if (abs_conc_tol /= NULL_R) then
          call TLS_info ('  ignoring ABS_CONC_TOL value; concentration is not a dependent variable')
        end if

        if (rel_conc_tol /= NULL_R) then
          call TLS_info ('  ignoring REL_CONC_TOL value; concentration is not a dependent variable')
        end if

      end select

      select case (ds_sys_type)
      case (DS_TEMP_SYS, DS_TEMP_SPEC_SYS)

        if (abs_temp_tol == NULL_R) then
          call TLS_fatal ('ABS_TEMP_TOL must be assigned a value')
        else if (abs_temp_tol < 0.0_r8) then
          call TLS_fatal ('ABS_TEMP_TOL must be >= 0')
        end if

        if (rel_temp_tol == NULL_R) then
          rel_temp_tol = 0.0_r8
          write(string,'(es8.2)') rel_temp_tol
          call TLS_info ('  using default REL_TEMP_TOL value: ' // string)
        else if (rel_temp_tol < 0.0_r8 .or. rel_temp_tol >= 1.0_r8) then
          call TLS_fatal ('REL_TEMP_TOL must be in [0,1)')
        end if

        if (abs_temp_tol == 0.0_r8) then
          if (rel_temp_tol == 0.0_r8) then
            call TLS_fatal ('ABS_TEMP_TOL and Rel_Temp_Tol are both 0')
          else
            call TLS_info ('  WARNING: using a pure relative error norm; temp must be bounded away from 0')
          end if
        end if

        call plist%set('abs-t-tol', abs_temp_tol)
        call plist%set('rel-t-tol', rel_temp_tol)

        if (abs_enthalpy_tol == NULL_R) then
          call TLS_fatal ('ABS_ENTHALPY_TOL must be assigned a value')
        else if (abs_enthalpy_tol < 0.0_r8) then
          call TLS_fatal ('ABS_ENTHALPY_TOL must be >= 0')
        end if

        if (rel_enthalpy_tol == NULL_R) then
          rel_enthalpy_tol = 0.0_r8
          write(string,'(es8.2)') rel_enthalpy_tol
          call TLS_info ('  using default REL_ENTHALPY_TOL value: ' // string)
        else if (rel_enthalpy_tol < 0.0_r8 .or. rel_enthalpy_tol >= 1.0_r8) then
          call TLS_fatal ('REL_ENTHALPY_TOL must be in [0,1)')
        end if

        if (abs_enthalpy_tol == 0.0_r8) then
          if (rel_enthalpy_tol == 0.0_r8) then
            call TLS_fatal ('ABS_ENTHALPY_TOL and REL_ENTHALPY_TOL are both 0')
          else
            call TLS_info ('  WARNING: using a pure relative error norm; enthalpy must be bounded away from 0')
          end if
        end if

        call plist%set('abs-h-tol', abs_enthalpy_tol)
        call plist%set('rel-h-tol', rel_enthalpy_tol)

      case default  ! system doesn't involve temperature

        if (abs_temp_tol /= NULL_R) then
          call TLS_info ('  ignoring ABS_TEMP_TOL value; temperature is not a dependent variable')
        end if

        if (rel_temp_tol /= NULL_R) then
          call TLS_info ('  ignoring REL_TEMP_TOL value; temperature is not a dependent variable')
        end if

        if (abs_enthalpy_tol /= NULL_R) then
          call TLS_info ('  ignoring ABS_ENTHALPY_TOL value; enthalpy is not a dependent variable')
        end if

        if (rel_enthalpy_tol /= NULL_R) then
          call TLS_info ('  ignoring REL_ENTHALPY_TOL value; enthalpy is not a dependent variable')
        end if

      end select

    case default

      if (abs_conc_tol /= NULL_R) then
        call TLS_info ('  ignoring ABS_CONC_TOL value; not relevant to STEPPING_METHOD choice')
      end if

      if (rel_conc_tol /= NULL_R) then
        call TLS_info ('  ignoring REL_CONC_TOL value; not relevant to STEPPING_METHOD choice')
      end if

      if (abs_temp_tol /= NULL_R) then
        call TLS_info ('  ignoring ABS_TEMP_TOL value; not relevant to STEPPING_METHOD choice')
      end if

      if (rel_temp_tol /= NULL_R) then
        call TLS_info ('  ignoring REL_TEMP_TOL value; not relevant to STEPPING_METHOD choice')
      end if

      if (abs_enthalpy_tol /= NULL_R) then
        call TLS_info ('  ignoring ABS_ENTHALPY_TOL value; not relevant to STEPPING_METHOD choice')
      end if

      if (rel_enthalpy_tol /= NULL_R) then
        call TLS_info ('  ignoring REL_ENTHALPY_TOL value; not relevant to STEPPING_METHOD choice')
      end if

    end select

    plist => params%sublist('precon')

    nlk_preconditioner = raise_case(adjustl(nlk_preconditioner))
    select case (nlk_preconditioner)
    case ('SSOR')
      ds_nlk_pc = DS_NLK_PC_SSOR
      call plist%set('method', 'SSOR')
    case ('HYPRE_AMG')
      ds_nlk_pc = DS_NLK_PC_HYPRE_AMG
      call plist%set('method', 'BoomerAMG')
    case (NULL_C)
      ds_nlk_pc = DS_NLK_PC_SSOR
    case default
      call TLS_fatal ('unknown value for NLK_PRECONDITIONER: ' // trim(nlk_preconditioner))
    end select

    plist => plist%sublist('params')

    if (ds_nlk_pc == DS_NLK_PC_SSOR) then

      if (pc_ssor_sweeps == NULL_I) then
        pc_ssor_sweeps = 4
        call TLS_info ('  using default PC_SSOR_SWEEPS value: ' // i_to_c(pc_ssor_sweeps))
      else if (pc_ssor_sweeps < 1) then
        call TLS_fatal ('PC_SSOR_SWEEPS must be > 0')
      end if
      call plist%set('num-cycles', pc_ssor_sweeps)

      if (pc_ssor_relax == NULL_R) then
        pc_ssor_relax = 1.0_r8
        write(string,'(es8.2)') pc_ssor_relax
        call TLS_info ('  using default PC_SSOR_RELAX value: ' // string)
      else if (pc_ssor_relax <= 0.0_r8 .or. pc_ssor_relax >= 2.0_r8) then
        call TLS_fatal ('PC_SSOR_RELAX must be in (0, 2)')
      end if
      call plist%set('omega', pc_ssor_relax)

    else

      if (pc_ssor_sweeps /= NULL_I) then
        call TLS_info ('  ignoring PC_SSOR_SWEEPS value; unused when NLK_PRECONDITIONER /= "SSOR"')
      end if

      if (pc_ssor_relax /= NULL_R) then
        call TLS_info ('  ignoring PC_SSOR_RELAX value; unused when NLK_PRECONDITIONER /= "SSOR"')
      end if

    end if

    if (ds_nlk_pc == DS_NLK_PC_HYPRE_AMG) then

      if (pc_amg_cycles == NULL_I) then
        pc_amg_cycles = 2
        call TLS_info ('  using default PC_AMG_CYCLES value: ' // i_to_c(pc_amg_cycles))
      else if (pc_amg_cycles < 1) then
        call TLS_fatal ('PC_AMG_CYCLES must be > 0')
      endif
      call plist%set('num-cycles', pc_amg_cycles)

      if (hypre_amg_coarsen_type /= NULL_I) &
        call plist%set('coarsen-type', hypre_amg_coarsen_type)

      if (hypre_amg_interp_type /= NULL_I) &
        call plist%set('interp-type', hypre_amg_interp_type)

      if (hypre_amg_relax_down_type /= NULL_I) &
        call plist%set('relax-down-type', hypre_amg_relax_down_type)

      if (hypre_amg_relax_up_type /= NULL_I) &
        call plist%set('relax-up-type', hypre_amg_relax_up_type)

      if (hypre_amg_strong_threshold /= NULL_R) &
        call plist%set('strong_threshold', hypre_amg_strong_threshold)

      if (hypre_amg_print_level == NULL_I) then
        hypre_amg_print_level = 0
      else if (hypre_amg_print_level < 0 .or. hypre_amg_print_level > 4) then
        call TLS_fatal ('HYPRE_AMG_PRINT_LEVEL must be >= 0 and <= 3')
      end if
      call plist%set('print-level', hypre_amg_print_level)

      if (hypre_amg_logging_level == NULL_I) then
        hypre_amg_logging_level = 0
      else if (hypre_amg_logging_level < 0) then
        call TLS_fatal ('HYPRE_AMG_LOGGING_LEVEL must be >= 0')
      end if
      call plist%set('logging-level', hypre_amg_logging_level)

      hypre_amg_debug_level = 0
      if (hypre_amg_debug) hypre_amg_debug_level = 1
      call plist%set('debug-level', hypre_amg_debug_level)

    else

      if (pc_amg_cycles /= NULL_I) then
        call TLS_info ('  ignoring PC_AMG_CYCLES value; unused when NLK_PRECONDITIONER /= "HYPRE_AMG"')
      end if

      if (hypre_amg_print_level /= NULL_I) then
        call TLS_info ('  ignoring HYPRE_AMG_PRINT_LEVEL value; unused when NLK_PRECONDITIONER /= "HYPRE_AMG"')
      end if

      if (hypre_amg_logging_level /= NULL_I) then
        call TLS_info ('  ignoring HYPRE_AMG_LOGGING_LEVEL value; unused when NLK_PRECONDITIONER /= "HYPRE_AMG"')
      end if

    end if

    if (integrator == DS_NONADAPTIVE_BDF1) then
      plist => params%sublist('norm')
      call plist%set('verbose', verbose_stepping)
      if (residual_atol == NULL_R) then
        residual_atol = 0.0d0
        call TLS_info ('  using default RESIDUAL_ATOL value: 0.0')
      else if (residual_atol < 0.0d0) then
        call TLS_fatal ('RESIDUAL_ATOL must be >= 0')
      end if
      call plist%set('abs-tol', residual_atol)
      if (residual_rtol == NULL_R) then
        call TLS_fatal ('RESIDUAL_RTOL must be assigned a value')
      else if (residual_rtol <= 0.0d0 .or. residual_rtol >= 1.0d0) then
        call TLS_fatal ('RESIDUAL_RTOL must be > 0 and < 1')
      end if
      call plist%set('rel-tol', residual_rtol)
      if (cond_vfrac_threshold == NULL_R) then
        cond_vfrac_threshold = 1.0d-3
      else if (cond_vfrac_threshold <= 0.0d0 .or. cond_vfrac_threshold >= 1.0d0) then
        call TLS_fatal ('COND_VFRAC_THRESHOLD must be > 0 and < 1')
      end if
      call params%set('epsilon', cond_vfrac_threshold)
    else
      if (residual_atol /= NULL_R) then
        call TLS_info ('  ignoring RESIDUAL_ATOL value; not relevant to STEPPING_METHOD choice')
      end if
      if (residual_rtol /= NULL_R) then
        call TLS_info ('  ignoring RESIDUAL_RTOL value; not relevant to STEPPING_METHOD choice')
      end if
      if (cond_vfrac_threshold /= NULL_R) then
        call TLS_info ('  ignoring COND_VFRAC_THRESHOLD value; not relevant to STEPPING_METHOD choice')
      end if
    end if

  end subroutine read_diffusion_solver_namelist

end module diffusion_solver_namelist
