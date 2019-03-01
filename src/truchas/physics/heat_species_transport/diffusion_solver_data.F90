!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module diffusion_solver_data

  use kinds
  use truchas_logging_services
  implicit none
  private
  save

  public :: read_ds_namelist

  !! These variables are defined manually prior to calling READ_DS_NAMELIST
  logical, public :: ds_enabled = .false.
  character(32), public :: system_type = ''
  integer,  public :: num_species = 0

  !! The remaining variables are defined by READ_DS_NAMELIST.
  integer, public :: ds_sys_type = 0
  integer, parameter, public :: DS_SPEC_SYS = 1
  integer, parameter, public :: DS_TEMP_SYS = 2
  integer, parameter, public :: DS_TEMP_SPEC_SYS = 3
  
  logical, public :: heat_eqn = .false.
  integer, public :: nconc = 0

  character(4), parameter, public :: mesh_name = 'MAIN'

  !! Time stepping control parameters.
  integer,  public :: max_step_tries
  real(r8), public :: hmin = tiny(1.0_r8)
  logical,  public :: verbose_stepping
  integer,  public :: integrator = 0
  integer, parameter, public :: DS_ADAPTIVE_BDF2 = 1
  integer, parameter, public :: DS_NONADAPTIVE_BDF1 = 2

  !! NLK nonlinear solver parameters for BDF2 step.
  integer,  public :: max_nlk_itr
  integer,  public :: max_nlk_vec
  real(r8), public :: nlk_tol
  real(r8), public :: nlk_vec_tol

  !! Preconditioner choices for the NLK nonlinear solver.
  integer, public :: ds_nlk_pc = 0
  integer, parameter, public :: DS_NLK_PC_SSOR = 1
  integer, parameter, public :: DS_NLK_PC_HYPRE_AMG = 2

  !! Preconditioner parameters
  integer,  public :: pc_ssor_sweeps
  real(r8), public :: pc_ssor_relax
  integer,  public :: pc_amg_cycles
  integer,  public :: hypre_amg_print_level
  integer,  public :: hypre_amg_debug_level
  integer,  public :: hypre_amg_logging_level

  !! Error tolerances
  real(r8), public :: abs_conc_tol, abs_temp_tol, abs_enthalpy_tol
  real(r8), public :: rel_conc_tol, rel_temp_tol, rel_enthalpy_tol
  real(r8), public :: residual_atol, residual_rtol
  
  real(r8), public :: cond_vfrac_threshold
  
  logical, public :: use_new_mfd  ! use the new MFD mass matrix (temporary)

contains

  subroutine read_ds_namelist (lun)

    use input_utilities, only: seek_to_namelist
    use parallel_communication
    use string_utilities
    
    integer, intent(in) :: lun

    integer :: ios
    logical :: found, hypre_amg_debug
    character(len=32) :: nlk_preconditioner, stepping_method
    character(len=8)  :: string

    !! Magic values used to detect variables not initialized by input
    character, parameter :: NULL_C = char(0)
    integer,   parameter :: NULL_I = huge(1)
    real(r8),  parameter :: NULL_R = huge(1.0_r8)

    namelist /diffusion_solver/ max_step_tries, abs_conc_tol, rel_conc_tol, &
                                abs_temp_tol, rel_temp_tol, abs_enthalpy_tol, rel_enthalpy_tol, &
                                max_nlk_itr, nlk_tol, max_nlk_vec, nlk_vec_tol, &
                                pc_ssor_sweeps, pc_ssor_relax, verbose_stepping, stepping_method, &
                                nlk_preconditioner, pc_amg_cycles, hypre_amg_print_level, &
                                hypre_amg_debug, hypre_amg_logging_level, &
                                cond_vfrac_threshold, residual_atol, residual_rtol, &
                                use_new_mfd

    !! Check that the manually-set module variables have good values before reading.
    INSIST(ds_enabled)
    system_type = raise_case(adjustl(system_type))
    select case (system_type)
    case ('SPECIES')
      ds_sys_type = DS_SPEC_SYS
      heat_eqn = .false.
      INSIST(num_species > 0)
    case ('THERMAL')
      ds_sys_type = DS_TEMP_SYS
      heat_eqn = .true.
      num_species = 0
    case ('THERMAL+SPECIES')
      ds_sys_type = DS_TEMP_SPEC_SYS
      heat_eqn = .true.
      INSIST(num_species > 0)
    case default
      INSIST(.false.)
    end select
    nconc = num_species
    
    call TLS_info ('')
    call TLS_info ('Reading DIFFUSION_SOLVER namelist ...')

    !! Locate the DIFFUSION_SOLVER namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'DIFFUSION_SOLVER', found, iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file')

    call broadcast (found)
    if (.not.found) call TLS_fatal ('DIFFUSION_SOLVER namelist not found')

    !! Read the namelist.
    if (is_IOP) then
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
      pc_ssor_sweeps = NULL_I
      pc_ssor_relax  = NULL_R
      pc_amg_cycles  = NULL_I
      hypre_amg_print_level = NULL_I
      hypre_amg_logging_level = NULL_I
      hypre_amg_debug = .false.
      verbose_stepping = .false.
      stepping_method = NULL_C
      cond_vfrac_threshold = NULL_R
      residual_atol = NULL_R
      residual_rtol = NULL_R
      use_new_mfd = .true.
      read(lun,nml=diffusion_solver,iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading DIFFUSION_SOLVER namelist')

    !! Broadcast the namelist variables.
    call broadcast (max_nlk_itr)
    call broadcast (nlk_tol)
    call broadcast (max_nlk_vec)
    call broadcast (nlk_vec_tol)
    call broadcast (max_step_tries)
    call broadcast (abs_conc_tol)
    call broadcast (rel_conc_tol)
    call broadcast (abs_temp_tol)
    call broadcast (rel_temp_tol)
    call broadcast (abs_enthalpy_tol)
    call broadcast (rel_enthalpy_tol)
    call broadcast (nlk_preconditioner)
    call broadcast (pc_ssor_sweeps)
    call broadcast (pc_ssor_relax)
    call broadcast (pc_amg_cycles)
    call broadcast (hypre_amg_print_level)
    call broadcast (hypre_amg_logging_level)
    call broadcast (hypre_amg_debug)
    call broadcast (verbose_stepping)
    call broadcast (stepping_method)
    call broadcast (cond_vfrac_threshold)
    call broadcast (residual_atol)
    call broadcast (residual_rtol)
    call broadcast (use_new_mfd)
    
    if (stepping_method == NULL_C) then
      stepping_method = 'Adaptive BDF2'
      call TLS_info ('  using default STEPPING_METHOD value: "' // trim(stepping_method) // '"')
    end if
    select case (stepping_method)
    case ('Adaptive BDF2')
      integrator = DS_ADAPTIVE_BDF2
    case ('Non-adaptive BDF1')
      integrator = DS_NONADAPTIVE_BDF1
    case default
      call TLS_fatal ('unknown value for STEPPING_METHOD: "' // trim(stepping_method) // '"')
    end select

    if (max_nlk_itr == NULL_I) then
      max_nlk_itr = 5
      call TLS_info ('  using default MAX_NLK_ITR value: ' // i_to_c(max_nlk_itr))
    else if (max_nlk_itr < 2) then
      call TLS_fatal ('MAX_NLK_ITR must be > 1')
    end if

    select case (integrator)
    case (DS_ADAPTIVE_BDF2)
      if (nlk_tol == NULL_R) then
        nlk_tol = 0.1_r8
        write(string,'(es8.2)') nlk_tol
        call TLS_info ('  using default NLK_TOL value: ' // string)
      else if (nlk_tol <= 0.0_r8 .or. nlk_tol >= 1.0_r8) then
        call TLS_fatal ('NLK_TOL must be positive and < 1')
      end if
      if (max_step_tries == NULL_I) then
        max_step_tries = 10
        call TLS_info ('  using default MAX_STEP_TRIES value: ' // i_to_c(max_step_tries))
      else if (max_step_tries < 1) then
        call TLS_fatal ('MAX_STEP_TRIES must be > 0')
      end if
    case default
      if (nlk_tol /= NULL_R) then
        call TLS_info ('  ignoring NLK_TOL value; not relevant to STEPPING_METHOD choice.')
      end if
      if (max_step_tries /= NULL_I) then
        call TLS_info ('  ignoring MAX_STEP_TRIES value; not relevant to STEPPING_METHOD choice.')
      end if
    end select

    if (max_nlk_vec == NULL_I) then
      max_nlk_vec = max_nlk_itr-1
      call TLS_info ('  using default MAX_NLK_VEC value: MAX_NLK_ITR - 1')
    else if (max_nlk_vec < 1) then
      call TLS_fatal ('MAX_NLK_VEC must be > 0')
    else if (max_nlk_vec >= max_nlk_itr) then
      max_nlk_vec = max_nlk_itr - 1
      call TLS_info ('  reducing MAX_NLK_VEC to MAX_NLK_ITR - 1')
    end if

    if (nlk_vec_tol == NULL_R) then
      nlk_vec_tol = 1.0d-3
      write(string,'(es8.2)') nlk_vec_tol
      call TLS_info ('  using default NLK_VEC_TOL value: ' // string)
    else if (nlk_vec_tol <= 0.0_r8) then
      call TLS_fatal ('NLK_VEC_TOL must be > 0')
    end if

    select case (integrator)
    case (DS_ADAPTIVE_BDF2)
    
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

    nlk_preconditioner = raise_case(adjustl(nlk_preconditioner))
    select case (nlk_preconditioner)
    case ('SSOR')
      ds_nlk_pc = DS_NLK_PC_SSOR
    case ('HYPRE_AMG')
      ds_nlk_pc = DS_NLK_PC_HYPRE_AMG
    case (NULL_C)
      ds_nlk_pc = DS_NLK_PC_SSOR
    case default
      call TLS_fatal ('unknown value for NLK_PRECONDITIONER: ' // trim(nlk_preconditioner))
    end select

    if (ds_nlk_pc == DS_NLK_PC_SSOR) then

      if (pc_ssor_sweeps == NULL_I) then
        pc_ssor_sweeps = 4
        call TLS_info ('  using default PC_SSOR_SWEEPS value: ' // i_to_c(pc_ssor_sweeps))
      else if (pc_ssor_sweeps < 1) then
        call TLS_fatal ('PC_SSOR_SWEEPS must be > 0')
      end if

      if (pc_ssor_relax == NULL_R) then
        pc_ssor_relax = 1.0_r8
        write(string,'(es8.2)') pc_ssor_relax
        call TLS_info ('  using default PC_SSOR_RELAX value: ' // string)
      else if (pc_ssor_relax <= 0.0_r8 .or. pc_ssor_relax >= 2.0_r8) then
        call TLS_fatal ('PC_SSOR_RELAX must be in (0, 2)')
      end if

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

      if (hypre_amg_print_level == NULL_I) then
        hypre_amg_print_level = 0
      else if (hypre_amg_print_level < 0 .or. hypre_amg_print_level > 4) then
        call TLS_fatal ('HYPRE_AMG_PRINT_LEVEL must be >= 0 and <= 3')
      end if

      if (hypre_amg_logging_level == NULL_I) then
        hypre_amg_logging_level = 0
      else if (hypre_amg_logging_level < 0) then
        call TLS_fatal ('HYPRE_AMG_LOGGING_LEVEL must be >= 0')
      end if
      
      hypre_amg_debug_level = 0
      if (hypre_amg_debug) hypre_amg_debug_level = 1

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
      if (residual_atol == NULL_R) then
        residual_atol = 0.0d0
        call TLS_info ('  using default RESIDUAL_ATOL value: 0.0')
      else if (residual_atol < 0.0d0) then
        call TLS_fatal ('RESIDUAL_ATOL must be >= 0')
      end if
      if (residual_rtol == NULL_R) then
        call TLS_fatal ('RESIDUAL_RTOL must be assigned a value')
      else if (residual_rtol <= 0.0d0 .or. residual_rtol >= 1.0d0) then
        call TLS_fatal ('RESIDUAL_RTOL must be > 0 and < 1')
      end if
      if (cond_vfrac_threshold == NULL_R) then
        cond_vfrac_threshold = 1.0d-3
      else if (cond_vfrac_threshold <= 0.0d0 .or. cond_vfrac_threshold >= 1.0d0) then
        call TLS_fatal ('COND_VFRAC_THRESHOLD must be > 0 and < 1')
      end if
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

  end subroutine read_ds_namelist

end module diffusion_solver_data
