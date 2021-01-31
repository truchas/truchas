!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HTSD_solver_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use HTSD_model_type
  use HTSD_solver_type
  use matl_mesh_func_type
  implicit none
  private
  
  public :: create_HTSD_solver
  
contains

  function create_HTSD_solver (mmf, model, stat, errmsg) result (solver)
  
    use enclosure_radiation_namelist, only: er_params => params
    use parallel_communication
    use diffusion_solver_data
    use truchas_env, only: output_file_name
    use parameter_list_type

    type(matl_mesh_func), intent(in), target :: mmf
    type(HTSD_model), intent(in), target :: model
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    type(HTSD_solver), pointer :: solver
    
    integer :: j, n, lun
    type(parameter_list) :: params
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    character(:), allocatable :: string
    character(16), allocatable :: vfr_precon_coupling(:)
    real(r8), allocatable :: rad_tol(:)
    
    stat = 0
    errmsg = ''

    !! BDF2 control parameters
    call params%set('hmin', hmin)
    call params%set('max-step-tries', max_step_tries)
    call params%set('verbose-stepping', verbose_stepping)
    call params%set('pc-freq', pc_freq)
    if (verbose_stepping) then
      lun = -1
      if (is_IOP) open(newunit=lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
      call params%set('output-unit', lun)
    end if
    !! BDF2 nonlinear solver parameters
    call params%set('nlk-max-iter', max_nlk_itr)
    call params%set('nlk-tol', nlk_tol)
    call params%set('nlk-max-vec', max_nlk_vec)
    call params%set('nlk-vec-tol', nlk_vec_tol)
    !! BDF2 error norm parameters
    plist => params%sublist('norm')
    if (associated(model%ht)) then
      call plist%set('abs-t-tol', abs_temp_tol)
      call plist%set('rel-t-tol', rel_temp_tol)
      call plist%set('abs-h-tol', abs_enthalpy_tol)
      call plist%set('rel-h-tol', rel_enthalpy_tol)
    end if
    if (associated(model%sd)) then
      call plist%set('abs-c-tol', spread(abs_conc_tol, dim=1, ncopies=model%num_comp))
      call plist%set('rel-c-tol', spread(rel_conc_tol, dim=1, ncopies=model%num_comp))
    end if
    call plist%set('verbose', verbose_stepping)
    if (verbose_stepping) call plist%set('unit', lun) ! defined above
    !! BDF2 nonlinear solver precon parameters.
    !! We use the same configuration for each of the diffusion equations.
    plist => params%sublist('precon')
    select case (ds_nlk_pc)
    case (DS_NLK_PC_SSOR)
      call plist%set('method', 'SSOR')
      plist => plist%sublist('params')
      call plist%set('num-sweeps', pc_ssor_sweeps)
      call plist%set('omega', pc_ssor_relax)
    case (DS_NLK_PC_HYPRE_AMG)
      call plist%set('method', 'BoomerAMG')
      plist => plist%sublist('params')
      call plist%set('num-cycles', pc_amg_cycles)
      call plist%set('print-level', hypre_amg_print_level)
      call plist%set('debug-level', hypre_amg_debug_level)
      call plist%set('logging-level', hypre_amg_logging_level)
    end select      
    if (associated(model%ht)) then
      if (associated(model%ht%vf_rad_prob)) then
        !TODO! Assumes model%ht%vf_rad_prob array was initialized under the same
        !TODO! loop so that that array and vfr_precon_coupling are in correspondence.
        !TODO! This is fragile and needs to be fixed.
        piter = parameter_list_iterator(er_params)
        n = piter%count()
        allocate(vfr_precon_coupling(n), rad_tol(n))
        do j = 1, n
          plist => piter%sublist()
          call plist%get('precon-coupling-method', string, default='BACKWARD GS')
          vfr_precon_coupling(j) = string
          call piter%next
        end do
        plist => params%sublist('precon')
        call plist%set('vfr-precon-coupling', vfr_precon_coupling)
      end if
    end if
    
    allocate(solver)
    call HTSD_solver_init (solver, mmf, model, params)

  end function create_HTSD_solver
  
end module HTSD_solver_factory
