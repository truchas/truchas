!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module FHT_solver_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use FHT_model_type
  use FHT_solver_type
  use matl_mesh_func_type
  implicit none
  private
  
  public :: create_FHT_solver
  
contains

  function create_FHT_solver (mmf, model, stat, errmsg) result (solver)
  
    use enclosure_radiation_namelist, only: er_params => params
    use cutoffs_module, only: cutvof
    use parallel_communication
    use diffusion_solver_data
    use truchas_env, only: output_file_name
    use parameter_list_type
    
    type(matl_mesh_func), intent(in), target :: mmf
    type(FHT_model), intent(in), target :: model
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    type(FHT_solver), pointer :: solver
    
    integer :: j, n, lun
    type(parameter_list) :: params
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    character(:), allocatable :: string
    character(16), allocatable :: vfr_precon_coupling(:)
    real(r8), allocatable :: rad_tol(:)
    
    stat = 0
    errmsg = ''
    
    !! Material volume fraction threshold.
    call params%set('epsilon', cond_vfrac_threshold)

    !! Parameters for the TofH solver.
    !FIXME! we need some user-settable parameters here, especially DELTA
    call params%set('tofh-tol', 0.0d0)
    call params%set('tofh-delta', 1.0d-3)
    call params%set('tofh-max-try', 50)
    
    !! HT system nonlinear solver parameters.
    call params%set('nlk-max-itr', max_nlk_itr)
    call params%set('nlk-max-vec', max_nlk_vec)
    call params%set('nlk-vec-tol', nlk_vec_tol)
    call params%set('verbose', verbose_stepping)
    if (verbose_stepping) then
      lun = -1
      if (is_IOP) open(newunit=lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
      call params%set('unit', lun)
    end if

    !! HT system nonlinear solver preconditioner parameters.
    plist => params%sublist('precon')
    select case (ds_nlk_pc)
    case (DS_NLK_PC_SSOR)
      call plist%set('method', 'SSOR')
      plist => plist%sublist('params')
      call plist%set('num-cycles', pc_ssor_sweeps)
      call plist%set('omega', pc_ssor_relax)
    case (DS_NLK_PC_HYPRE_AMG)
      call plist%set('method', 'BoomerAMG')
      plist => plist%sublist('params')
      call plist%set('num-cycles', pc_amg_cycles)
      call plist%set('print-level', hypre_amg_print_level)
      call plist%set('debug-level', hypre_amg_debug_level)
      call plist%set('logging-level', hypre_amg_logging_level)
    end select      

    !! HT system nonlinear solver error norm parameters.
    plist => params%sublist('norm')
    call plist%set('abs-tol', residual_atol)
    call plist%set('rel-tol', residual_rtol)
    call plist%set('verbose', verbose_stepping)
    if (verbose_stepping) call plist%set('unit', lun) ! defined above
    
    !! Nonlinear solver and preconditioner parameters connected to enclosure radiation.
    if (associated(model%vf_rad_prob)) then
      !TODO! Assumes model%vf_rad_prob array was initialized under the same
      !TODO! loop so that that array and vfr_precon_coupling are in correspondence.
      !TODO! This is fragile and needs to be fixed.
      piter = parameter_list_iterator(er_params, sublists_only=.true.)
      n = piter%count()
      allocate(vfr_precon_coupling(n), rad_tol(n))
      do j = 1, n
        plist => piter%sublist()
        call plist%get('precon-coupling-method', string, default='BACKWARD GS')
        vfr_precon_coupling(j) = string
        call plist%get('error-tol', rad_tol(j), default=1d-3)
        call piter%next
      end do
      plist => params%sublist('precon')
      call plist%set('vfr-precon-coupling', vfr_precon_coupling)
      plist => params%sublist('norm')
      call plist%set('rad-tol', rad_tol)
    end if
    
    allocate(solver)
    call FHT_solver_init (solver, mmf, model, params)
    
  end function create_FHT_solver

end module FHT_solver_factory
