!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module FHT_solver_factory

  use FHT_model_type
  use FHT_solver_type
  use material_mesh_function
  implicit none
  private
  
  public :: create_FHT_solver
  
contains

  function create_FHT_solver (mmf, model, stat, errmsg) result (solver)
  
    use ER_input
    use cutoffs_module, only: cutvof
    use parallel_communication
    use diffusion_solver_data
#ifdef SUPPORTS_NEWUNIT
    use truchas_env, only: output_file_name
#else
    use truchas_env, only: output_file_name, new_unit
#endif
    
    type(mat_mf), intent(in), target :: mmf
    type(FHT_model), intent(in), target :: model
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    type(FHT_solver), pointer :: solver
    
    integer :: j, n, lun
    type(FHT_solver_params) :: params
    character(len=31), allocatable :: encl_name(:)
    
    stat = 0
    errmsg = ''
    
    !! Material volume fraction threshold.
    !FIXME! we need a user-settable parameter for epsilon
    params%epsilon = cond_vfrac_threshold  ! conduction threshold

    !! Parameters for the TofH solver.
    !FIXME! we need some user-settable parameters here, especially DELTA
    params%TofH_tol = 0.0d0       ! absolute temperature convergence tolerance, >= 0
    params%TofH_delta  = 1.0d-3   ! initial endpoint shift when seeking bracket, > 0
    params%TofH_max_try = 50      ! max tries at seeking a bracketing interval, >= 0
    
    !! HT system nonlinear solver parameters.
    params%nlk_max_itr = max_nlk_itr
    params%nlk_max_vec = max_nlk_vec
    params%nlk_vec_tol = nlk_vec_tol
    params%verbose = verbose_stepping
    if (verbose_stepping) then
      lun = -1
#ifdef SUPPORTS_NEWUNIT
      if (is_IOP) open(newunit=lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
#else
      if (is_IOP) then
        call new_unit (lun)
        open(lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
      end if
#endif
      params%unit = lun
    end if

    !! HT system nonlinear solver preconditioner parameters.
    select case (ds_nlk_pc)
    case (DS_NLK_PC_SSOR)
      params%precon_params%HC_precon_params%solver = 'SSOR'
      params%precon_params%HC_precon_params%ssor_params%num_iter = pc_ssor_sweeps
      params%precon_params%HC_precon_params%ssor_params%omega = pc_ssor_relax
    case (DS_NLK_PC_HYPRE_AMG)
      params%precon_params%HC_precon_params%solver = 'BoomerAMG'
      params%precon_params%HC_precon_params%bamg_params%max_iter = pc_amg_cycles
      params%precon_params%HC_precon_params%bamg_params%print_level = hypre_amg_print_level
      params%precon_params%HC_precon_params%bamg_params%debug_level = hypre_amg_debug_level
      params%precon_params%HC_precon_params%bamg_params%logging_level = hypre_amg_logging_level
    end select      

    !! HT system nonlinear solver error norm parameters.
    params%norm_params%abs_tol = residual_atol
    params%norm_params%rel_tol = residual_rtol
    params%norm_params%verbose = verbose_stepping
    if (verbose_stepping) params%norm_params%unit = lun	! defined above
    
    !! Nonlinear solver and preconditioner parameters connected to enclosure radiation.
    if (associated(model%vf_rad_prob)) then
      !TODO! Assumes model%vf_rad_prob array was initialized under the same
      !TODO! loop so that that array and vfr_precon_coupling are in correspondence.
      !TODO! This is fragile and needs to be fixed.
      n = ERI_num_enclosures()
      allocate(encl_name(n))
      allocate(params%precon_params%vfr_precon_coupling(n), params%norm_params%rad_tol(n))
      call ERI_get_names (encl_name)
      do j = 1, n
        call ERI_get_precon_coupling_method (encl_name(j), params%precon_params%vfr_precon_coupling(j))
        call ERI_get_error_tolerance (encl_name(j), params%norm_params%rad_tol(j))
      end do
      deallocate(encl_name)
    end if
    
    allocate(solver)
    call FHT_solver_init (solver, mmf, model, params)
    
  end function create_FHT_solver

end module FHT_solver_factory
