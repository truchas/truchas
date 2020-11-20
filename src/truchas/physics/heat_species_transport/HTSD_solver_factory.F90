!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HTSD_solver_factory

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
    type(HTSD_solver_params) :: params
    !character(len=31), allocatable :: encl_name(:)
    type(diff_precon_params) :: precon_params
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    character(:), allocatable :: string
    
    stat = 0
    errmsg = ''

    !! BDF2 control parameters
    params%hmin = hmin
    params%max_step_tries = max_step_tries
    params%verbose_stepping = verbose_stepping
    params%pc_freq = pc_freq
    if (verbose_stepping) then
      lun = -1
      if (is_IOP) open(newunit=lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
      params%output_unit = lun
    end if
    !! BDF2 nonlinear solver parameters
    params%max_nlk_itr = max_nlk_itr
    params%nlk_tol = nlk_tol
    params%max_nlk_vec = max_nlk_vec
    params%nlk_vec_tol = nlk_vec_tol
    !! BDF2 error norm parameters
    if (associated(model%ht)) then
      params%norm_params%abs_T_tol = abs_temp_tol
      params%norm_params%rel_T_tol = rel_temp_tol
      params%norm_params%abs_H_tol = abs_enthalpy_tol
      params%norm_params%rel_H_tol = rel_enthalpy_tol
    end if
    if (associated(model%sd)) then
      allocate(params%norm_params%abs_C_tol(model%num_comp))
      params%norm_params%abs_C_tol = abs_conc_tol
      allocate(params%norm_params%rel_C_tol(model%num_comp))
      params%norm_params%rel_C_tol = rel_conc_tol
    end if
    params%norm_params%verbose = verbose_stepping
    if (verbose_stepping) params%norm_params%unit = lun ! defined above
    !! BDF2 nonlinear solver precon parameters.
    !! We use the same configuration for each of the diffusion equations.
    select case (ds_nlk_pc)
    case (DS_NLK_PC_SSOR)
      precon_params%solver = 'SSOR'
      precon_params%ssor_params%num_iter = pc_ssor_sweeps
      precon_params%ssor_params%omega = pc_ssor_relax
    case (DS_NLK_PC_HYPRE_AMG)
      precon_params%solver = 'BoomerAMG'
      precon_params%bamg_params%max_iter = pc_amg_cycles
      precon_params%bamg_params%print_level = hypre_amg_print_level
      precon_params%bamg_params%debug_level = hypre_amg_debug_level
      precon_params%bamg_params%logging_level = hypre_amg_logging_level
    end select      
    if (associated(model%ht)) then
      params%precon_params%htprecon_params%hcprecon_params = precon_params
      if (associated(model%ht%vf_rad_prob)) then
        !TODO! Assumes model%ht%vf_rad_prob array was initialized under the same
        !TODO! loop so that that array and vfr_precon_coupling are in correspondence.
        !TODO! This is fragile and needs to be fixed.
        piter = parameter_list_iterator(er_params)
        n = piter%count()
        allocate(params%precon_params%htprecon_params%vfr_precon_coupling(n))
        do j = 1, n
          plist => piter%sublist()
          call plist%get('precon-coupling-method', string, default='BACKWARD GS')
          params%precon_params%htprecon_params%vfr_precon_coupling(j) = string
          call piter%next
        end do
      end if
    end if
    if (associated(model%sd)) then
      params%precon_params%sdprecon_params%precon_params = precon_params
    end if
    
    allocate(solver)
    call HTSD_solver_init (solver, mmf, model, params)

  end function create_HTSD_solver
  
end module HTSD_solver_factory
