! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: UbikSolve.F,v 1.18 2002/12/12 01:03:40 turner Exp $
!
! Primary interface to UbikSolve.  Host codes should access UbikSolve data types,
! routines, etc. only through this module (or one like it).
!
! The prefix Ubik_ is prepended to all publicly exposed entities in order to
! help guard against namespace collision.  If it should be a problem in
! spite of this, alternative interfaces to UbikSolve can easily be created by
! copying this module and replacing Ubik_ with some other Fortran-legal string.
!
! The "only" form of "use" is used, so all entities are listed, even if not
! renamed (some have Ubik_ prepended already).  This makes 
!
! Note that we need not worry about public/private here, since that is
! controlled within the modules comprising UbikSolve - everything they expose
! externally can be exposed to host codes.
!
!*********************************************************************
module UbikSolve_module
  use UbikVector_module

  use ubik_constants_module, only: &
       Ubik_method_min => method_min, &
       Ubik_precond_none => precond_none, &
       Ubik_method_LU => method_LU, &
       Ubik_method_CG => method_CG, &
       Ubik_method_GMRES => method_GMRES, &
       Ubik_method_GMRES_Left_Radicati => method_GMRES_Left_Radicati, &
       Ubik_method_GMRES_Left_Saad => method_GMRES_Left_Saad, &
       Ubik_method_FGMRES => method_FGMRES, &
       Ubik_method_BCGS => method_BCGS, &
       Ubik_method_TFQMR => method_TFQMR, &
       Ubik_method_BCGSTAB => method_BCGSTAB, &
       Ubik_method_SOR => method_SOR, &
       Ubik_method_SSOR => method_SSOR, &
       Ubik_method_Jacobi => method_Jacobi, &
       Ubik_method_IC => method_IC, &
       Ubik_method_ILU => method_ILU, &
       Ubik_method_ApproxInv => method_ApproxInv, &
       Ubik_method_max => method_max

  ! derived type containing control params
  use ubik_control_type_module, only: &
       Ubik_control_type => control_type, &
#ifdef dec
       Ubik_set_control => set_control, &
#else
       Ubik_set => set, &
#endif
       Ubik_set_defaults => set_defaults, &
       Ubik_initialize => initialize, &
       Ubik_create => create, &
       Ubik_nullify, &
       Ubik_destroy => destroy

  ! routines to set components
  use ubik_control_type_module, only: &
       Ubik_set_iter => set_iter, &
       Ubik_set_itmax => set_itmax, &
       Ubik_set_luout => set_luout, &
       Ubik_set_luerr => set_luerr, &
       Ubik_set_lutty => set_lutty, &
       Ubik_set_method => set_method, &
       Ubik_set_subspace => set_subspace, &
       Ubik_set_nold => set_nold, &
       Ubik_set_norm_type => set_norm_type, &
       Ubik_set_outmode_none => set_outmode_none, &
       Ubik_set_outmode_errors => set_outmode_errors, &
       Ubik_set_outmode_warnings => set_outmode_warnings, &
       Ubik_set_outmode_summary => set_outmode_summary, &
       Ubik_set_outmode_system => set_outmode_system, &
       Ubik_set_outmode_iterates => set_outmode_iterates, &
       Ubik_set_outmode_full => set_outmode_full, &
       Ubik_set_output_frequency => set_output_frequency

  use ubik_control_type_module, only: &
       Ubik_set_residual => set_residual, &
         Ubik_set_residual_default => set_residual_default, &
         Ubik_set_residual_compute_only => set_residual_compute_only, &
         Ubik_set_residual_update => set_residual_update, &
         Ubik_set_residual_compute_upd => set_residual_compute_upd, &
       Ubik_set_scope => set_scope, &
         Ubik_set_scope_global => set_scope_global, &
         Ubik_set_scope_local => set_scope_local, &
       Ubik_set_steps => set_steps, &
       Ubik_set_stopping_relchg => set_stopping_relchg, &
       Ubik_set_stopping_axb => set_stopping_axb, &
       Ubik_set_stopping_b => set_stopping_b, &
       Ubik_set_stopping_x => set_stopping_x, &
       Ubik_set_stopping_r0 => set_stopping_r0, &
       Ubik_set_stopping_r => set_stopping_r, &
       Ubik_set_stopping_true => set_stopping_true, &
       Ubik_set_stopping_recursive => set_stopping_recursive

  use ubik_control_type_module, only: &
       Ubik_set_anorm => set_anorm, &
       Ubik_set_bnorm => set_bnorm, &
       Ubik_set_eps => set_eps, &
       Ubik_set_omega => set_omega, &
       Ubik_set_small => set_small, &
       Ubik_set_smallest => set_smallest, &
       Ubik_set_cpu => set_cpu, &
       Ubik_set_rnormt => set_rnormt, &
       Ubik_set_errt => set_errt, &
       Ubik_set_rnorm => set_rnorm, &
       Ubik_set_err => set_err

  ! routines to access params
  use ubik_control_type_module, only: &
       Ubik_freq => freq, &
       Ubik_iter => iter, &
       Ubik_itmax => itmax, &
       Ubik_luout => luout, &
       Ubik_luerr => luerr, &
       Ubik_matvecs => matvecs, &
       Ubik_method => method, &
       Ubik_subspace => subspace, &
       Ubik_nold => nold, &
       Ubik_norm_type => norm_type, &
       Ubik_outmode_none => outmode_none, &
       Ubik_outmode_errors => outmode_errors, &
       Ubik_outmode_warnings => outmode_warnings, &
       Ubik_outmode_summary => outmode_summary, &
       Ubik_outmode_system => outmode_system, &
       Ubik_outmode_iterates => outmode_iterates, &
       Ubik_outmode_full => outmode_full, &
       Ubik_output_frequency => output_frequency

  use ubik_control_type_module, only: &
       Ubik_preconds => preconds, &
       Ubik_residual => residual, &
         Ubik_residual_default => residual_default, &
         Ubik_residual_compute_only => residual_compute_only, &
         Ubik_residual_update => residual_update, &
         Ubik_residual_compute_update => residual_compute_update, &
       Ubik_scope => scope, &
         Ubik_scope_global => scope_global, &
         Ubik_scope_local => scope_local, &
       Ubik_steps => steps, &
       Ubik_stopping_true => stopping_true, &
       Ubik_stopping_relchg => stopping_relchg, &
       Ubik_stopping_axb => stopping_axb, &
       Ubik_stopping_b => stopping_b, &
       Ubik_stopping_x => stopping_x, &
       Ubik_stopping_r0 => stopping_r0, &
       Ubik_stopping_r => stopping_r

  use ubik_control_type_module, only: &
       Ubik_anorm => anorm, &
       Ubik_bnorm => bnorm, &
       Ubik_eps => eps, &
       Ubik_omega => omega, &
       Ubik_small => small, &
       Ubik_smallest => smallest, &
       Ubik_cpu => cpu, &
       Ubik_rnormt => rnormt, &
       Ubik_errt => errt, &
       Ubik_rnorm => rnorm, &
       Ubik_err => err

  use ubik_control_type_module, only: &
       Ubik_converged => converged, &
       Ubik_invalid_input => invalid_input, &
       Ubik_alloc_failure => alloc_failure, &
       Ubik_internal_error => internal_error, &
       Ubik_itmax_exceeded => itmax_exceeded, &
       Ubik_breakdown => breakdown, &
       Ubik_matmul_error => matmul_error, &
       Ubik_precond_error => precond_error, &
       Ubik_next_iter => next_iter, &
       Ubik_set_status_converged => set_status_converged, &
       Ubik_set_status_invalid_input => set_status_invalid_input, &
       Ubik_set_status_alloc_failure => set_status_alloc_failure, &
       Ubik_set_status_internal_error => set_status_internal_error, &
       Ubik_set_status_itmax_exceeded => set_status_itmax_exceeded, &
       Ubik_set_status_breakdown => set_status_breakdown, &
       Ubik_set_status_matmul_error => set_status_matmul_error, &
       Ubik_set_status_precond_error => set_status_precond_error, &
       Ubik_set_status_next_iter => set_status_next_iter

  ! Utilities.

  use Ubik_Convergence_module, only: &
       Ubik_StoppingTest => StoppingTest, &
       Ubik_CheckConvergence => CheckConvergence

  ! Iterative solver routines.

  use Ubik_BCGS_module, only: &
       Ubik_BCGS => BCGS

  use Ubik_BCGSTAB_module, only: &
       Ubik_BCGSTAB => BCGSTAB

  use Ubik_CG_module, only: &
       Ubik_CG => CG

  use Ubik_GMRES_module, only: &
       Ubik_GMRES => GMRES, &
       Ubik_GMRES_Left_Radicati => GMRES_Left_Radicati, &
       Ubik_GMRES_Left_Saad => GMRES_Left_Saad, &
       Ubik_FGMRES => FGMRES

  use Ubik_Stationary_module, only: &
       Ubik_Stationary => Stationary

  use Ubik_TFQMR_module, only: &
       Ubik_TFQMR => TFQMR

end module UbikSolve_module
