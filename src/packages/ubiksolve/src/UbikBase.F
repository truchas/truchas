! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: UbikBase.F,v 1.2 2002/03/30 16:48:40 turner Exp $
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
module UbikBase_module

  use ubik_constants_module, only: &
       Ubik_no => no, &
       Ubik_yes => yes, &
       Ubik_false => false, &
       Ubik_true => true, &
       Ubik_success__ => success__, &
       Ubik_converged__ => converged__, &
       Ubik_invalid_input__ => invalid_input__, &
       Ubik_alloc_failure__ => alloc_failure__, &
       Ubik_internal_error__ => internal_error__, &
       Ubik_itmax_exceeded__ => itmax_exceeded__, &
       Ubik_breakdown__ => breakdown__, &
       Ubik_matmul_error__ => matmul_error__, &
       Ubik_precond_error__ => precond_error__, &
       Ubik_next_iter__ => next_iter__

  use ubik_types_module, only: &
       Ubik_logical_type => logical_type, &
       Ubik_int_type => int_type, &
       Ubik_real_type => real_type

  use Ubik_ParallelInfo_module, only: &
       Ubik_on_IO_PE => on_IO_PE, &
       Ubik_initialize_parallel => initialize_parallel, &
       Ubik_initialize_scope => initialize_scope

end module UbikBase_module
