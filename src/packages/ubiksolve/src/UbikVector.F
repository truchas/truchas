! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: UbikVector.F,v 1.8 2002/07/01 06:19:11 turner Exp $
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
module UbikVector_module
  use UbikUtil_module

  use ubik_vector_type_module, only: &
       Ubik_vector_type => vector_type, &
       Ubik_create => create, &
       Ubik_destroy => destroy, &
       assignment(=), &
       Ubik_set_values_ptr => set_values_ptr, &
       Ubik_set_aux1_ptr => set_aux1_ptr, &
       Ubik_values => values, &
       Ubik_values_ptr => values_ptr, &
       Ubik_aux1_ptr => aux1_ptr, &
       Ubik_aux1_elem => aux1_elem, &
       Ubik_spew => spew, &
       SIZE, &
       Ubik_Norm => Norm, &
       Ubik_set_status_alloc_failure => set_status_alloc_failure, &
       Ubik_alloc_failure => alloc_failure, &
#ifdef dec
       Ubik_vector_minus_vector_type => vector_minus_vector_type, &
       Ubik_vector_type_minus_vector => vector_type_minus_vector
#else
       operator(-)
#endif

  use Ubik_Vector_module, only: &
       Ubik_Dot_Product, &
       Ubik_Dump => Dump, &
       Ubik_Load => Load, &
       Ubik_MaxVal, &
       Ubik_Norm => Norm, &
       Ubik_Spew => Spew, &
       Ubik_Sum, &
       Ubik_Swap => Swap

end module UbikVector_module
