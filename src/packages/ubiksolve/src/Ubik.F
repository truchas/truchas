! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: Ubik.F,v 1.1 2002/03/22 03:59:13 turner Exp $
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
module Ubik_module

  use UbikMatrix_module
  use UbikSolve_module

end module Ubik_module
