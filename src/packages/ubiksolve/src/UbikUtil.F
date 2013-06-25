! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: UbikUtil.F,v 1.2 2002/05/09 05:40:20 turner Exp $
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
module UbikUtil_module
  use UbikBase_module

  ! Input processing routines.
  use Ubik_ParseArg_module, only: &
       Ubik_ParseArg => ParseArg

  ! Output routines.
  use Ubik_FlushUnit_module, only: &
       Ubik_FlushUnit => FlushUnit
  use Ubik_WriteString_module, only: &
       Ubik_Spew => Spew

  ! Utilities.
  use Ubik_ErrorFunction_module, only: &
       Ubik_ErrorFunction => ErrorFunction

end module UbikUtil_module
