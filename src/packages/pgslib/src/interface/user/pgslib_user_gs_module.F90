MODULE PGSLib_User_GS_MODULE
  !======================================================================
  ! PURPOSE
  !   User level gather and scatter routines.  This is the mother module.
  !
  !======================================================================

  ! $Id: pgslib_user_gs_module.F,v 1.1.1.1 2000/10/11 22:44:28 ferrell Exp $

  USE PGSLib_Gather_MODULE,  ONLY : PGSLib_Gather
  USE PGSLib_Scatter_MODULE, ONLY : PGSLib_Scatter_SUM,  &
       &                               PGSLib_Scatter_MAX,  &
       &                               PGSLib_Scatter_MIN,  &
       &                               PGSLib_Scatter_AND,  &
       &                               PGSLib_Scatter_OR

END MODULE PGSLib_User_GS_MODULE
