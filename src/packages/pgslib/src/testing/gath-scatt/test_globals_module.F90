MODULE TEST_Globals_MODULE
  USE PGSLib_MODULE
  IMPLICIT NONE
  SAVE

  ! Contains the global variables for the users interface to PGSLib.  These
  ! variables are fully in the control of the user, not under the control of PGSLib.

  ! $Id: test_globals_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $


  ! The trace for the Element <-> Node communication

  Type (PGSLib_GS_Trace), POINTER :: EN_Trace

END MODULE TEST_Globals_MODULE
