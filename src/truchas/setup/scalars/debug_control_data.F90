!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DEBUG_CONTROL_DATA
  !=======================================================================
  ! Purpose:
  !
  !    Define data which are used to turn on debugging and verbose
  !    capabilities.
  !
  ! Author(s): Robert C. Ferrell (ferrell@cpca.com)
  !=======================================================================

  implicit none
  private

  public :: verbose, VERBOSE_QUIET, VERBOSE_NORMAL, VERBOSE_NOISY, VERBOSE_DEFAULT, VERBOSE_DEFAULT_SET
  public :: debug, DEBUG_NONE, DEBUG_QUIET, DEBUG_NOISY, DEBUG_DEFAULT, DEBUG_DEFAULT_SET

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  ! these two sets of constants should use the same terms for their levels [lally]

  ! parameters which define levels of verbosity
  Integer, Parameter :: VERBOSE_QUIET       = 0
  Integer, Parameter :: VERBOSE_NORMAL      = 1
  Integer, Parameter :: VERBOSE_NOISY       = 2
  ! if set on command line, with no level specified
  Integer, Parameter :: VERBOSE_DEFAULT_SET = VERBOSE_NOISY
  ! if not set on command line (program default)
  Integer, Parameter :: VERBOSE_DEFAULT     = VERBOSE_NORMAL
  
  ! parameters which define levels of debugging
  Integer, Parameter :: DEBUG_NONE          = 0 
  Integer, Parameter :: DEBUG_QUIET         = 1
  Integer, Parameter :: DEBUG_NOISY         = 2
  ! if set on command line, with no level specified
  Integer, Parameter :: DEBUG_DEFAULT_SET   = DEBUG_QUIET
  ! if not set on command line (program default)
  Integer, Parameter :: DEBUG_DEFAULT       = DEBUG_NONE

  ! global flags for controlling verbosity of output
  Integer, Save :: debug
  Integer, Save :: verbose

END MODULE DEBUG_CONTROL_DATA
