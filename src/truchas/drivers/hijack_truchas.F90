!!
!!  HIJACKING TRUCHAS
!!
!!  To hijack Truchas, find the line "!#define HIJACK" in drivers.F90 and
!!  uncomment it.  Recompile.  Now this subroutine will be called *instead* of
!!  CYCLE_DRIVER.  All the usual initialization will have been done and all
!!  of that data is available to this subroutine to use.  When this subroutine
!!  returns, truchas is shutdown as usual.
!!
!!  NB: Do not commit your customized version of this file over the stub file,
!!  nor the modified version of drivers.F90.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

subroutine hijack_truchas

  use truchas_logging_services
  implicit none

  call TLS_info ('*** HIJACKED ***')
  
  ! ADD YOUR CODE HERE
  
end subroutine
