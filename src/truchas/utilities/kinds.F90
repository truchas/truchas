!!
!! KINDS
!!
!! This module provides kind parameter values for integer and real types
!! that are identified by byte sizes.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module kinds

  public

  integer, parameter :: r4 = selected_real_kind(6)   ! 4-byte IEEE float
  integer, parameter :: r8 = selected_real_kind(15)  ! 8-byte IEEE float

  integer, parameter :: i4 = selected_int_kind(9)    ! 4-byte IEEE integer
  integer, parameter :: i8 = selected_int_kind(18)   ! 8-byte IEEE integer

  !! The ISO_C_BINDING intrinsic module defines the integer kind parameter
  !! C_INTPTR_T that is equivalent to the C typedef intptr_t.  We don't yet
  !! want to require that feature so for know we define it ourselves.  If
  !! the macro M64 is defined its an 8-byte integer, otherwise a 4-byte
  !! integer.  This is fragile.
#ifdef M64
  integer, parameter :: c_intptr_t = i8
#else
  integer, parameter :: c_intptr_t = i4
#endif

end module kinds
