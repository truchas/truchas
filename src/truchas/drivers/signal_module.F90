!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module signal_module
  !-----------------------------------------------------------------------------
  ! purpose
  !
  ! provide explicit interfaces for signal handling code
  !
  ! note: the integers used as arguments here are _not_ of type
  ! INT_KIND.  The are standard default integers.  They
  ! _must_ agree with what is in signal.C.
  !-----------------------------------------------------------------------------

  Implicit None
  Private

  Public SignalSet, SignalInquire

  Interface
     Subroutine SignalSet() bind(c, name='signalset')
     End Subroutine SignalSet
  End Interface

  Interface
     Subroutine SignalInquire(hup, usr2, urg) bind(c, name='signalinquire')
       use,intrinsic :: iso_c_binding, only: c_int
       Integer(c_int), intent(out) :: hup, usr2, urg
     End Subroutine SignalInquire
  End Interface

End Module signal_module
