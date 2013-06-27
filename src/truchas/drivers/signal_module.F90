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
     Subroutine SignalSet ()
     End Subroutine SignalSet
  End Interface

  Interface
     Subroutine SignalInquire (hup, usr2, urg)
       Integer :: hup
       Integer :: usr2
       Integer :: urg
     End Subroutine SignalInquire
  End Interface

End Module signal_module
