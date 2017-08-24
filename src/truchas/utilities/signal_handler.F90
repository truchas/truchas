!!
!! SIGNAL_HANDLER
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! CALL INIT_SIGNAL_HANDLER(SIGNUM) will register a signal handler for the
!! signal SIGNUM that sets an associated global logical flag to true whenever
!! the signal is received. Use READ_SIGNAL to get the value of that flag.
!! SIGNUM should be one of the parameters defined by this module:
!! SIGHUP, SIGUSR1, SIGUSR2, or SIGURG.
!!
!! CALL READ_SIGNAL(SIGNUM, SIG_RCVD) returns the status of the SIGNUM flag
!! in the logical argument SIG_RCVD, and sets the flag to false.
!!
!! NB: This uses the old deprecated ANSI signal function. It works under Linux,
!! but may not be portable. Its behavior in a threaded environment is uncertain.
!! See 'man 2 signal' for all the cautions.
!!

module signal_handler

  use,intrinsic :: iso_c_binding, only: c_int
  implicit none
  private

  public :: init_signal_handler, read_signal

  !! From /usr/include/bits/signum.h on Linux
  integer(c_int), parameter, public :: SIGHUP  = 1
  integer(c_int), parameter, public :: SIGUSR1 = 10
  integer(c_int), parameter, public :: SIGUSR2 = 12
  integer(c_int), parameter, public :: SIGURG  = 23

  abstract interface
    subroutine sighandler_t(signum) bind(c)
      import c_int
      integer(c_int), value :: signum
    end subroutine
  end interface

  interface
    subroutine signal(signum, sighandler) bind(c,name='signal')
      import c_int, sighandler_t
      integer(c_int), value :: signum
      procedure(sighandler_t) :: sighandler
    end subroutine
  end interface

  logical, save :: hup_rcvd  = .false.
  logical, save :: usr1_rcvd = .false.
  logical, save :: usr2_rcvd = .false.
  logical, save :: urg_rcvd  = .false.

contains

  subroutine init_signal_handler(signum)
    integer(c_int), intent(in) :: signum
    call signal(signum, sighandler)
  end subroutine init_signal_handler

  subroutine sighandler(signum) bind(c,name='')
    integer(c_int), value :: signum
    select case (signum)
    case (SIGHUP)
      hup_rcvd = .true.
    case (SIGUSR1)
      usr1_rcvd = .true.
    case (SIGUSR2)
      usr2_rcvd = .true.
    case (SIGURG)
      urg_rcvd = .true.
    end select
  end subroutine sighandler

  subroutine read_signal(signum, sig_rcvd)
    integer(c_int), intent(in) :: signum
    logical, intent(out) :: sig_rcvd
    select case (signum)
    case (SIGHUP)
      sig_rcvd = hup_rcvd
      hup_rcvd = .false.
    case (SIGUSR1)
      sig_rcvd = usr1_rcvd
      usr1_rcvd = .false.
    case (SIGUSR2)
      sig_rcvd = usr2_rcvd
      usr2_rcvd = .false.
    case (SIGURG)
      sig_rcvd = urg_rcvd
      urg_rcvd = .false.
    case default
      sig_rcvd = .false.
    end select
  end subroutine read_signal

end module signal_handler
