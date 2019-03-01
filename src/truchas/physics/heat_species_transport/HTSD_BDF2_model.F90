!!
!! HTSD_BDF2_MODEL
!!
!! This module provides the call-back procedures passed to the BDF2 DAE solver
!! that implement the heat transfer/species diffusion system.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The public procedures BDF2_PCFUN, BDF2_UPDPC, BDF2_ENORM implement the
!!  call-back procedures required by the BDF2 DAE integrator.  These require
!!  look-aside data (or a context) to perform their work.  This context data
!!  is stored as a private data structure variable in this module and it must
!!  be defined before invoking any procedure these call-backs are passed to.
!!  CALL SET_CONTEXT (MODEL, PRECON, NORM) sets this context data:
!!    TYPE(HTSD_MODEL),  TARGET :: MODEL
!!    TYPE(HTSD_PRECON), TARGET :: PRECON
!!    TYPE(HTSD_NORM),   TARGET :: NORM
!!  The actual arguments must either be pointers or have the target attribute.
!!
!! IMPLEMENTATION NOTES
!!
!!  Transitioning to F2003.  The design is intended to make the change to an
!!  OO implementation straightforward:
!!  * The CONTEXT type becomes a public HTSD_BDF2_MODEL type with the call-backs
!!    as type-bound procedures.  This type should extend an abstract base type
!!    that is defined by the BDF2 DAE solver, which should take arguments of
!!    this type in place of the bare call-back procedures.
!!  * The call-back interfaces are modified to have THIS as the first argument.
!!  * Because instances of the new type carry the context, the private context
!!    module data is no longer needed.
!!

#include "f90_assert.fpp"

module HTSD_BDF2_model

  use kinds, only: r8
  use HTSD_model_type
  use HTSD_precon_type
  use HTSD_norm_type
  implicit none
  private

  public :: set_context
  public :: bdf2_pcfun, bdf2_updpc, bdf2_enorm!, bdf2_enorm2
  !public :: bdf1_fun, bdf1_pc, bdf1_updpc, bdf1_fnorm

  type :: context
    type(HTSD_model),  pointer :: model  => null()
    type(HTSD_precon), pointer :: precon => null()
    type(HTSD_norm),   pointer :: norm   => null()
  end type context
  type(context), save :: this

contains

  subroutine set_context (model, precon, norm)
    type(HTSD_model),  target :: model
    type(HTSD_precon), target :: precon
    type(HTSD_norm),   target :: norm
    this%model  => model
    this%precon => precon
    this%norm   => norm
  end subroutine set_context

  subroutine bdf2_pcfun (t, u, udot, f)
    real(r8), intent(in) :: t, u(:), udot(:)
    real(r8), intent(out) :: f(:)
    call HTSD_model_compute_f (this%model, t, u, udot, f)
    call HTSD_precon_apply (this%precon, t, u, f)
  end subroutine bdf2_pcfun

  subroutine bdf2_updpc (t, u, h, errc)
    real(r8), intent(in) :: t, u(:), h
    integer, intent(out) :: errc
    call HTSD_precon_compute (this%precon, t, u, h, errc)
  end subroutine bdf2_updpc

  function bdf2_enorm (u, du) result (enorm)
    real(r8), intent(in) :: u(:), du(:)
    real(r8) :: enorm
    call HTSD_norm_compute (this%norm, u, du, enorm)
  end function bdf2_enorm

!  function bdf2_enorm2 (u, du) result (enorm)
!    real(r8), intent(in) :: u(:), du(:)
!    real(r8) :: enorm
!    call HT_norm_compute2 (this%norm, u, du, enorm)
!  end function bdf2_enorm2
!
!  subroutine bdf1_fun (t, u, udot, f)
!    real(r8), intent(in) :: t, u(:), udot(:)
!    real(r8), intent(out) :: f(:)
!    call HT_model_compute_f (this%model, t, u, udot, f)
!  end subroutine bdf1_fun
!
!  subroutine bdf1_pc (t, u, udot, f)
!    real(r8), intent(in) :: t, u(:), udot(:)
!    real(r8), intent(inout) :: f(:)
!    call HT_precon_apply (this%precon, t, u, f)
!  end subroutine bdf1_pc
!
!  subroutine bdf1_updpc (t, u, h, errc)
!    real(r8), intent(in) :: t, u(:), h
!    integer, intent(out) :: errc
!    call HT_precon_compute (this%precon, t, u, h, errc)
!  end subroutine bdf1_updpc
!
!  subroutine bdf1_fnorm (t, u, udot, f, error)
!    real(r8), intent(in) :: t, u(:), udot(:), f(:)
!    real(r8), intent(out), optional :: error
!    call HT_norm_fnorm (this%norm, t, u, udot, f, error)
!  end subroutine bdf1_fnorm

end module HTSD_BDF2_model
