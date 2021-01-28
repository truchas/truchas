!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_nlsol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use nlsol_type
  use sm_model_type
  use sm_ds_precon_type
  !use sm_norm_class
  implicit none
  private

  type, extends(nlsol_model), public :: sm_nlsol_model
    private
    type(sm_model), pointer :: model => null() ! unowned reference
    type(sm_ds_precon), pointer :: precon => null() ! unowned reference
    !type(sm_norm), pointer :: norm => null() ! unowned reference
    integer :: model_size
  contains
    procedure :: init
    !! Deferred procedures from nlsol_model
    procedure :: size => model_size
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
  end type sm_nlsol_model

contains

  subroutine init(this, model, precon)
    class(sm_nlsol_model), intent(out) :: this
    type(sm_model), intent(in), target :: model
    type(sm_ds_precon), intent(in), target :: precon
    !type(sm_norm), intent(in), target :: norm
    this%model => model
    this%precon => precon
    !this%norm => norm
    ASSERT(associated(this%model, precon%model))
    this%model_size = model%size()
  end subroutine init

  integer function model_size(this)
    class(sm_nlsol_model), intent(in) :: this
    model_size = this%model_size
  end function model_size

  subroutine compute_f(this, t, u, udot, f)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:), udot(:)
    real(r8), intent(out), contiguous, target :: f(:)
    real(r8), pointer :: u2(:,:), f2(:,:)
    u2(1:3, 1:this%model%mesh%nnode_onP) => u
    f2(1:3, 1:this%model%mesh%nnode_onP) => f
    call this%model%compute_residual(t, u2, f2)
  end subroutine

  subroutine apply_precon(this, t, u, f)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:)
    real(r8), intent(inout), contiguous, target :: f(:)
    real(r8), pointer :: u2(:,:)
    u2(1:3, 1:this%model_size) => u
    call this%precon%apply(u2, f)
  end subroutine

  subroutine compute_precon(this, t, u, dt)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), contiguous, target :: u(:)
    real(r8), pointer :: u2(:,:)
    u2(1:3, 1:this%model%mesh%nnode_onP) => u
    call this%precon%compute(t, dt, u2)
  end subroutine

  real(r8) function du_norm(this, t, u, du)
    use parallel_communication, only: global_maxval
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:), du(:)
    real(r8) :: l
    du_norm = global_maxval(abs(du))
    l = global_maxval(abs(u))
    if (l > 0) du_norm = du_norm / l
    ! du_norm = 1
    ! l = global_maxval(abs(u))
    ! if (l > 0) du_norm = global_maxval(abs(du)) / l
    ! du_norm = maxval(abs(du) / (this%atol + this%rtol*abs(u)))
    ! du_norm = global_maxval(du_norm)
    ! error = this%norm%compute(u, du)
  end function

end module sm_nlsol_model_type
