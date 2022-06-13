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

    real(r8) :: atol, rtol, ftol
  contains
    procedure :: init
    !! Deferred procedures from nlsol_model
    procedure :: size => model_size
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
    procedure :: is_converged
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
    this%atol = this%model%atol
    this%rtol = this%model%rtol
    this%ftol = this%model%ftol
  end subroutine init

  integer function model_size(this)
    class(sm_nlsol_model), intent(in) :: this
    model_size = 3*this%model%mesh%nnode
  end function model_size

  subroutine compute_f(this, t, u, udot, f)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:), udot(:)
    real(r8), intent(out), contiguous, target :: f(:)
    real(r8), pointer :: u2(:,:), f2(:,:)
    u2(1:3, 1:this%model%mesh%nnode) => u
    f2(1:3, 1:this%model%mesh%nnode) => f
    call this%model%compute_residual(t, u2, f2)
    f2(:,this%model%mesh%nnode_onP+1:) = 0 ! clear out halo so it doesn't screw with norms
  end subroutine

  subroutine apply_precon(this, t, u, f)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:)
    real(r8), intent(inout), contiguous, target :: f(:)
    real(r8), pointer :: u2(:,:), f2(:,:)
    u2(1:3, 1:this%model%mesh%nnode) => u
    f2(1:3, 1:this%model%mesh%nnode) => f
    call this%precon%apply(u2, f2)
    f2(:,this%model%mesh%nnode_onP+1:) = 0 ! clear out halo so it doesn't screw with norms
  end subroutine

  subroutine compute_precon(this, t, u, dt)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), contiguous, target :: u(:)
    real(r8), pointer :: u2(:,:)
    u2(1:3, 1:this%model%mesh%nnode) => u
    call this%precon%compute(t, dt, u2)
  end subroutine

  real(r8) function du_norm(this, t, u, du)

    use parallel_communication, only: global_maxval
    use,intrinsic :: ieee_arithmetic, only: ieee_is_finite

    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:), du(:)
    integer :: i
    real(r8) :: err

    du_norm = 0
    do i = 1, size(u)
      if ((this%atol == 0 .and. abs(u(i)) == 0)) then
        err = huge(1.0_r8)
      else
        err = abs(du(i)) / (this%atol + this%rtol*abs(u(i)))
      end if
      du_norm = max(du_norm, err)
    end do
    du_norm = global_maxval(du_norm)

    ! du_norm = global_maxval(abs(du))
    ! l = global_maxval(abs(u))
    ! if (l > 0) du_norm = du_norm / l
    ! du_norm = 1
    ! l = global_maxval(abs(u))
    ! if (l > 0) du_norm = global_maxval(abs(du)) / l
    ! du_norm = maxval(abs(du) / (this%atol + this%rtol*abs(u)))
    ! du_norm = global_maxval(du_norm)
    ! error = this%norm%compute(u, du)

  end function


  logical function is_converged(this, itr, t, u, du, f_lnorm, tol)
    class(sm_nlsol_model) :: this
    integer, intent(in) :: itr
    real(r8), intent(in) :: t, tol
    real(r8), intent(in), contiguous, target :: u(:), du(:), f_lnorm(:)
    is_converged = this%du_norm(t, u, du) < tol .and. f_lnorm(3) < this%ftol
    !is_converged = f_lnorm(3) < tol
  end function is_converged

end module sm_nlsol_model_type
