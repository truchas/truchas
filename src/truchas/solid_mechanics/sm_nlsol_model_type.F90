module sm_nlsol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use nonlinear_solver_type
  use sm_model_type
  implicit none
  private

  type, extends(nlsol_model) :: sm_nlsol_model
    private
    type(sm_model), pointer :: model => null() ! unowned reference
    type(sm_precon), pointer :: precon => null() ! unowned reference
    type(sm_norm), pointer :: norm => null() ! unowned reference
  contains
    procedure :: init
    !! Deferred procedures from nlsol_model
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
  end type sm_nlsol_model

contains

  subroutine init(this, model, precon, norm)
    class(sm_nlsol_model), intent(out) :: this
    type(sm_model), intent(in), target :: model
    type(sm_precon), intent(in), target :: precon
    type(sm_norm), intent(in), target :: norm
    this%model => model
    this%precon => precon
    this%norm => norm
    ASSERT(associated(model, precon%model))
  end subroutine init

  subroutine compute_f(this, t, u, udot, f)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:), udot(:)
    real(r8), intent(out), contiguous :: f(:)
    call this%model%compute_residual(t, reshape(u, [3,this%model%mesh%nnode]), f)
  end subroutine

  subroutine apply_precon(this, t, u, f)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:)
    real(r8), intent(inout), contiguous :: f(:)
    call this%precon%apply(t, reshape(u, [3,this%model%mesh%nnode]), f)
  end subroutine

  subroutine compute_precon(this, t, u)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:)
    call this%precon%compute(t, reshape(u, [3,this%model%mesh%nnode]))
  end subroutine

  real(r8) function du_norm(this, t, u, du)
    class(sm_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:), du(:)
    real(r8), intent(out) :: error
    error = this%norm%compute(u, du)
  end function

end module sm_nlsol_model_type
