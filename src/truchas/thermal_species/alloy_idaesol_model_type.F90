!!
!! ALLOY_IDAESOL_MODEL_TYPE
!!
!! Adapter between the alloy phase-change model and the common IDAESOL
!! interface used by thermal/species transport solvers.
!!

#include "f90_assert.fpp"

module alloy_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use new_idaesol_type, only: idaesol_model
  use vector_class
  use alloy_vector_type
  use alloy_model_type
  use alloy_precon_type
  use alloy_norm_type
  use truchas_timers, only: start_timer, stop_timer
  implicit none
  private

  type, extends(idaesol_model), public :: alloy_idaesol_model
    type(alloy_model), pointer :: model => null() ! unowned reference
    type(alloy_precon), pointer :: precon => null() ! unowned reference
    type(alloy_norm), pointer :: norm => null() ! unowned reference
    real(r8), pointer :: conc(:,:) => null() ! unowned reference
    real(r8), pointer :: conc_dot(:,:) => null() ! unowned reference
  contains
    procedure :: init
    procedure :: alloc_vector
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
  end type

contains

  subroutine init(this, model, precon, norm, conc, conc_dot)
    class(alloy_idaesol_model), intent(out) :: this
    type(alloy_model), intent(in), target :: model
    type(alloy_precon), intent(in), target :: precon
    type(alloy_norm), intent(in), target :: norm
    real(r8), intent(in), target :: conc(:,:), conc_dot(:,:)
    this%model => model
    this%precon => precon
    this%norm => norm
    this%conc => conc
    this%conc_dot => conc_dot
    ASSERT(associated(this%model, precon%model))
  end subroutine

  subroutine alloc_vector(this, vec)
    class(alloy_idaesol_model), intent(in) :: this
    class(vector), allocatable, intent(out) :: vec
    type(alloy_vector), allocatable :: tmp
    allocate(tmp)
    call this%model%init_vector(tmp)
    call move_alloc(tmp, vec)
  end subroutine

  subroutine compute_f(this, t, u, udot, f)
    class(alloy_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u, udot
    class(vector), intent(inout) :: f
    call start_timer('residual')
    select type (u)
    class is (alloy_vector)
      select type (udot)
      class is (alloy_vector)
        select type (f)
        class is (alloy_vector)
          call this%model%compute_f(this%conc, this%conc_dot, t, u, udot, f)
        end select
      end select
    end select
    call stop_timer('residual')
  end subroutine

  subroutine apply_precon(this, t, u, f)
    class(alloy_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u, f
    call start_timer('precon apply')
    select type (u)
    class is (alloy_vector)
      select type (f)
      class is (alloy_vector)
        call this%precon%apply(t, u, f)
      end select
    end select
    call stop_timer('precon apply')
  end subroutine

  subroutine compute_precon(this, t, u, dt)
    class(alloy_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    class(vector), intent(inout) :: u
    type(alloy_vector) :: udot
    call start_timer('precon compute')
    select type (u)
    class is (alloy_vector)
      !! The refactored IDAESOL interface deliberately does not expose an
      !! iterate derivative here.  The alloy Jacobian only needs it for the
      !! optional back-diffusion algebraic block; a zero derivative gives a
      !! consistent first-order preconditioner for that block.
      call udot%init(u)
      call udot%setval(0.0_r8)
      call this%precon%compute(this%conc, this%conc_dot, t, u, udot, dt)
    end select
    call stop_timer('precon compute')
  end subroutine

  subroutine du_norm(this, t, u, du, error)
    class(alloy_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u, du
    real(r8), intent(out) :: error
    select type (u)
    class is (alloy_vector)
      select type (du)
      class is (alloy_vector)
        call this%norm%compute(t, u, du, error)
      end select
    end select
  end subroutine

end module alloy_idaesol_model_type
