!! HT_2D_IDAESOL_MODEL_TYPE
!!
!! Vector-based IDAESOL adapter for the 2D heat-transfer model.

#include "f90_assert.fpp"

module ht_2d_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use new_idaesol_type, only: idaesol_model
  use vector_class
  use ht_2d_vector_type
  use ht_2d_model_type
  use ht_2d_precon_type
  use ht_2d_norm_type
  use timer_tree_type, only: timer_tree
  implicit none
  private

  type, extends(idaesol_model), public :: ht_2d_idaesol_model
    type(ht_2d_model), pointer :: model => null() ! unowned reference
    type(ht_2d_precon), pointer :: precon => null() ! unowned reference
    type(ht_2d_norm), pointer :: norm => null() ! unowned reference
    type(timer_tree), pointer :: timer => null() ! unowned reference
  contains
    procedure :: init
    procedure :: alloc_vector
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
  end type ht_2d_idaesol_model

contains

  subroutine init(this, model, precon, norm, timer)
    class(ht_2d_idaesol_model), intent(out) :: this
    type(ht_2d_model), intent(in), target :: model
    type(ht_2d_precon), intent(in), target :: precon
    type(ht_2d_norm), intent(in), target :: norm
    type(timer_tree), pointer, intent(in) :: timer

    this%model => model
    this%precon => precon
    this%norm => norm
    this%timer => timer
    ASSERT(associated(this%model, precon%model))
    ASSERT(associated(this%timer))
  end subroutine init


  subroutine alloc_vector(this, vec)
    class(ht_2d_idaesol_model), intent(in) :: this
    class(vector), allocatable, intent(out) :: vec

    type(ht_2d_vector), allocatable :: tmp

    allocate(tmp)
    call this%model%init_vector(tmp)
    call move_alloc(tmp, vec)
  end subroutine alloc_vector


  subroutine compute_f(this, t, u, udot, f)
    class(ht_2d_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u, udot
    class(vector), intent(inout) :: f

    call this%timer%start('residual')
    select type (u)
    class is (ht_2d_vector)
      select type (udot)
      class is (ht_2d_vector)
        select type (f)
        class is (ht_2d_vector)
          call this%model%residual(t, u, udot, f)
        end select
      end select
    end select
    call this%timer%stop('residual')
  end subroutine compute_f


  subroutine apply_precon(this, t, u, f)
    class(ht_2d_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u
    class(vector), intent(inout) :: f

    call this%timer%start('precon apply')
    select type (u)
    class is (ht_2d_vector)
      select type (f)
      class is (ht_2d_vector)
        call this%precon%apply(t, u, f)
      end select
    end select
    call this%timer%stop('precon apply')
  end subroutine apply_precon


  subroutine compute_precon(this, t, u, dt)
    class(ht_2d_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    class(vector), intent(inout) :: u

    call this%timer%start('precon compute')
    select type (u)
    class is (ht_2d_vector)
      call this%precon%compute(t, u, dt)
    end select
    call this%timer%stop('precon compute')
  end subroutine compute_precon


  subroutine du_norm(this, t, u, du, error)
    class(ht_2d_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u, du
    real(r8), intent(out) :: error

    select type (u)
    class is (ht_2d_vector)
      select type (du)
      class is (ht_2d_vector)
        call this%norm%compute(t, u, du, error)
      end select
    end select
  end subroutine du_norm

end module ht_2d_idaesol_model_type
