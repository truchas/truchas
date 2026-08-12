!! HT_2D_IDAESOL_MODEL_TYPE
!!
!! Vector-based IDAESOL adapter for the 2D heat-transfer model.

#include "f90_assert.fpp"

module ht_2d_idaesol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use new_idaesol_type, only: idaesol_model
  use vector_class
  use ht_2d_vector_type
  use HT_2d_model_type
  use HT_2d_precon_type
  use HT_2d_norm_type
  implicit none
  private

  type, extends(idaesol_model), public :: ht_2d_idaesol_model
    type(HT_2d_model), pointer :: model => null() ! unowned reference
    type(HT_2d_precon), pointer :: precon => null() ! unowned reference
    type(HT_2d_norm), pointer :: norm => null() ! unowned reference
  contains
    procedure :: init
    procedure :: alloc_vector
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
  end type ht_2d_idaesol_model

contains

  subroutine init(this, model, precon, norm)
    class(ht_2d_idaesol_model), intent(out) :: this
    type(HT_2d_model), intent(in), target :: model
    type(HT_2d_precon), intent(in), target :: precon
    type(HT_2d_norm), intent(in), target :: norm

    this%model => model
    this%precon => precon
    this%norm => norm
    ASSERT(associated(this%model, precon%model))
  end subroutine init


  subroutine alloc_vector(this, vec)
    class(ht_2d_idaesol_model), intent(in) :: this
    class(vector), allocatable, intent(out) :: vec

    allocate(ht_2d_vector :: vec)
    select type (vec)
    class is (ht_2d_vector)
      call vec%init(this%model%mesh)
    end select
  end subroutine alloc_vector


  subroutine compute_f(this, t, u, udot, f)
    class(ht_2d_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u, udot
    class(vector), intent(inout) :: f

    select type (u)
    class is (ht_2d_vector)
      select type (udot)
      class is (ht_2d_vector)
        select type (f)
        class is (ht_2d_vector)
          call this%model%compute_f(t, u, udot, f)
        end select
      end select
    end select
  end subroutine compute_f


  subroutine apply_precon(this, t, u, f)
    class(ht_2d_idaesol_model) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u
    class(vector), intent(inout) :: f

    select type (u)
    class is (ht_2d_vector)
      select type (f)
      class is (ht_2d_vector)
        call this%precon%apply(f)
      end select
    end select
  end subroutine apply_precon


  subroutine compute_precon(this, t, u, dt)
    class(ht_2d_idaesol_model) :: this
    real(r8), intent(in) :: t, dt
    class(vector), intent(inout) :: u

    select type (u)
    class is (ht_2d_vector)
      call this%precon%compute(t, u, dt)
    end select
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
        call this%norm%compute(u, du, error)
      end select
    end select
  end subroutine du_norm

end module ht_2d_idaesol_model_type
