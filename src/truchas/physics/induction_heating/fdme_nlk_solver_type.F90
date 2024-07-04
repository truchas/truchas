#include "f90_assert.fpp"

module fdme_nlk_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nlk_solver_class
  use fdme_vector_type
  use fdme_model_type
  use fdme_precon_class
  use vector_class
  implicit none
  private

  type, extends(nlk_solver) :: my_nlk_solver
    type(fdme_model), pointer :: model => null() ! unowned reference
    class(fdme_precon), pointer :: precon => null() ! unowned reference
    real(r8), allocatable :: b(:,:)
  contains
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
  end type

  type, public :: fdme_nlk_solver
    type(my_nlk_solver) :: nlk
  contains
    procedure :: init
    procedure :: solve
  end type

contains

  subroutine init(this, vec, model, precon, params, stat, errmsg)
    use parameter_list_type
    class(fdme_nlk_solver), intent(out) :: this
    type(fdme_vector), intent(in) :: vec
    type(fdme_model), intent(in), target :: model
    class(fdme_precon), intent(in), target :: precon
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%nlk%init(vec, params, stat, errmsg)
    if (stat /= 0) return
    this%nlk%model => model
    this%nlk%precon => precon
  end subroutine

  subroutine solve(this, efield, stat)
    class(fdme_nlk_solver), intent(inout) :: this
    type(fdme_vector), intent(inout) :: efield
    integer, intent(out) :: stat
    call this%nlk%solve(efield, stat)
  end subroutine

  subroutine compute_f(this, u, f)
    class(my_nlk_solver) :: this
    class(vector), intent(inout) :: u, f
    select type (u)
    type is (fdme_vector)
      select type (f)
      type is (fdme_vector)
        !call this%model%compute_f(u%array, f%array, ax=.false.)
        call this%model%compute_f(u, f, ax=.false.)
      end select
    end select
  end subroutine

  subroutine apply_precon(this, u, f)
    class(my_nlk_solver) :: this
    class(vector), intent(inout) :: u, f
    select type (u)
    type is (fdme_vector)
      select type (f)
      type is (fdme_vector)
        this%b = f%array  !this allocates b
        f%array = 0.0_r8
        ASSERT(all(ieee_is_finite(this%b)))
        call this%precon%apply(this%b, f%array)
        ASSERT(all(ieee_is_finite(f%array)))
      end select
    end select
  end subroutine

  subroutine compute_precon(this, u)
    class(my_nlk_solver) :: this
    class(vector), intent(inout) :: u
    call this%precon%setup
  end subroutine

end module fdme_nlk_solver_type
