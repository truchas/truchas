#include "f90_assert.fpp"

module fdme_gmres_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use gmres_left_solver_class
  use fdme_vector_type
  use fdme_model_type
  use fdme_precon_class
  use vector_class
  implicit none
  private

  type, extends(gmres_left_solver) :: my_gmres_solver
    type(fdme_model), pointer :: model => null() ! unowned reference
    class(fdme_precon), pointer :: precon => null() ! unowned reference
    real(r8), allocatable :: b(:,:)
  contains
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
  end type

  type, public :: fdme_gmres_solver
    type(my_gmres_solver) :: gmres
    type(fdme_model), pointer :: model => null() ! unowned reference
  contains
    procedure :: init
    procedure :: solve
  end type

contains

  subroutine init(this, vec, model, precon, params, stat, errmsg)
    use parameter_list_type
    class(fdme_gmres_solver), intent(out) :: this
    type(fdme_vector), intent(in) :: vec
    type(fdme_model), pointer :: model !TODO: don't make a pointer
    class(fdme_precon), intent(in), target :: precon
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    this%model => model
    call this%gmres%init(vec, params, stat, errmsg)
    this%gmres%model => model
    this%gmres%precon => precon
  end subroutine

  subroutine solve(this, efield, stat)
    class(fdme_gmres_solver), intent(inout) :: this
    type(fdme_vector), intent(inout) :: efield
    integer, intent(out) :: stat
    select type (rhs => this%gmres%rhs)
    type is (fdme_vector)
      rhs%array(:,:) = this%model%rhs
    end select
    call this%gmres%solve(efield, stat)
  end subroutine

  subroutine compute_f(this, u, f, ax)
    class(my_gmres_solver), intent(inout) :: this
    class(vector), intent(in) :: u
    class(vector), intent(inout) :: f
    logical, intent(in), optional :: ax
    select type (u)
    type is (fdme_vector)
      select type (f)
      type is (fdme_vector)
        call this%model%compute_f(u%array, f%array, ax)
      end select
    end select
  end subroutine

  subroutine apply_precon(this, u, f)
    class(my_gmres_solver), intent(inout) :: this
    class(vector), intent(in) :: u  ! not used
    class(vector), intent(inout) :: f
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
    class(my_gmres_solver), intent(inout) :: this
    class(vector), intent(in) :: u
    call this%precon%setup
  end subroutine

end module fdme_gmres_solver_type
