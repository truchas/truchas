#include "f90_assert.fpp"

module fdme_gmres_solver_type

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
  contains
    procedure :: matvec
    procedure :: apply_precon
  end type

  type, public :: fdme_gmres_solver
    type(fdme_model), pointer :: model => null() ! unowned reference
    type(my_gmres_solver) :: gmres
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
    call this%gmres%solve(this%model%rhs, efield, stat)
  end subroutine

  subroutine matvec(this, x, Ax)
    class(my_gmres_solver), intent(inout) :: this
    class(vector), intent(inout) :: x, Ax
    select type (x)
    type is (fdme_vector)
      select type (Ax)
      type is (fdme_vector)
        call this%model%matvec(x, Ax)
      end select
    end select
  end subroutine

  subroutine apply_precon(this, x)
    class(my_gmres_solver), intent(inout) :: this
    class(vector), intent(inout) :: x
    select type (x)
    type is (fdme_vector)
      call this%precon%apply(x%array)
    end select
  end subroutine

end module fdme_gmres_solver_type
