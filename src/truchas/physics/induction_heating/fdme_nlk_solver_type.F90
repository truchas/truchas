#include "f90_assert.fpp"

module fdme_nlk_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
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
  contains
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
  end type

  type, public :: fdme_nlk_solver
    type(my_nlk_solver) :: nlk
    type(fdme_vector) :: efield
  contains
    procedure :: init
    procedure :: solve
  end type

contains

  subroutine init(this, model, precon, params, stat, errmsg)
    use parameter_list_type
    class(fdme_nlk_solver), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    class(fdme_precon), intent(in), target :: precon
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%efield%init(model%mesh) ! initialized to 0
    call this%nlk%init(this%efield, params, stat, errmsg)
    if (stat /= 0) return
    this%nlk%model => model
    this%nlk%precon => precon
  end subroutine

  subroutine solve(this, efield, stat)
    class(fdme_nlk_solver), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:)
    integer, intent(out) :: stat
    call this%nlk%solve(this%efield, stat)
    efield%re = this%efield%array(1,:)
    efield%im = this%efield%array(2,:)
  end subroutine

  subroutine compute_f(this, u, f)
    class(my_nlk_solver) :: this
    class(vector), intent(inout) :: u, f
    select type (u)
    type is (fdme_vector)
      select type (f)
      type is (fdme_vector)
        call u%gather_offp
        call this%model%A2%matvec(u%array, f%array)
        f%array(1,:) = f%array(1,:) - this%model%rhs%re
        f%array(2,:) = f%array(2,:) - this%model%rhs%im
      end select
    end select
  end subroutine

  subroutine apply_precon(this, u, f)
    class(my_nlk_solver) :: this
    class(vector), intent(inout) :: u ! not used
    class(vector), intent(inout) :: f
    select type (f)
    type is (fdme_vector)
      call this%precon%apply(f%array)
    end select
  end subroutine

  subroutine compute_precon(this, u)
    class(my_nlk_solver) :: this
    class(vector), intent(inout) :: u
    call this%precon%setup
  end subroutine

end module fdme_nlk_solver_type
