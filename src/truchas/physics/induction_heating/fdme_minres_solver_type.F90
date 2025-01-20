#include "f90_assert.fpp"

module fdme_minres_solver2_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_lin_op2_class
  use fdme_zvector_type
  use fdme_model_type
  use fdme_precon_class
  use cs_minres_solver2_type
  implicit none
  private

  type, extends(complex_lin_op2) :: fdme_lin_op
    type(fdme_model), pointer :: model => null() ! unowned reference
    class(fdme_precon), allocatable :: my_precon
  contains
    procedure :: matvec
    procedure :: precon
  end type

  type, public :: fdme_minres_solver2
    type(fdme_model), pointer :: model => null() ! unowned reference
    type(cs_minres_solver2) :: minres
    type(fdme_lin_op) :: lin_op
    type(fdme_zvector) :: efield, rhs
  contains
    procedure :: init
    procedure :: solve
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type
    use fdme_precon_gs_type
    use fdme_precon_pcsr_type

    class(fdme_minres_solver2), intent(out) :: this
    type(fdme_model), pointer :: model !TODO: don't make a pointer
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    character(:), allocatable :: precon_type

    this%model => model
    call this%efield%init(model%mesh) ! initialized to 0
    call this%rhs%init(model%mesh)
    call this%minres%init(params)
    this%lin_op%model => model

    plist => params%sublist('precon')
    call plist%get('type', precon_type, stat, errmsg, default='gs')
    if (stat /= 0) return
    select case (precon_type)
    case ('gs') ! block Gauss-Seidel
      allocate(fdme_precon_gs :: this%lin_op%my_precon)
      call this%lin_op%my_precon%init(model, plist, stat, errmsg)
    case ('boomer', 'ilu') ! Hypre BoomerAMG and ILU
      allocate(fdme_precon_pcsr :: this%lin_op%my_precon)
      call this%lin_op%my_precon%init(model, plist, stat, errmsg)
    case ('none')
    case default
      stat = 1
      errmsg = 'invalid type value: ' // precon_type
      return
    end select
    if (stat /= 0) return

  end subroutine

  subroutine solve(this, efield, stat)
    class(fdme_minres_solver2), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:)
    integer, intent(out) :: stat
    this%rhs%w1(:) = this%model%rhs
    if (allocated(this%lin_op%my_precon)) call this%lin_op%my_precon%setup
    call this%minres%solve(this%lin_op, this%rhs, this%efield)
    efield(:) = this%efield%w1
    stat = 0 !FIXME: need to extract from minres
  end subroutine

  subroutine matvec(this, x, y)
    use zvector_class
    class(fdme_lin_op), intent(inout) :: this
    class(zvector) :: x, y
    select type (x)
    type is (fdme_zvector)
      select type (y)
      type is (fdme_zvector)
        call x%gather_offp
        call this%model%A%matvec(x%w1, y%w1)
      end select
    end select
  end subroutine

  subroutine precon(this, x, y)
    use zvector_class
    class(fdme_lin_op), intent(inout) :: this
    class(zvector) :: x, y
    if (allocated(this%my_precon)) then
      select type (x)
      type is (fdme_zvector)
        select type (y)
        type is (fdme_zvector)
          call x%gather_offp
          call this%my_precon%apply(x%w1, y%w1)
          call y%gather_offp ! necessary?
        end select
      end select
    else ! no preconditioning
      call y%copy(x)
    end if
  end subroutine

end module fdme_minres_solver2_type
