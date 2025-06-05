#include "f90_assert.fpp"

module fdme_minres_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_lin_op_class
  use imap_zvector_type
  use fdme_model_type
  use fdme_precon_class
  use cs_minres_solver_type
  implicit none
  private

  type, extends(complex_lin_op), public :: fdme_minres_solver
    type(fdme_model), pointer :: model => null() ! unowned reference
    class(fdme_precon), allocatable :: my_precon
    type(cs_minres_solver) :: minres
    type(imap_zvector) :: efield, rhs
  contains
    procedure :: init
    procedure :: solve
    ! deferred procedures from complex_lin_op class
    procedure :: matvec
    procedure :: precon
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type
    use fdme_precon_pcsr_type

    class(fdme_minres_solver), intent(out) :: this
    type(fdme_model), pointer :: model !TODO: don't make a pointer
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    character(:), allocatable :: precon_type

    this%model => model
    call this%efield%init(model%mesh%edge_imap) ! initialized to 0
    call this%rhs%init(model%mesh%edge_imap)
    call this%minres%init(params)

    plist => params%sublist('precon')
    call plist%get('type', precon_type, stat, errmsg, default='gs')
    if (stat /= 0) return
    select case (precon_type)
    case ('boomer','ssor') ! Hypre BoomerAMG
      allocate(fdme_precon_pcsr :: this%my_precon)
      call this%my_precon%init(model, plist, stat, errmsg)
    case ('none')
    case default
      stat = 1
      errmsg = 'invalid type value: ' // precon_type
      return
    end select
    if (stat /= 0) return

  end subroutine

  subroutine solve(this, efield, stat, errmsg)
    use string_utilities, only: i_to_c
    class(fdme_minres_solver), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    character(:), allocatable :: msg
    this%rhs%v(:) = this%model%rhs
    if (allocated(this%my_precon)) call this%my_precon%setup
    call this%minres%solve(this, this%rhs, this%efield, stat, msg)
    !TODO: add solver summary output
    stat = merge(0, 1, stat >= 0) ! for minres stat >= 0 is success and < 0 failure
    if (stat /= 0) errmsg = msg
    efield(:) = this%efield%v
  end subroutine

  subroutine matvec(this, x, y)
    use zvector_class
    class(fdme_minres_solver), intent(inout) :: this
    class(zvector) :: x, y
    select type (x)
    type is (imap_zvector)
      select type (y)
      type is (imap_zvector)
        call x%gather_offp
        call this%model%A%matvec(x%v, y%v)
      end select
    end select
  end subroutine

  subroutine precon(this, x, y)
    use zvector_class
    class(fdme_minres_solver), intent(inout) :: this
    class(zvector) :: x, y
    if (allocated(this%my_precon)) then
      select type (x)
      type is (imap_zvector)
        select type (y)
        type is (imap_zvector)
          call x%gather_offp
          call this%my_precon%apply(x%v, y%v)
          call y%gather_offp ! necessary?
        end select
      end select
    else ! no preconditioning
      call y%copy(x)
    end if
  end subroutine

end module fdme_minres_solver_type
