!!
!! FDME_MINRES_SOLVER
!!
!! A solver for the frequency-domain Maxwell equations that solves the linear
!! system iteratively using the preconditioned CS-MINRES method for complex
!! symmetric systems.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fdme_minres_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_lin_op_class
  use fdme_model_type
  use cs_minres_solver_type
  use pcsr_precon_class
  use pcsr_matrix_type
  use imap_zvector_type
  use truchas_timers
  implicit none
  private

  type, extends(complex_lin_op), public :: fdme_minres_solver
    type(fdme_model), pointer :: model => null() ! unowned reference
    type(cs_minres_solver) :: minres
    class(pcsr_precon), allocatable :: pc
    type(pcsr_matrix), pointer :: pc_mtx => null() ! pointer to avoid dangling pointer
    type(imap_zvector) :: efield, rhs
    real(r8) :: beta
  contains
    procedure :: init
    procedure :: solve
    ! deferred procedures from complex_lin_op class
    procedure :: matvec
    procedure :: precon
    final :: fdme_minres_solver_delete
  end type

contains

  impure elemental subroutine fdme_minres_solver_delete(this)
    type(fdme_minres_solver), intent(inout) :: this
    if (associated(this%pc_mtx)) deallocate(this%pc_mtx)
  end subroutine

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type
    use pcsr_precon_factory, only: alloc_pcsr_precon

    class(fdme_minres_solver), intent(out) :: this
    type(fdme_model), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(parameter_list), pointer :: plist
    character(:), allocatable :: precon_type

    call start_timer('cs-minres-solver')

    this%model => model
    call this%efield%init(model%mesh%edge_imap) ! initialized to 0
    call this%rhs%init(mold=this%efield)
    call params%get('print-level', n, stat, errmsg, default=0)
    if (stat /= 0) return
    call params%set('verbose', n > 0)
    call this%minres%init(params)

    plist => params%sublist('precon')
    call plist%get('type', precon_type, stat, errmsg, default='ssor')
    if (stat /= 0) return
    if (precon_type == 'none') return

    !! Create the preconditioner
    call start_timer('preconditioner')
    allocate(this%pc_mtx)
    call this%pc_mtx%init(this%model%A%graph, take_graph=.false.)
    call plist%get('beta', this%beta, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return
    call alloc_pcsr_precon(this%pc, this%pc_mtx, plist, stat, errmsg)
    if (stat /= 0) return
    call stop_timer('preconditioner')

    call stop_timer('cs-minres-solver')

  end subroutine init

  subroutine solve(this, efield, stat, errmsg)

    class(fdme_minres_solver), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: msg

    call start_timer('cs-minres-solver')

    !! Setup the preconditioner
    if (allocated(this%pc)) then
      call start_timer('preconditioner')
      this%pc_mtx%values(:) = this%model%A%values%re + this%beta * this%model%M%values%re
      call this%pc%compute
      call stop_timer('preconditioner')
    end if

    this%rhs%v(:) = this%model%rhs

    call this%minres%solve(this, this%rhs, this%efield, stat, msg)
    stat = merge(0, 1, stat >= 0) ! for minres stat >= 0 is success and < 0 failure
    if (stat /= 0) errmsg = msg
    call this%efield%gather_offp
    efield(:) = this%efield%v

    block
      use truchas_logging_services
      character(80) :: msg
      if (stat == 0) then
        write(msg,'(a,i0,a)') 'CS-MINRES converged: ', this%minres%num_iter, ' iterations'
        call TLS_info(msg)
      else
        write(msg,'(a,": ",i0,a)') 'CS-MINRES failed: '//errmsg, this%minres%num_iter, ' iterations'
        call TLS_info(msg)
        errmsg = 'CS-MINRES: ' // errmsg
      end if
    end block

    call stop_timer('cs-minres-solver')

  end subroutine solve

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

  end subroutine matvec

  subroutine precon(this, x, y)

    use zvector_class

    class(fdme_minres_solver), intent(inout) :: this
    class(zvector) :: x, y

    call y%copy(x)
    if (.not.allocated(this%pc)) return

    call start_timer('preconditioner')

    select type (x)
    type is (imap_zvector)
      select type (y)
      type is (imap_zvector)
        call y%gather_offp
#if defined(INTEL_BUG20250605) || defined(GNU_PR119986)
        call workaround(y%v)
#else
        call this%pc%apply(y%v%re)
        call this%pc%apply(y%v%im)
#endif
      end select
    end select

    call stop_timer('preconditioner')

#if defined(INTEL_BUG20250605) || defined(GNU_PR119986)
  contains

    subroutine workaround(z)
      complex(r8), intent(inout) :: z(:)
      call this%pc%apply(z%re)
      call this%pc%apply(z%im)
    end subroutine
#endif

  end subroutine precon

end module fdme_minres_solver_type
