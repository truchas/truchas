!!
!! FDME_MUMPS_SOLVER
!!
!! A solver for thefrequency-domain Maxwell equations that solves the linear
!! system directly using the MUMPS library.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fdme_mumps_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use zmumps_solver_type
  use dmumps_solver_type
  use fdme_model_type
  use truchas_timers
  use index_map_type
  implicit none
  private

  type, public :: fdme_mumps_solver
    type(fdme_model), pointer :: model => null() ! unowned reference
    type(zmumps_solver) :: zmumps
    type(dmumps_solver) :: dmumps

    ! used for converting pbsr to pcsr
    type(index_map), pointer :: imap => null()
    complex(r8), allocatable :: bc(:), xc(:)
    real(r8), allocatable :: br(:), xr(:)
    logical :: use_complex_mumps = .false.
  contains
    procedure :: init
    procedure :: solve
    procedure, private :: pbsr_to_pcsr
    procedure, private :: solve_nonmixed
    procedure, private :: solve_mixed
    final :: fdme_mumps_solver_delete
  end type

contains

  elemental subroutine fdme_mumps_solver_delete(this)
    type(fdme_mumps_solver), intent(inout) :: this
    if (associated(this%imap)) deallocate(this%imap)
  end subroutine fdme_mumps_solver_delete


  subroutine init(this, model, params, stat, errmsg)
    use parameter_list_type
    class(fdme_mumps_solver), intent(out) :: this
    type(fdme_model), pointer :: model !TODO: don't make a pointer
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    stat = 0
    this%model => model
    if (this%use_complex_mumps) then
      call this%zmumps%init(0, stat) !, 3)
    else
      call this%dmumps%init(0, stat) !, 3)
    end if
  end subroutine

  subroutine solve(this, efield, stat)
    class(fdme_mumps_solver), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:)
    integer, intent(out) :: stat
    if (this%model%use_mixed_form) then
      call this%solve_mixed(this%model%rhs, efield, stat)
    else
      call this%solve_nonmixed(this%model%rhs, efield, stat)
    end if
  end subroutine solve


  ! Solves the mixed formulation
  subroutine solve_mixed(this, rhs, efield, stat)

    use truchas_logging_services

    class(fdme_mumps_solver), intent(inout) :: this
    complex(r8), intent(inout) :: rhs(:), efield(:)
    integer, intent(out) :: stat

    integer :: Ne, Nn

    ASSERT(.not.this%use_complex_mumps)

    call start_timer("mumps-setup")
    Ne = this%model%mesh%nedge_onp
    Nn = this%model%mesh%nnode_onp
    if (.not.allocated(this%br)) allocate(this%br(2*Ne + 2*Nn))
    if (.not.allocated(this%xr)) allocate(this%xr(2*Ne + 2*Nn))

    call this%dmumps%setup(this%model%Am, stat)
    if (stat /= 0) return
    call stop_timer("mumps-setup")

    call start_timer("mumps-solve")
    this%br(1:Ne) = rhs(1:Ne)%re
    this%br(Ne+1:2*Ne) = rhs(1:Ne)%im
    this%br(2*Ne+1:) = 0
    call this%dmumps%solve(this%br, this%xr, stat)
    efield(1:Ne)%re = this%xr(1:Ne)
    efield(1:Ne)%im = this%xr(Ne+1:2*Ne)
    call stop_timer("mumps-solve")

  end subroutine solve_mixed


  subroutine pbsr_to_pcsr(this, A, Ac)

    use parallel_communication, only: is_IOP
    use pbsr_matrix_type
    use pcsr_matrix_type

    class(fdme_mumps_solver), intent(inout) :: this
    type(pbsr_matrix), intent(in) :: A
    type(pcsr_matrix), intent(out) :: Ac

    type(pcsr_graph), pointer :: g => null()
    integer, allocatable :: nvars(:)
    integer :: i, ii, ic, xj, j, jj, jc

    if (.not.associated(this%imap)) allocate(this%imap)

    ! convert pbsr matrix to pcsr
    allocate(g, nvars(merge(A%graph%row_imap%global_size, 0, is_IOP)))
    nvars = A%bsize
    call this%imap%init(A%graph%row_imap, nvars)

    call g%init(this%imap)
    do i = 1, A%nrow
      do xj = A%graph%xadj(i), A%graph%xadj(i+1)-1
        j = A%graph%adjncy(xj)
        do jj = 1, A%bsize
          jc = A%bsize * (j - 1) + jj
          do ii = 1, A%bsize
            ic = A%bsize * (i - 1) + ii
            call g%add_edge(ic, jc)
          end do
        end do
      end do
    end do
    call g%add_complete

    call Ac%init(g, take_graph=.true.)
    call Ac%set_all(0.0_r8)
    do i = 1, A%nrow
      do xj = A%graph%xadj(i), A%graph%xadj(i+1)-1
        j = A%graph%adjncy(xj)
        do jj = 1, A%bsize
          jc = A%bsize * (j - 1) + jj
          do ii = 1, A%bsize
            ic = A%bsize * (i - 1) + ii
            call Ac%set(ic, jc, A%values(ii, jj, xj))
          end do
        end do
      end do
    end do

  end subroutine pbsr_to_pcsr


  subroutine solve_nonmixed(this, rhs, efield, stat)

    use pbsr_matrix_type
    use pcsr_matrix_type

    class(fdme_mumps_solver), intent(inout) :: this
    complex(r8), intent(inout) :: rhs(:), efield(:)
    integer, intent(out) :: stat

    type(pcsr_matrix) :: Ac

    ! ASSERT(size(b) >= this%mumps%nloc_rhs)
    ! ASSERT(size(x) >= this%mumps%lsol_loc)

    call start_timer("mumps-setup")
    if (this%use_complex_mumps) then
      call this%zmumps%setup(this%model%A, stat)
    else
      call this%pbsr_to_pcsr(this%model%A2, Ac)
      call this%dmumps%setup(Ac, stat)
    end if
    if (stat /= 0) return
    call stop_timer("mumps-setup")

    call start_timer("mumps-solve")
    if (this%use_complex_mumps) then
      call this%zmumps%solve(rhs, efield, stat)
    else
      if (.not.allocated(this%br)) allocate(this%br(this%model%A2%bsize*this%model%A2%nrow))
      if (.not.allocated(this%xr)) allocate(this%xr(this%model%A2%bsize*this%model%A2%nrow))
      this%br(1::2) = rhs(:)%re
      this%br(2::2) = rhs(:)%im
      call this%dmumps%solve(this%br, this%xr, stat)
      efield(:)%re = this%xr(1::2)
      efield(:)%im = this%xr(2::2)
    end if
    call stop_timer("mumps-solve")

  end subroutine solve_nonmixed

end module fdme_mumps_solver_type
