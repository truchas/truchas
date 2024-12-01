#include "f90_assert.fpp"

module fdme_mumps_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use zmumps_solver_type
  use dmumps_solver_type
  use fdme_vector_type
  use fdme_model_type
  use truchas_timers
  use index_map_type
  implicit none
  private

  type, public :: fdme_mumps_solver
    type(fdme_model), pointer :: model => null() ! unowned reference
    type(zmumps_solver) :: zmumps
    type(dmumps_solver) :: dmumps
    type(fdme_vector) :: efield, rhs

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
    call this%efield%init(model%mesh) ! initialized to 0
    call this%rhs%init(model%mesh)
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
    this%rhs%array(1,:) = this%model%rhs%re
    this%rhs%array(2,:) = this%model%rhs%im
    if (this%model%use_mixed_form) then
      call this%solve_mixed(this%rhs, this%efield, stat)
    else
      call this%solve_nonmixed(this%rhs, this%efield, stat)
    end if
    efield%re = this%efield%array(1,:)
    efield%im = this%efield%array(2,:)
  end subroutine solve


  ! Solves the mixed formulation
  subroutine solve_mixed(this, rhs, efield, stat)

    use truchas_logging_services

    class(fdme_mumps_solver), intent(inout) :: this
    type(fdme_vector), intent(inout) :: rhs, efield
    integer, intent(out) :: stat

    integer :: Ne, Nn

    INSIST(.false.)

    ! call start_timer("mumps-setup")
    ! Ne = this%model%mesh%nedge_onp
    ! Nn = this%model%mesh%nnode_onp
    ! if (.not.allocated(this%bc)) allocate(this%bc(2*Ne + 2*Nn))
    ! if (.not.allocated(this%xc)) allocate(this%xc(2*Ne + 2*Nn))

    ! ! ! DEBUGGING
    ! ! if (.not.allocated(this%xc)) allocate(this%xc(2*this%model%mesh%nedge + 2*this%model%mesh%nnode))
    ! ! block
    ! !   use parallel_communication, only: is_iop, this_pe
    ! !   real(r8) :: norm, dx(3)
    ! !   real(r8), allocatable :: xt(:)
    ! !   integer :: i, j, xi, e1, e2, xe1, xe2
    ! !   allocate(xt(size(this%xc)))
    ! !   this%xc = 0
    ! !   do i = 1, this%model%mesh%nedge_onp
    ! !     dx = this%model%mesh%x(:,this%model%mesh%enode(2,i)) - this%model%mesh%x(:,this%model%mesh%enode(1,i))
    ! !     !if (dot_product(dx, [1.0_r8, 1.0_r8, 1.0_r8]) < 0) dx = -dx
    ! !     !dx = dx / norm2(dx)
    ! !     this%xc(i) = sign(1.0_r8, dot_product(dx, [1.0_r8, 1.0_r8, 1.0_r8]))

    ! !     ! e2 = this%model%ic%eval(i,1)
    ! !     ! this%xc(e2) = sign(1.0_r8, dot_product(dx, [1.0_r8, 1.0_r8, 1.0_r8]))

    ! !     this%xc(Ne+i) = this%xc(i)
    ! !   end do

    ! !   call this%model%Am%graph%row_imap%gather_offp(this%xc)
    ! !   call this%model%Am%matvec(this%xc, xt)

    ! !   ! do j = 1, this%model%mesh%ncell
    ! !   !   if (21 /= this%model%mesh%xcell(j)) cycle

    ! !   !   if (this_pe == 3 .or. this_pe == 1) then
    ! !   !     do xe1 = 1, 1 !6
    ! !   !       e1 = this%model%mesh%cedge(xe1,j)
    ! !   !       if (e1 > this%model%mesh%nedge_onp) cycle
    ! !   !       do xe2 = 1, 6
    ! !   !         e2 = this%model%mesh%cedge(xe2,j)
    ! !   !         ! xi = this%model%A2%graph%index(e1, e2)
    ! !   !         ! if (xi > 0) print *, "debugA: ", e1, e2, this%model%A2%values(xi)
    ! !   !         xi = this%model%Am%graph%index(this%model%ic%eval(e1,1), this%model%ic%eval(e2,1))
    ! !   !         if (xi > 0) print *, "debugAm: ", e1, e2, this%model%Am%values(xi)
    ! !   !       end do
    ! !   !     end do
    ! !   !   end if

    ! !   !   do xi = 1, 6
    ! !   !     i = this%model%mesh%cedge(xi,j)
    ! !   !     e2 = this%model%ic%eval(i,1)
    ! !   !     print '(a,5i6,2es15.5)', "debugXI: ", this_pe, i, this%model%mesh%edge_imap%global_index(i), &
    ! !   !         e2, this%model%am%graph%row_imap%global_index(e2), &
    ! !   !         this%xc(e2), xt(e2)
    ! !   !   end do
    ! !   ! end do

    ! !   norm = global_norm2(xt(:2*this%model%mesh%nedge_onp))
    ! !   if (is_iop) print '(a,es20.10)', "debug: ", norm
    ! !   INSIST(.false.)
    ! ! end block

    ! !call tls_info("mumps setup")
    ! call this%mumps%setup(this%model%Am, stat)
    ! if (stat /= 0) return
    ! call stop_timer("mumps-setup")

    ! call start_timer("mumps-solve")
    ! this%bc(1:Ne) = rhs%array(1,1:Ne)
    ! this%bc(Ne+1:2*Ne) = rhs%array(2,1:Ne)
    ! this%bc(2*Ne+1:) = 0
    ! !call tls_info("mumps solve")
    ! call this%mumps%solve(this%bc, this%xc, stat)
    ! efield%array(1,1:Ne) = this%xc(1:Ne)
    ! efield%array(2,1:Ne) = this%xc(Ne+1:2*Ne)
    ! call stop_timer("mumps-solve")

  end subroutine solve_mixed


  ! TODO: needs a complex pbsr
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
    type(fdme_vector), intent(inout) :: rhs, efield
    integer, intent(out) :: stat

    type(pcsr_matrix) :: Ac
    integer :: i, ii, ic

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

    ! ! DEBUGGING
    ! block
    !   use parallel_communication, only: is_iop
    !   real(r8) :: norm, dx(3)
    !   real(r8), allocatable :: xt(:)
    !   integer :: i, i1, i2
    !   allocate(xt(2*this%model%mesh%nedge))
    !   do i = 1, this%model%mesh%nedge_onp
    !     i1 = 2 * (i - 1) + 1
    !     i2 = 2 * (i - 1) + 2
    !     dx = this%model%mesh%x(:,this%model%mesh%enode(2,i)) - this%model%mesh%x(:,this%model%mesh%enode(1,i))
    !     !dx = dx / norm2(dx)
    !     this%xc(i1) = sign(1.0_r8, dot_product(dx, [1.0_r8, 1.0_r8, 1.0_r8]))
    !     this%xc(i2) = this%xc(i1)
    !   end do
    !   call Ac%graph%row_imap%gather_offp(this%xc)
    !   call Ac%matvec(this%xc, xt)
    !   norm = global_norm2(xt(:2*this%model%mesh%nedge_onp))
    !   if (is_iop) print '(a,es20.10)', "debug: ", norm
    !   INSIST(.false.)
    ! end block

    call start_timer("mumps-solve")
    if (this%use_complex_mumps) then
      if (.not.allocated(this%xc)) allocate(this%xc(this%model%A%nrow))
      this%xc(:)%re = efield%array(1,:)
      this%xc(:)%im = efield%array(2,:)
      call this%zmumps%solve(this%model%rhs, this%xc, stat)
      efield%array(1,:) = this%xc(:)%re
      efield%array(2,:) = this%xc(:)%im
    else
      if (.not.allocated(this%br)) allocate(this%br(this%model%A2%bsize*this%model%A2%nrow))
      if (.not.allocated(this%xr)) allocate(this%xr(this%model%A2%bsize*this%model%A2%nrow))

      do i = 1, this%model%A2%nrow
        do ii = 1, this%model%A2%bsize
          ic = this%model%A2%bsize * (i - 1) + ii
          this%br(ic) = rhs%array(ii, i)
        end do
      end do

      call this%dmumps%solve(this%br, this%xr, stat)

      do i = 1, this%model%A2%nrow
        do ii = 1, this%model%A2%bsize
          ic = this%model%A2%bsize * (i - 1) + ii
          efield%array(ii, i) = this%xr(ic)
        end do
      end do
    end if
    call stop_timer("mumps-solve")

  end subroutine solve_nonmixed

end module fdme_mumps_solver_type
