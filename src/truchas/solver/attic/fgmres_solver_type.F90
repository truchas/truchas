!!
!! FGMRES_SOLVER_TYPE
!!
!! This module defines a preconditioned flexible GMRES solver type, and a
!! preconditioner class for users to extend into their own concrete types.
!!
!! Zach Jibben <zjibben@lanl.gov>
!!
!! References:
!!
!! 1. Saad, A Flexible Inner-Outer Preconditioned GMRES Algorithm, 1993.
!!
!! 2. Saad, Iterative Methods for Sparse Linear Systems, 2003.
!!
!! 3. van der Vorst, Iterative Krylov Methods for Large Linear Systems, 2003.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fgmres_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use parallel_communication, only: global_dot_product
  use nlsol_type, only: nlsol_model
  !use truchas_logging_services
  !use truchas_timers
  implicit none
  private

  public :: nlsol_model ! re-export

  type, public :: fgmres_solver
    private
    class(nlsol_model), pointer :: model => null() ! unowned reference
    integer :: nrows_onP, krylov_dim, iter_pc, max_iter, iter
    real(r8) :: res_norm, tol
    real(r8), allocatable :: e1(:), h(:,:), r(:), w(:), v(:,:), z(:,:), work(:), udot(:)
  contains
    procedure :: init => fgmres_solver_init
    !procedure :: setup => fgmres_solver_setup
    procedure :: solve => fgmres_solver_solve
    procedure :: metrics_string
  end type fgmres_solver

contains

  subroutine fgmres_solver_init(this, model, params, stat, errmsg)

    use parameter_list_type
    external :: dgels

    class(fgmres_solver), intent(out) :: this
    class(nlsol_model), intent(in), target :: model
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context
    integer :: lwork

    stat = 0 ! TODO: error handling / argument verification
    this%model => model
    this%nrows_onP = model%size()

    call params%get('gmres-krylov-dim', this%krylov_dim, default=5)
    call params%get('abs-tol', this%tol, default=1d-8)
    call params%get('max-iter', this%max_iter, default=20)

    allocate(this%e1(this%krylov_dim+1), this%h(this%krylov_dim+1,this%krylov_dim), &
        this%r(this%nrows_onP), this%w(this%nrows_onP), &
        this%v(this%nrows_onP,this%krylov_dim+1), this%z(this%nrows_onP,this%krylov_dim), &
        this%udot(this%nrows_onP))
    this%e1 = 0
    this%e1(1) = 1
    this%h = 0
    this%udot = 0

    ! get the LAPACK work size
    allocate(this%work(1))
    lwork = -1
    call dgels('N', size(this%h, dim=1), size(this%h, dim=2), 1, this%h, size(this%h, dim=1), &
        this%e1, size(this%h, dim=1), &
        this%work, lwork, stat)
    ASSERT(stat == 0)
    lwork = this%work(1)
    deallocate(this%work)
    allocate(this%work(lwork))

  end subroutine fgmres_solver_init


  subroutine fgmres_solver_solve(this, t, h, u0, u, errc)

    class(fgmres_solver), intent(inout) :: this
    real(r8), intent(in) :: t, h, u0(:)
    real(r8), intent(inout) :: u(:)
    integer,  intent(out) :: errc

    integer :: iter, i, j
    real(r8) :: y(this%krylov_dim), rnorm

    ASSERT(this%max_iter > 0)
    ASSERT(associated(this%model))
    ASSERT(size(u) >= this%nrows_onP)

    errc = 0
    this%iter_pc = 0
    this%z = 0
    this%h = 0
    call this%model%compute_precon(t, u, h)

    do iter = 1, this%max_iter+1
      !! compute the residual
      call this%model%compute_f(t, u, this%udot, this%r)

      this%res_norm = global_norm2(this%r)
      rnorm = global_norm2(this%r / this%model%rhs)
      print '("FGMRES iter, res anorm, rnorm, tol: ",i6,3es13.3)', &
          iter, this%res_norm, rnorm, this%tol
      if (this%res_norm < this%tol) then
        this%iter = iter
        return
      end if
      this%v(:,1) = this%r / this%res_norm

      do j = 1, this%krylov_dim
        this%iter_pc = this%iter_pc + 1
        this%z(:,j) = this%v(:,j)
        call this%model%apply_precon(t, u, this%z(:,j))
        call this%model%compute_f(t, this%z(:,j), this%udot, this%w, ax=.true.)

        do i = 1, j
          this%h(i,j) = global_dot_product(this%w, this%v(:,i))
          this%w = this%w - this%h(i,j) * this%v(:,i)
        end do
        this%h(j+1,j) = global_norm2(this%w)
        this%v(:,j+1) = this%w / this%h(j+1,j)

        !! TODO: van der Vorst has some extra logic for early exit in this loop.
        !! Should consider implementing this.
      end do

      this%e1 = 0
      this%e1(1) = this%res_norm
      call argmin(this%h, this%e1, this%work, y, errc)
      if (errc /= 0) then
        errc = 3
        return
      end if
      u = u + matmul(this%z, y)
    end do

    errc = 1

  end subroutine fgmres_solver_solve


  real(r8) function global_norm2(x)
    real(r8), intent(in) :: x(:)
    global_norm2 = sqrt(global_dot_product(x, x))
  end function global_norm2


  !! Return the y that minimizes ||beta*e1 - A*y||_2.
  subroutine argmin(A, b, work, y, status)

    external dgels

    real(r8), intent(in) :: A(:,:)
    real(r8), intent(inout) :: b(:), work(:)
    real(r8), intent(out) :: y(:)
    integer, intent(out) :: status

    ASSERT(size(A, dim=2) == size(y))
    ASSERT(size(A, dim=1) == size(b))
    ASSERT(size(b) == size(y) + 1)

    call dgels('N', size(A, dim=1), size(A, dim=2), 1, A, size(A, dim=1), &
        b, size(A, dim=1), &
        work, size(work), status)
    y = b(:size(y))

  end subroutine argmin


  pure function metrics_string(this) result(string)
    class(fgmres_solver), intent(in) :: this
    character(:), allocatable :: string
    character(128) :: buffer
    write(buffer,'(i4," (PC), ",i4," (FGMRES), ",es14.4," (|r|)")') &
        this%iter_pc, this%iter, this%res_norm
    string = trim(buffer)
  end function metrics_string

end module fgmres_solver_type
