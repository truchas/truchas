!!
!! GMRES_LEFT_SOLVER_TYPE
!!
!! This module defines a left-preconditioned GMRES solver type, and a
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

module gmres_left_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use vector_class
  use parameter_list_type
  use truchas_timers
  use parallel_communication, only: global_dot_product
  !use nlsol_type, only: nlsol_model
  use nlk_solver_type, only: nlk_solver_model
  implicit none
  private

  !public :: nlsol_model ! re-export

  type, public :: gmres_left_solver
    private
    class(nlk_solver_model), pointer :: model => null() ! unowned reference
    integer :: krylov_dim, iter_pc, max_iter, iter
    real(r8) :: res_norm, tol, rtol
    class(vector), allocatable :: r, w, v(:)
    real(r8), allocatable :: e1(:), h(:,:), work(:), udot(:)
  contains
    procedure :: init => gmres_left_solver_init
    !procedure :: setup => gmres_left_solver_setup
    procedure :: solve => gmres_left_solver_solve
    procedure :: metrics_string
  end type gmres_left_solver

contains

  subroutine gmres_left_solver_init(this, model, params, stat, errmsg)

    use parameter_list_type
    external :: dgels

    class(gmres_left_solver), intent(out) :: this
    class(nlk_solver_model), intent(in), target :: model
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context
    integer :: lwork

    stat = 0 ! TODO: error handling / argument verification
    this%model => model

    call params%get('gmres-krylov-dim', this%krylov_dim, default=5)
    call params%get('abs-tol', this%tol, default=1d-8)
    call params%get('rel-tol', this%rtol, default=0.0_r8)
    call params%get('max-iter', this%max_iter, default=20)

    call model%alloc_vector(this%r)
    call this%r%clone(this%w)
    call this%r%clone(this%v, this%krylov_dim+1)

    allocate(this%e1(this%krylov_dim+1), this%h(this%krylov_dim+1,this%krylov_dim))
    this%e1 = 0
    this%e1(1) = 1
    this%h = 0

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

  end subroutine gmres_left_solver_init


  subroutine gmres_left_solver_solve(this, u, errc)

    class(gmres_left_solver), intent(inout) :: this
    class(vector), intent(inout) :: u
    integer,  intent(out) :: errc

    integer :: iter, i, j
    real(r8) :: y(this%krylov_dim), anorm, rnorm, n2_b

    ASSERT(this%max_iter > 0)
    ASSERT(associated(this%model))

    call start_timer("GMRES solve")

    errc = 0
    this%iter_pc = 0
    this%h = 0
    call this%r%setval(0.0_r8)
    call this%model%compute_precon(u)
    n2_b = this%model%rhs%norm2()

    do iter = 1, this%max_iter+1
      !! compute the residual
      call this%model%compute_f(u, this%r)
      anorm = this%r%norm2()
      call this%model%apply_precon(u, this%r)
!ASSERT(all(ieee_is_finite(this%r)))

      rnorm = anorm / n2_b
      this%res_norm = this%r%norm2()
      print '("iter, res norm, anorm, rnorm, tol: ",i6,4es13.3)', &
          iter, this%res_norm, anorm, rnorm, this%tol
print *, "FOO", this%res_norm
      if (this%res_norm < this%tol .or. anorm < this%tol .or. rnorm < this%rtol) exit
      !this%v(:,1) = this%r / this%res_norm
      call this%v(1)%update(1.0_r8/this%res_norm, this%r, 0.0_r8)

      do j = 1, this%krylov_dim
        call this%model%compute_f(this%v(j), this%w, ax=.true.)
        call this%model%apply_precon(u, this%w)
        this%iter_pc = this%iter_pc + 1

        do i = 1, j
          this%h(i,j) = this%w%dot(this%v(i))
          !this%w = this%w - this%h(i,j) * this%v(:,i)
          call this%w%update(-this%h(i,j), this%v(i))
        end do
        this%h(j+1,j) = this%w%norm2()
        !this%v(:,j+1) = this%w / this%h(j+1,j)
        call this%v(j+1)%update(1.0_r8/this%h(j+1,j), this%w, 0.0_r8)

        !! TODO: van der Vorst has some extra logic for early exit in this loop.
        !! Should consider implementing this.
      end do

      this%e1 = 0
      this%e1(1) = this%res_norm
      call argmin(this%h, this%e1, this%work, y, errc)
      if (errc /= 0) exit
      !u = u + matmul(this%v(:,:this%krylov_dim), y)
      do i = 1, this%krylov_dim
        call u%update(y(i), this%v(i))
      end do
    end do

    this%iter = iter
    if (errc == 0 .and. iter > this%max_iter+1) errc = 1
    call stop_timer("GMRES solve")

  end subroutine gmres_left_solver_solve


  !! Return the y that minimizes ||b - A*y||_2.
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
    class(gmres_left_solver), intent(in) :: this
    character(:), allocatable :: string
    character(128) :: buffer
    write(buffer,'(i4," (PC), ",i4," (GMRES_LEFT), ",es14.4," (|r|)")') &
        this%iter_pc, this%iter, this%res_norm
    string = trim(buffer)
  end function metrics_string

end module gmres_left_solver_type
