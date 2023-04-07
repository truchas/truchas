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

module gmres_left_solver_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  implicit none
  private

  type, abstract, public :: gmres_left_solver
    private
    integer, public :: num_iter = 0
    real(r8), public :: r0norm, rel_rnorm
    integer :: krylov_dim, max_iter, print_level
    real(r8) :: atol, rtol
    integer, public :: iter, iter_pc
    real(r8), public :: res_norm
    class(vector), allocatable, public :: r
    class(vector), allocatable :: w, v(:)
    real(r8), allocatable :: e1(:), h(:,:), work(:)
  contains
    procedure :: init
    procedure :: solve
    procedure(matvec), deferred :: matvec
    procedure(apply_precon), deferred :: apply_precon
  end type

  abstract interface
    subroutine matvec(this, x, Ax)
      import gmres_left_solver, vector
      class(gmres_left_solver), intent(inout) :: this
      class(vector), intent(inout) :: x, Ax
    end subroutine
    subroutine apply_precon(this, x)
      import gmres_left_solver, vector
      class(gmres_left_solver), intent(inout) :: this
      class(vector), intent(inout) :: x
    end subroutine
  end interface

contains

  subroutine init(this, vec, params, stat, errmsg)

    use parameter_list_type
    external :: dgels

    class(gmres_left_solver), intent(out) :: this
    class(vector), intent(in) :: vec
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: lwork

    stat = 0 ! TODO: error handling / argument verification

    call params%get('krylov-dim', this%krylov_dim, default=5)
    call params%get('abs-tol', this%atol, default=1d-8)
    call params%get('rel-tol', this%rtol, default=0.0_r8)
    call params%get('max-iter', this%max_iter, default=100)
    call params%get('print-level', this%print_level, default=0)

    call vec%clone(this%r)
    call vec%clone(this%w)
    call vec%clone(this%v, this%krylov_dim+1)

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

  end subroutine init


  subroutine solve(this, b, x, stat)

    use truchas_logging_services

    class(gmres_left_solver), intent(inout) :: this
    class(vector), intent(in) :: b
    class(vector), intent(inout) :: x
    integer,  intent(out) :: stat

    integer :: iter, i, j
    real(r8) :: y(this%krylov_dim), rnorm
    character(80) :: message

    stat = 0
    this%h = 0

    !! Initial residual
    call this%matvec(x, this%r)
    call this%r%update(1.0_r8, b, -1.0_r8) ! r = b - r
    this%r0norm = this%r%norm2()

    !! Preconditioned initial residual
    call this%apply_precon(this%r)
    this%res_norm = this%r%norm2()
    if (this%print_level > 0) then
      write(message,'(i0,": |r|_2, |Pr|_2 =",2(es10.3,:,","))') 0, this%r0norm, this%res_norm
      call TLS_info(message)
    end if

    do iter = this%krylov_dim, this%max_iter, this%krylov_dim

      !this%v(:,1) = this%r / this%res_norm
      call this%v(1)%update(1.0_r8/this%res_norm, this%r, 0.0_r8)

      do j = 1, this%krylov_dim
        call this%matvec(this%v(j), this%w)
        call this%apply_precon(this%w)

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
      call argmin(this%h, this%e1, this%work, y, stat)
      if (stat /= 0) exit
      !x = x + matmul(this%v(:,:this%krylov_dim), y)
      do i = 1, this%krylov_dim
        call x%update(y(i), this%v(i))
      end do

      !! Compute the residual
      call this%matvec(x, this%r)
      call this%r%update(1.0_r8, b, -1.0_r8) ! r = b - r
      rnorm = this%r%norm2()

      !! Preconditioned residual
      call this%apply_precon(this%r)
      this%res_norm = this%r%norm2()

      if (this%print_level > 0) then
        write(message,'(t8,a,i4,*(a,es10.3))') 'gmres iterate', iter, &
            ': |r|=', rnorm, ', |Pr|=', this%res_norm
        call TLS_info(message)
      end if

      if (rnorm < max(this%atol, this%rtol * this%r0norm)) exit
      ! drop the norm on the preconditioned r in order for comparison to NLK
    end do

    this%num_iter = iter
    this%rel_rnorm = rnorm/this%r0norm
    if (stat == 0 .and. iter > this%max_iter) stat = 1

  end subroutine solve


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

end module gmres_left_solver_class
