!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module smold_nlsol_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use nlsol_type
  use var_vector_types, only: real_var_vector
  implicit none
  private

  type, extends(nlsol_model), public :: smold_nlsol_model
    private
    procedure(residual), nopass, pointer :: residualf => null()

    ! preconditioner data & parameters
    type(real_var_vector), pointer :: precon_a(:) => null() ! unowned reference
    real(r8) :: omega = 1.0_r8
    integer :: niter = 1
  contains
    procedure :: init
    !! Deferred procedures from nlsol_model
    procedure :: size => model_size
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
  end type smold_nlsol_model

  abstract interface
    subroutine residual(x_old, x, r)
      import r8
      real(r8), intent(in)  :: x_old(:), x(:)
      real(r8), intent(out) :: r(:)
    end subroutine
  end interface

contains

  subroutine init(this, residualf, precon_a, stat, errmsg)

    class(smold_nlsol_model), intent(out) :: this
    procedure(residual) :: residualf
    type(real_var_vector), intent(in), target :: precon_a(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context

    this%residualf => residualf
    this%precon_a => precon_a

  end subroutine init

  integer function model_size(this)
    class(smold_nlsol_model), intent(in) :: this
    model_size = size(this%precon_a)
  end function model_size

  subroutine compute_f(this, t, u, udot, f)
    class(smold_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:), udot(:)
    real(r8), intent(out), contiguous :: f(:)
    real(r8) :: u_old(size(u))
    u_old = u - udot ! h = 1 was passed to the solver
    call this%residualf(u_old, u, f)
  end subroutine

  ! Diagonal scaling
  subroutine apply_precon(this, t, u, f)

    class(smold_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous, target :: u(:)
    real(r8), intent(inout), contiguous :: f(:)

    integer :: i, j
    real(r8) :: d

    do j = 1, size(this%precon_a)
      d = f(j) / this%precon_a(j)%v(1)
      do i = 1, this%niter
        f(j) = f(j) + this%omega*(d-f(j))
      end do
    end do

  end subroutine

  subroutine compute_precon(this, t, u, dt)
    class(smold_nlsol_model) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), contiguous, target :: u(:)
    ! no op
  end subroutine

  real(r8) function du_norm(this, t, u, du)
    use parallel_communication, only: global_maxval
    class(smold_nlsol_model) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), contiguous :: u(:), du(:)
    real(r8) :: l
    du_norm = 1
    l = global_maxval(abs(u))
    if (l > 0) du_norm = global_maxval(abs(du)) / l
    ! du_norm = maxval(abs(du) / (this%atol + this%rtol*abs(u)))
    ! du_norm = global_maxval(du_norm)
  end function

end module smold_nlsol_model_type
