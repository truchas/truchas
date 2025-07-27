#include "f90_assert.fpp"

module alloy_back_diff_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  use inverse_func_class
  implicit none
  private

  type, extends(inverse_func) :: invf
    procedure(func), nopass, pointer :: func => null()
  contains
    procedure :: g => invf_g
  end type

  abstract interface
    function func(x) result (gx)
      import r8
      real(r8), intent(in) :: x
      real(r8) :: gx
    end function
  end interface

  type, public :: alloy_back_diff
    real(r8), allocatable :: liq_slope(:), part_coef(:), C0(:)
    real(r8) :: T_f, Teut, dHdTmax, gamma, Hdot
    integer :: num_comp
    class(scalar_func), allocatable :: h_sol, h_liq
    class(scalar_func), allocatable :: h_sol_deriv, h_liq_deriv
    type(invf) :: my_invf
  contains
    procedure :: init
    procedure :: compute_f1, compute_f1_jac
    procedure :: compute_f2, compute_f2_jac
    procedure :: compute_f3, compute_f3_jac
    procedure, private :: T_sol, T_liq
  end type

contains

  real(r8) function invf_g(this, x) result(gx)
    class(invf), intent(in) :: this
    real(r8), intent(in) :: x
    gx = this%func(x)
  end function

  subroutine init(this, matl, params, stat, errmsg)

    use material_class
    use parameter_list_type

    class(alloy_back_diff), intent(out) :: this
    class(material), intent(in) :: matl
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: k
    type(phase), pointer :: phi

    if (matl%num_phase() == 2) then
      phi => matl%phase_ref(1)
      call phi%get_prop('enthalpy', this%h_sol)
      phi => matl%phase_ref(2)
      call phi%get_prop('enthalpy', this%h_liq)
    else
      stat = 1
      errmsg = 'not a 2-phase material'
      return
    end if

    block
      use scalar_func_tools
      call alloc_scalar_func_deriv(this%h_sol, 1, this%h_sol_deriv)
      call alloc_scalar_func_deriv(this%h_liq, 1, this%h_liq_deriv)
    end block

    call params%get('num-comp', this%num_comp, stat, errmsg, default=1)
    if (stat /= 0) return

    !! Fusion temperature of the pure solvent
    call params%get('temp-fusion', this%T_f, stat, errmsg)
    if (stat /= 0) return

    !! Eutectic temperature (or 0)
    call params%get('temp-eutectic', this%Teut, stat, errmsg, default=-huge(1.0_r8))
    if (stat /= 0) return
    if (this%Teut >= this%T_f) then
      stat = 1
      errmsg = 'temp-eutectic is >= temp-fusion'
      return
    end if

    !! Component liquidus slopes
    call params%get('liq-slope', this%liq_slope, stat, errmsg)
    if (stat /= 0) return
    if (size(this%liq_slope) /= this%num_comp) then
      stat = 1
      errmsg = 'invalid liq-slope size'
      return
    else if (any(this%liq_slope == 0)) then
      stat = 1
      errmsg = 'liq-slope must not be 0.0'
      return
    end if

    !! Component partition coefficients
    call params%get('part-coef', this%part_coef, stat, errmsg)
    if (stat /= 0) return
    if (size(this%part_coef) /= this%num_comp) then
      stat = 1
      errmsg = 'invalid part-coef size'
      return
    else if (any(this%part_coef <= 0)) then
      stat = 1
      errmsg = 'part-coef must be > 0.0'
    end if

    !! Check compatibility
    do k = 1, this%num_comp
      if (this%liq_slope(k)*(1-this%part_coef(k)) >= 0) then
        stat = 1
        if (this%liq_slope(k) > 0) then
          errmsg = 'part-coef must be > 1.0 when liq-slope is > 0.0'
        else
          errmsg = 'part-coef must be < 1.0 when liq-slope is < 0.0'
        end if
        return
      end if
    end do

    !! Solute concentrations (uniform, no advection); TODO: spatially variable
    call params%get('concentration', this%C0, stat, errmsg)
    if (stat /= 0) return
    if (size(this%C0) /= this%num_comp) then
      stat = 1
      errmsg = 'invalid concentration size'
      return
    else if (any(this%C0 < 0) .or. sum(this%C0) >= 1) then
      stat = 1
      errmsg = 'invalid concentration values'
      return
    end if

    call params%get('dhdtmax', this%dHdTmax, stat, errmsg, default = 100.0_r8)
    if (stat /= 0) return

    call params%get('gamma', this%gamma, stat, errmsg)
    if (stat /= 0) return

    call params%get('Hdot', this%Hdot, stat, errmsg)
    if (stat /= 0) return

    !! Ridders solver to invert functions
    call this%my_invf%init(eps=1d-9,maxadj=10)

  end subroutine init

  pure real(r8) function T_liq(this, C)
    class(alloy_back_diff), intent(in) :: this
    real(r8), intent(in) :: C(:)
    integer :: k
    T_liq = this%T_f
    do k = 1, this%num_comp
      T_liq = T_liq + this%liq_slope(k)*C(k)
    end do
  end function

  pure real(r8) function T_sol(this, C)
    class(alloy_back_diff), intent(in) :: this
    real(r8), intent(in) :: C(:)
    integer :: k
    T_sol = this%T_f
    do k = 1, this%num_comp
      T_sol = T_sol + this%liq_slope(k)*C(k)/this%part_coef(k)
    end do
  end function

  subroutine compute_f1(this, u, udot, f)

    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in)  :: u(:), udot(:)
    real(r8), intent(out) :: f(:)

    integer  :: n
    real(r8) :: Tmin, Tmax, Hmin, Hmax

    n = this%num_comp
    associate (v => u(1:n), g => u(n+1), H => u(n+2), T => u(n+3), &
               vdot => udot(1:n), gdot => udot(n+1), Hdot => udot(n+2))

      Tmax = this%T_liq(this%C0)
      Hmax = this%H_liq%eval([Tmax])
      if (H >= Hmax) then
        f(1:n) = v - this%C0
        f(n+1) = g - 1
        f(n+2) = this%H_liq%eval([T]) - H
        f(n+3) = Hdot - this%Hdot
        return
      end if

      !Tmin = this%T_sol(this%C0)
      !Hmin = this%H_sol%eval([Tmin])
      !if (H <= Hmin) then
      if (g <= 0) then
        f(1:n) = v
        f(n+1) = g
        f(n+2) = this%H_sol%eval([T]) - H
        f(n+3) = Hdot - this%Hdot
        return
      end if

      f(1:n) = vdot - this%part_coef*(v/g)*gdot + &
               this%gamma*(this%part_coef*(v/g) - (this%C0 - v)/(1-g + 1.0d-10))
      f(n+1) = (1-g)*this%H_sol%eval([this%T_liq(v/g)]) + g*this%H_liq%eval([this%T_liq(v/g)]) - H
      f(n+2) = (1-g)*this%H_sol%eval([T]) + g*this%H_liq%eval([T]) - H
      f(n+3) = Hdot - this%Hdot

    end associate

  end subroutine compute_f1

  subroutine compute_f1_jac(this, u, udot, dt, jac)

    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in)  :: u(:), udot(:), dt
    real(r8), intent(out) :: jac(:,:)

    integer  :: n, k
    real(r8) :: Tmin, Tmax, Hmin, Hmax

    n = this%num_comp
    associate (v => u(1:n), g => u(n+1), H => u(n+2), T => u(n+3), gdot => udot(n+1))

      jac = 0

      Tmax = this%T_liq(this%C0)
      Hmax = this%H_liq%eval([Tmax])
      if (H >= Hmax) then
        do k = 1, n
          jac(k,k) = 1
        end do
        jac(n+1,n+1) = 1
        jac(n+2,n+2) = -1
        jac(n+2,n+3) = this%H_liq_deriv%eval([T])
        jac(n+3,n+2) = 1/dt
        return
      end if

      !Tmin = this%T_sol(this%C0)
      !Hmin = this%H_sol%eval([Tmin])
      !if (H <= Hmin) then
      if (g <= 0) then
        do k = 1, n
          jac(k,k) = 1
        end do
        jac(n+1,n+1) = 1
        jac(n+2,n+2) = -1
        jac(n+2,n+3) = this%H_sol_deriv%eval([T])
        jac(n+3,n+2) = 1/dt
        return
      end if

      do k = 1, n
        jac(k,k) = 1/dt - this%part_coef(k)*(1/g)*gdot + &
                   this%gamma*(this%part_coef(k)*(1/g) + 1/(1-g+1d-10))
      end do
      jac(:n,n+1) = - this%part_coef*(v/g)/dt - this%part_coef*gdot*(-v/g**2) &
                    + this%gamma*(this%part_coef*(-v/g**2) - (this%C0 - v)/(1-g+1d-10)**2)

      jac(n+1,:n) = ((1-g)*this%H_sol_deriv%eval([this%T_liq(v/g)]) + &
                         g*this%H_liq_deriv%eval([this%T_liq(v/g)])) * (this%liq_slope/g)
      jac(n+1,n+1) = ((1-g)*this%H_sol_deriv%eval([this%T_liq(v/g)]) + &
                          g*this%H_liq_deriv%eval([this%T_liq(v/g)])) * &
                          (-dot_product(this%liq_slope, v)/g**2) &
                     + (this%H_liq%eval([this%T_liq(v/g)]) - this%H_sol%eval([this%T_liq(v/g)]))
      jac(n+1,n+2) = -1

      jac(n+2,n+1) = this%H_liq%eval([T]) - this%H_sol%eval([T])
      jac(n+2,n+2) = -1
      jac(n+2,n+3) = (1-g)*this%H_sol_deriv%eval([T]) + g*this%H_liq_deriv%eval([T])

      jac(n+3,n+2) = 1/dt

    end associate

  end subroutine compute_f1_jac

  subroutine compute_f2(this, u, udot, f)

    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in)  :: u(:), udot(:)
    real(r8), intent(out) :: f(:)

    integer  :: n
    real(r8) :: Tmin, Tmax, Hmin, Hmax

    n = this%num_comp
    associate (v => u(1:n), g => u(n+1), H => u(n+2), T => u(n+3), &
               vdot => udot(1:n), gdot => udot(n+1), Hdot => udot(n+2))

      Tmax = this%T_liq(this%C0)
      Hmax = this%H_liq%eval([Tmax])
      if (H >= Hmax) then
        f(1:n) = v - this%C0
        f(n+1) = g - 1
        f(n+2) = this%H_liq%eval([T]) - H
        f(n+3) = Hdot - this%Hdot
        return
      end if

      !Tmin = this%T_sol(this%C0)
      !Hmin = this%H_sol%eval([Tmin])
      !if (H <= Hmin) then
      if (g <= 0) then
        f(1:n) = v
        f(n+1) = g
        f(n+2) = this%H_sol%eval([T]) - H
        f(n+3) = Hdot - this%Hdot
        return
      end if

      f(1:n) = vdot - this%part_coef*(v/g)*gdot + &
               this%gamma*(this%part_coef*(v/g) - (this%C0 - v)/(1-g + 1.0d-10))
      !f(n+1) = (1-g)*this%H_sol%eval([this%T_liq(v/g)]) + g*this%H_liq%eval([this%T_liq(v/g)]) - H
      f(n+1) = g*(T-this%T_f) - dot_product(this%liq_slope, v)
      f(n+2) = (1-g)*this%H_sol%eval([T]) + g*this%H_liq%eval([T]) - H
      f(n+3) = Hdot - this%Hdot

    end associate

  end subroutine compute_f2

  subroutine compute_f2_jac(this, u, udot, dt, jac)

    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in)  :: u(:), udot(:), dt
    real(r8), intent(out) :: jac(:,:)

    integer  :: n, k
    real(r8) :: Tmin, Tmax, Hmin, Hmax

    n = this%num_comp
    associate (v => u(1:n), g => u(n+1), H => u(n+2), T => u(n+3), gdot => udot(n+1))

      jac = 0

      Tmax = this%T_liq(this%C0)
      Hmax = this%H_liq%eval([Tmax])
      if (H >= Hmax) then
        do k = 1, n
          jac(k,k) = 1
        end do
        jac(n+1,n+1) = 1
        jac(n+2,n+2) = -1
        jac(n+2,n+3) = this%H_liq_deriv%eval([T])
        jac(n+3,n+2) = 1/dt
        return
      end if

      !Tmin = this%T_sol(this%C0)
      !Hmin = this%H_sol%eval([Tmin])
      !if (H <= Hmin) then
      if (g <= 0) then
        do k = 1, n
          jac(k,k) = 1
        end do
        jac(n+1,n+1) = 1
        jac(n+2,n+2) = -1
        jac(n+2,n+3) = this%H_sol_deriv%eval([T])
        jac(n+3,n+2) = 1/dt
        return
      end if

      do k = 1, n
        jac(k,k) = 1/dt - this%part_coef(k)*(1/g)*gdot + &
                   this%gamma*(this%part_coef(k)*(1/g) + 1/(1-g+1d-10))
      end do
      jac(:n,n+1) = - this%part_coef*(v/g)/dt - this%part_coef*gdot*(-v/g**2) &
                    + this%gamma*(this%part_coef*(-v/g**2) - (this%C0 - v)/(1-g+1d-10)**2)

      !jac(n+1,:n) = ((1-g)*this%H_sol_deriv%eval([this%T_liq(v/g)]) + &
      !                   g*this%H_liq_deriv%eval([this%T_liq(v/g)])) * (this%liq_slope/g)
      !jac(n+1,n+1) = ((1-g)*this%H_sol_deriv%eval([this%T_liq(v/g)]) + &
      !                    g*this%H_liq_deriv%eval([this%T_liq(v/g)])) * &
      !                    (-dot_product(this%liq_slope, v)/g**2) &
      !               + (this%H_liq%eval([this%T_liq(v/g)]) - this%H_sol%eval([this%T_liq(v/g)]))
      !jac(n+1,n+2) = -1

      jac(n+1,:n) = -this%liq_slope
      jac(n+1,n+1) = T - this%T_f
      jac(n+1,n+3) = g

      jac(n+2,n+1) = this%H_liq%eval([T]) - this%H_sol%eval([T])
      jac(n+2,n+2) = -1
      jac(n+2,n+3) = (1-g)*this%H_sol_deriv%eval([T]) + g*this%H_liq_deriv%eval([T])

      jac(n+3,n+2) = 1/dt

    end associate

  end subroutine compute_f2_jac

  subroutine compute_f3(this, u, udot, f)

    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in)  :: u(:), udot(:)
    real(r8), intent(out) :: f(:)

    integer  :: n
    real(r8) :: Tmin, Tmax, Hmin, Hmax

    n = this%num_comp
    associate (v => u(1:n), g => u(n+1), H => u(n+2), T => u(n+3), &
               vdot => udot(1:n), gdot => udot(n+1), Hdot => udot(n+2))

      Tmax = this%T_liq(this%C0)
      Hmax = this%H_liq%eval([Tmax])
      if (H >= Hmax) then
        f(1:n) = v - this%C0
        f(n+1) = g - 1
        f(n+2) = this%H_liq%eval([T]) - H
        f(n+3) = Hdot - this%Hdot
        return
      end if

      !Tmin = this%T_sol(this%C0)
      !Hmin = this%H_sol%eval([Tmin])
      !if (H <= Hmin) then
      if (g <= 0) then
        f(1:n) = v
        f(n+1) = g
        f(n+2) = this%H_sol%eval([T]) - H
        f(n+3) = Hdot - this%Hdot
        return
      end if

      f(1:n) = g*vdot - this%part_coef*v*gdot + &
               this%gamma*(this%part_coef*v - (g/(1-g+1d-10))*(this%C0 - v))
      !f(n+1) = (1-g)*this%H_sol%eval([this%T_liq(v/g)]) + g*this%H_liq%eval([this%T_liq(v/g)]) - H
      f(n+1) = g*(T-this%T_f) - dot_product(this%liq_slope, v)
      f(n+2) = (1-g)*this%H_sol%eval([T]) + g*this%H_liq%eval([T]) - H
      f(n+3) = Hdot - this%Hdot

    end associate

  end subroutine compute_f3

  subroutine compute_f3_jac(this, u, udot, dt, jac)

    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in)  :: u(:), udot(:), dt
    real(r8), intent(out) :: jac(:,:)

    integer  :: n, k
    real(r8) :: Tmin, Tmax, Hmin, Hmax

    n = this%num_comp
    associate (v => u(1:n), g => u(n+1), H => u(n+2), T => u(n+3), vdot => udot(1:n), gdot => udot(n+1))

      jac = 0

      Tmax = this%T_liq(this%C0)
      Hmax = this%H_liq%eval([Tmax])
      if (H >= Hmax) then
        do k = 1, n
          jac(k,k) = 1
        end do
        jac(n+1,n+1) = 1
        jac(n+2,n+2) = -1
        jac(n+2,n+3) = this%H_liq_deriv%eval([T])
        jac(n+3,n+2) = 1/dt
        return
      end if

      !Tmin = this%T_sol(this%C0)
      !Hmin = this%H_sol%eval([Tmin])
      !if (H <= Hmin) then
      if (g <= 0) then
        do k = 1, n
          jac(k,k) = 1
        end do
        jac(n+1,n+1) = 1
        jac(n+2,n+2) = -1
        jac(n+2,n+3) = this%H_sol_deriv%eval([T])
        jac(n+3,n+2) = 1/dt
        return
      end if

      do k = 1, n
        jac(k,k) = g/dt - this%part_coef(k)*gdot + this%gamma*(this%part_coef(k) + g/(1-g+1d-10))
      end do
      jac(:n,n+1) = vdot - this%part_coef*v/dt - this%gamma*(this%C0 - v)/(1-g+1d-10)**2

      !jac(n+1,:n) = ((1-g)*this%H_sol_deriv%eval([this%T_liq(v/g)]) + &
      !                   g*this%H_liq_deriv%eval([this%T_liq(v/g)])) * (this%liq_slope/g)
      !jac(n+1,n+1) = ((1-g)*this%H_sol_deriv%eval([this%T_liq(v/g)]) + &
      !                    g*this%H_liq_deriv%eval([this%T_liq(v/g)])) * &
      !                    (-dot_product(this%liq_slope, v)/g**2) &
      !               + (this%H_liq%eval([this%T_liq(v/g)]) - this%H_sol%eval([this%T_liq(v/g)]))
      !jac(n+1,n+2) = -1

      jac(n+1,:n) = -this%liq_slope
      jac(n+1,n+1) = T - this%T_f
      jac(n+1,n+3) = g

      jac(n+2,n+1) = this%H_liq%eval([T]) - this%H_sol%eval([T])
      jac(n+2,n+2) = -1
      jac(n+2,n+3) = (1-g)*this%H_sol_deriv%eval([T]) + g*this%H_liq_deriv%eval([T])

      jac(n+3,n+2) = 1/dt

    end associate

  end subroutine compute_f3_jac

!
!  subroutine compute_f1(this, u, udot, f)
!
!    class(alloy_back_diff), intent(inout) :: this
!    real(r8), intent(in)  :: u(:), udot(:)
!    real(r8), intent(out) :: f(:)
!
!    integer  :: n
!    real(r8) :: Tmin, Tmax, Hmin, Hmax
!
!    n = this%num_comp
!    associate (v => u(1:n), g => u(n+1), H => u(n+2), T => u(n+3), &
!               vdot => udot(1:n), gdot => udot(n+1), Hdot => udot(n+2))
!
!      Tmax = this%T_liq(this%C0)
!      Hmax = this%H_liq%eval([Tmax])
!      if (H >= Hmax) then
!        f(1:n) = v - this%C0
!        f(n+1) = g - 1
!        f(n+2) = this%H_liq%eval([T]) - H
!        f(n+3) = Hdot - this%Hdot
!        return
!      end if
!
!      Tmin = this%T_sol(this%C0)
!      Hmin = this%H_sol%eval([Tmin])
!      if (H <= Hmin) then
!        f(1:n) = v
!        f(n+1) = g
!        f(n+2) = this%H_sol%eval([T]) - H
!        f(n+3) = Hdot - this%Hdot
!        return
!      end if
!
!      f(1:n) = vdot - this%part_coef*(v/g)*gdot + &
!               this%gamma*(this%part_coef*(v/g) - (this%C0 - v)/(1-g + 1.0d-10))
!      f(n+1) = (1-g)*this%H_sol%eval([this%T_liq(v/g)]) + g*this%H_liq%eval([this%T_liq(v/g)]) - H
!      f(n+2) = (1-g)*this%H_sol%eval([T]) + g*this%H_liq%eval([T]) - H
!      f(n+3) = Hdot - this%Hdot
!
!    end associate
!
!  end subroutine compute_f1
!
!  subroutine compute_f1_jac(this, u, udot, dt, jac)
!
!    class(alloy_back_diff), intent(inout) :: this
!    real(r8), intent(in)  :: u(:), udot(:), dt
!    real(r8), intent(out) :: jac(:,:)
!
!    integer  :: n, k
!    real(r8) :: Tmin, Tmax, Hmin, Hmax
!
!    n = this%num_comp
!    associate (v => u(1:n), g => u(n+1), H => u(n+2), T => u(n+3), gdot => udot(n+1))
!
!      jac = 0
!
!      Tmax = this%T_liq(this%C0)
!      Hmax = this%H_liq%eval([Tmax])
!      if (H >= Hmax) then
!        do k = 1, n
!          jac(k,k) = 1
!        end do
!        jac(n+1,n+1) = 1
!        jac(n+2,n+2) = -1
!        jac(n+2,n+3) = this%H_liq_deriv%eval([T])
!        jac(n+3,n+2) = 1/dt
!        return
!      end if
!
!      Tmin = this%T_sol(this%C0)
!      Hmin = this%H_sol%eval([Tmin])
!      if (H <= Hmin) then
!        do k = 1, n
!          jac(k,k) = 1
!        end do
!        jac(n+1,n+1) = 1
!        jac(n+2,n+2) = -1
!        jac(n+2,n+3) = this%H_sol_deriv%eval([T])
!        jac(n+3,n+2) = 1/dt
!        return
!      end if
!
!      do k = 1, n
!        jac(k,k) = 1/dt - this%part_coef(k)*(1/g)*gdot + &
!                   this%gamma*(this%part_coef(k)*(1/g) + 1/(1-g+1d-10))
!      end do
!      jac(:n,n+1) = - this%part_coef*(v/g)/dt - this%part_coef*gdot*(-v/g**2) &
!                    + this%gamma*(this%part_coef*(-v/g**2) - (this%C0 - v)/(1-g+1d-10)**2)
!
!      jac(n+1,:n) = ((1-g)*this%H_sol_deriv%eval([this%T_liq(v/g)]) + &
!                         g*this%H_liq_deriv%eval([this%T_liq(v/g)])) * (this%liq_slope/g)
!      jac(n+1,n+1) = ((1-g)*this%H_sol_deriv%eval([this%T_liq(v/g)]) + &
!                          g*this%H_liq_deriv%eval([this%T_liq(v/g)])) * &
!                          (-dot_product(this%liq_slope, v)/g**2) &
!                     + (this%H_liq%eval([this%T_liq(v/g)]) - this%H_sol%eval([this%T_liq(v/g)]))
!      jac(n+1,n+2) = -1
!
!      jac(n+2,n+1) = this%H_liq%eval([T]) - this%H_sol%eval([T])
!      jac(n+2,n+2) = -1
!      jac(n+2,n+3) = (1-g)*this%H_sol_deriv%eval([T]) + g*this%H_liq_deriv%eval([T])
!
!      jac(n+3,n+2) = 1/dt
!
!    end associate
!
!  end subroutine compute_f1_jac

end module
