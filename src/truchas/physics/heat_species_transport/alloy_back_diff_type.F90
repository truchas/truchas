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
    real(r8), allocatable :: liq_slope(:), part_coef(:)
    real(r8) :: T_f, Teut, dHdTmax, gamma, Hdot
    integer :: num_comp
    class(scalar_func), allocatable :: h_sol, h_liq
    class(scalar_func), allocatable :: h_sol_deriv, h_liq_deriv
    type(invf) :: my_invf
  contains
    procedure :: init
    procedure :: H_of_g_T, T_of_g_H
    procedure :: compute_f, compute_f_jac
    procedure, private :: T_sol, T_liq
  end type

  type, public :: alloy_back_diff_jac
    real(r8), allocatable :: dfpdp(:), dfpdg(:), dfgdp(:)
    real(r8) :: dfgdg, dfgdT, dfHdg, dfHdH, dfHdT
  contains
    procedure :: init => jac_init
    procedure :: lu_factor, lower_solve, upper_solve
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

    call params%get('dhdtmax', this%dHdTmax, stat, errmsg, default = 100.0_r8)
    if (stat /= 0) return

    call params%get('gamma', this%gamma, stat, errmsg)
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

  function H_of_g_T(this, g, T) result(H)
    class(alloy_back_diff), intent(in) :: this
    real(r8), intent(in) :: g, T
    real(r8) :: H
    H = (1-g)*this%H_sol%eval([T]) + g*this%H_liq%eval([T])
  end function

  function T_of_g_H(this, g, H, T1, T2) result(T)
    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in) :: g, H, T1, T2
    real(r8) :: T
    integer :: stat
    character(:), allocatable :: errmsg
    !solve H = (1-g)*this%H_sol([T]) + g*this%H_liq([T]) for T
    this%my_invf%func => H_of_T
    call this%my_invf%compute(H, T1, T2, T, stat, errmsg) !TODO: [T1,T2] bracket
    INSIST(stat == 0)
  contains
    real(r8) function H_of_T(x) result(H)
      real(r8), intent(in) :: x
      H = (1-g)*this%H_sol%eval([x]) + g*this%H_liq%eval([x])
    end function
  end function

  subroutine compute_f(this, C, Cdot, p, g, H, T, pdot, gdot, Hdot, fp, fg, fH)

    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in)  :: C(:), Cdot(:) ! parameters
    real(r8), intent(in)  :: p(:), g, H, T, pdot(:), gdot, Hdot
    real(r8), intent(out) :: fp(:), fg, fH

    real(r8) :: Tmax, Hmax

    Tmax = this%T_liq(C)
    Hmax = this%H_liq%eval([Tmax])
    if (H >= Hmax) then
      fp = p - C
      fg = g - 1
      fH = this%H_liq%eval([T]) - H
      return
    end if

    !Tmin = this%T_sol(this%C)
    !Hmin = this%H_sol%eval([Tmin])
    !if (H <= Hmin) then
    if (g <= 0) then
      fp = p
      fg = g
      fH = this%H_sol%eval([T]) - H
      return
    end if

    fp = g*(pdot - Cdot) - this%part_coef*p*gdot + this%gamma*(this%part_coef*p - (g/(1-g+1d-6))*(C - p))
    fg = g*(T-this%T_f) - dot_product(this%liq_slope, p)
    fH = (1-g)*this%H_sol%eval([T]) + g*this%H_liq%eval([T]) - H

  end subroutine compute_f

  subroutine compute_f_jac(this, C, Cdot, p, g, H, T, pdot, gdot, dt, jac)

    class(alloy_back_diff), intent(inout) :: this
    real(r8), intent(in) :: C(:), Cdot(:) ! parameters
    real(r8), intent(in) :: p(:), g, H, T, pdot(:), gdot, dt
    type(alloy_back_diff_jac), intent(inout) :: jac

    real(r8) :: Tmax, Hmax

    Tmax = this%T_liq(C)
    Hmax = this%H_liq%eval([Tmax])
    if (H >= Hmax) then
      jac%dfpdp = 1
      jac%dfpdg = 0
      jac%dfgdp = 0
      jac%dfgdg = 1
      jac%dfgdT = 0
      jac%dfHdg = 0
      jac%dfHdH = -1
      jac%dfHdT = this%H_liq_deriv%eval([T])
      return
    end if

    !Tmin = this%T_sol(C)
    !Hmin = this%H_sol%eval([Tmin])
    !if (H <= Hmin) then
    if (g <= 0) then
      jac%dfpdp = 1
      jac%dfpdg = 0
      jac%dfgdp = 0
      jac%dfgdg = 1
      jac%dfgdT = 0
      jac%dfHdg = 0
      jac%dfHdH = -1
      jac%dfHdT = this%H_sol_deriv%eval([T])
      return
    end if

    jac%dfpdp = g/dt - this%part_coef*gdot + this%gamma*(this%part_coef + g/(1-g+1d-6))
    jac%dfpdg = (pdot - Cdot) - this%part_coef*p/dt - this%gamma*(C - p)/(1-g+1d-6)**2

    jac%dfgdp = -this%liq_slope
    jac%dfgdg = T - this%T_f
    jac%dfgdT = g

    jac%dfHdg = this%H_liq%eval([T]) - this%H_sol%eval([T])
    jac%dfHdH = -1
    jac%dfHdT = (1-g)*this%H_sol_deriv%eval([T]) + g*this%H_liq_deriv%eval([T])

  end subroutine compute_f_jac

  subroutine jac_init(this, n)
    class(alloy_back_diff_jac), intent(out) :: this
    integer, intent(in) :: n
    allocate(this%dfpdp(n), this%dfpdg(n), this%dfgdp(n))
  end subroutine

  subroutine lu_factor(this)
    class(alloy_back_diff_jac), intent(inout) :: this
    this%dfgdp = this%dfgdp/this%dfpdp
    this%dfgdg = this%dfgdg - dot_product(this%dfgdp, this%dfpdg)
    this%dfHdg = this%dfHdg/this%dfgdg
    this%dfHdT = this%dfHdT - this%dfHdg*this%dfgdT
  end subroutine

  subroutine lower_solve(this, fp, fg, fH)
    class(alloy_back_diff_jac), intent(in) :: this
    real(r8), intent(in) :: fp(:)
    real(r8), intent(inout) :: fg, fH
    fg = fg - dot_product(this%dfgdp, fp)
    fH = fH - this%dfHdg * fg
  end subroutine

  subroutine upper_solve(this, fp, fg, fH, fT)
    class(alloy_back_diff_jac), intent(in) :: this
    real(r8), intent(inout) :: fp(:), fg, fH
    real(r8), intent(in) :: fT
    fH = (fH - this%dfHdT*fT) / this%dfHdH
    fg = (fg - this%dfgdT*fT) / this%dfgdg
    fp = (fp - this%dfpdg*fg) / this%dfpdp
  end subroutine

end module
