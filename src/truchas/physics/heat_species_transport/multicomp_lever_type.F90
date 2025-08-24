#include "f90_assert.fpp"
!#define OLD_EUTECTIC

module multicomp_lever_type

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

  type, public :: multicomp_lever
    real(r8), allocatable :: liq_slope(:), part_coef(:)
    real(r8) :: T_f, Teut, eutectic_eps
    integer :: num_comp
    class(scalar_func), allocatable :: h_sol, h_liq
    class(scalar_func), allocatable :: h_sol_deriv, h_liq_deriv
    type(invf) :: my_invf
  contains
    procedure :: init
    procedure :: solve
    procedure :: solve_for_H_g
    procedure :: compute_f, compute_f_jac
    procedure :: compute_C_liq, compute_C_sol
    procedure, private :: temp
  end type

  type, public :: multicomp_lever_jac
    real(r8) :: dfgdg, dfgdH, dfgdT, dfHdg, dfHdH, dfHdT
  contains
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

    class(multicomp_lever), intent(out) :: this
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

    call params%get('eutectic-eps', this%eutectic_eps, stat, errmsg)
    if (stat /= 0) return

    !! Ridders solver to invert functions
    call this%my_invf%init(eps=1d-9,maxadj=10)

  end subroutine init

  subroutine get_liq_conc(this, c, g, c_liq)
    type(multicomp_lever), intent(in) :: this
    real(r8), intent(in)  :: c(:), g
    real(r8), intent(out) :: c_liq(:)
    ASSERT(size(c) == this%num_comp)
    ASSERT(size(c_liq) == this%num_comp)
    !ASSERT(g >= 0 .and. g <= 1)
    c_liq = c / ((1-g)*this%part_coef + g)
  end subroutine

  pure real(r8) function temp(this, g, c)
    class(multicomp_lever), intent(in) :: this
    real(r8), intent(in) :: g, c(:)
    integer :: k
    temp = this%T_f
    do k = 1, this%num_comp
      temp = temp + this%liq_slope(k)*c(k) / ((1-g)*this%part_coef(k) + g)
    end do
  end function

  pure real(r8) function temp_deriv(this, g, c)
    class(multicomp_lever), intent(in) :: this
    real(r8), intent(in) :: g, c(:)
    integer :: k
    temp_deriv = 0
    do k = 1, this%num_comp
      temp_deriv = temp_deriv + (this%part_coef(k)-1)*this%liq_slope(k)*c(k) &
                                  / ((1-g)*this%part_coef(k) + g)**2
    end do
  end function

  !! Given the enthalpy H and solute concentrations C, this solves for the
  !! consistent temperature T and liquid volume fraction g according to the
  !! lever rule. [T1,T2] is the initial bracket for finding T.

  subroutine solve(this, H, C, T1, T2, T, g)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: h, C(:), T1, T2
    real(r8), intent(out) :: g, T

    integer :: stat
    real(r8) :: Tmax, Hmax, Tmin, Hmin, dT, Heut, geut
    character(:), allocatable :: errmsg

    Tmax = this%temp(1.0_r8, C)
    Hmax = this%H_liq%eval([Tmax])
    if (H >= Hmax) then ! pure liquid
      !solve H = this%H_liq([T]) for T
      this%my_invf%func => H_liq_of_T
      call this%my_invf%compute(H, max(T1,Tmax), T2, T, stat, errmsg) !TODO: [T1,T2] bracket
      INSIST(stat == 0)
      g = 1
      return
    end if

    Tmin = this%temp(0.0_r8, C)
    if (Tmin < this%Teut) then ! eutectic transformation may be in play
      !solve this%Teut = this%temp(geut, C) for geut
      this%my_invf%func => T_of_g
      call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
      INSIST(stat == 0)
      dT = geut * this%eutectic_eps
      Tmin = this%Teut - dT
      Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
    end if

    Hmin = this%H_sol%eval([Tmin])
    if (H <= Hmin) then ! pure solid
      !solve H = this%H_sol([T]) for T
      this%my_invf%func => H_sol_of_T
      call this%my_invf%compute(H, T1, min(T2,Tmin), T, stat, errmsg) !TODO: [T1,T2] bracket
      INSIST(stat == 0)
      g = 0
      return
    end if

    if (H > merge(Heut, Hmin, Tmin < this%Teut)) then ! 2-phase region
      !solve H = (1-g)*this%H_sol%eval([this%temp(g,C)]) + g*this%H_liq%eval([this%temp(g,C)) for g
      this%my_invf%func => H_of_g
      call this%my_invf%compute(H, 0.0_r8, 1.0_r8, g, stat, errmsg)
      INSIST(stat == 0)
      T = this%temp(g, C)
    else ! smeared eutectic transformation
      !solve H = (1-g(T))*this%H_sol%eval([T]) + g(T)*this%H_liq%eval([T]) for T
      this%my_invf%func => H_of_T
      call this%my_invf%compute(H, Tmin, this%Teut, T, stat, errmsg)
      INSIST(stat == 0)
      g = (geut/dT)*(T - Tmin)
    end if

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C)
    end function

    real(r8) function H_sol_of_T(x) result(H)
      real(r8), intent(in) :: x
      H = this%H_sol%eval([x])
    end function

    real(r8) function H_liq_of_T(x) result(H)
      real(r8), intent(in) :: x
      H = this%H_liq%eval([x])
    end function

    real(r8) function H_of_g(x) result(H)
      real(r8), intent(in) :: x
      H = (1-x)*this%H_sol%eval([this%temp(x,C)]) + x*this%H_liq%eval([this%temp(x,C)])
    end function

    real(r8) function H_of_T(x) result(H)
      real(r8), intent(in) :: x
      real(r8) :: g
      g = (geut/dT)*(x - Tmin)
      H = (1-g)*this%H_sol%eval([x]) + g*this%H_liq%eval([x])
    end function

  end subroutine solve

  !! Given the temperature T and solute concentrations C, this solves for
  !! the consistent enthalpy H and liquid volume fraction according to the
  !! lever rule.

  subroutine solve_for_H_g(this, C, T, H, g)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: C(:,:), T(:)
    real(r8), intent(out) :: H(:), g(:)

    integer  :: j, stat
    real(r8) :: Tmin, Tmax, dT, geut
    character(:), allocatable :: errmsg

    ASSERT(size(H) == size(T))
    ASSERT(size(g) == size(T))

    do j = 1, size(T)
      Tmax = this%temp(1.0_r8, C(:,j))
      if (T(j) >= Tmax) then ! pure liquid
        g(j) = 1
        H(j) = this%H_liq%eval([T(j)])
        cycle
      end if

      Tmin = this%temp(0.0_r8, C(:,j))
      if (Tmin < this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, C) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        dT = geut * this%eutectic_eps
        Tmin = this%Teut - dT
      end if

      if (T(j) <= Tmin) then ! pure solid
        g(j) = 0
        H(j) = this%H_sol%eval([T(j)])
        cycle
      end if

      if (Tmin < this%Teut) then
        Tmin = this%Teut
      end if

      if (T(j) >= this%Teut) then ! 2-phase region
        ! solve T = this%temp(g, C(:,j)) for g
        this%my_invf%func => T_of_g
        call this%my_invf%compute(T(j), 0.0_r8, 1.0_r8, g(j), stat, errmsg)
        INSIST(stat == 0)
        H(j) = (1-g(j))*this%H_sol%eval([T(j)]) + g(j)*this%H_liq%eval([T(j)])
      else ! eutectic transformation
        g(j) = (geut/dT)*(T(j) - Tmin)
        H(j) = (1-g(j))*this%H_sol%eval([T(j)]) + g(j)*this%H_liq%eval([T(j)])
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C(:,j))
    end function

  end subroutine solve_for_H_g

  subroutine compute_C_liq(this, C, H, n, C_liq)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in) :: C(:,:), H(:)
    integer, intent(in) :: n
    real(r8), intent(out) :: C_liq(:)

    integer :: j, stat
    real(r8) :: Tmax, Hmax, Tmin, Hmin, g, geut, Heut, Ceut, dT
    character(:), allocatable :: errmsg

    do j = 1, size(C_liq)
      Tmax = this%temp(1.0_r8, C(:,j))
      Hmax = this%H_liq%eval([Tmax])
      if (H(j) >= Hmax) then ! pure liquid
        C_liq(j) = C(n,j)
        !g = 1
        cycle
      end if

      Tmin = this%temp(0.0_r8, C(:,j))
      if (Tmin < this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, C) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        dT = geut * this%eutectic_eps
        Tmin = this%Teut - dT
        Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
        Ceut = C(n,j) / ((1-geut)*this%part_coef(n) + geut)
      end if

      Hmin = this%H_sol%eval([Tmin])
      if (H(j) <= Hmin) then ! pure solid
        C_liq(j) = 0 ! actually undefined
        !g = 0
        cycle
      end if

      if (H(j) > merge(Heut, Hmin, Tmin < this%Teut)) then ! 2-phase region
        !solve H = (1-g)*this%H_sol%eval([this%temp(g,C)]) + g*this%H_liq%eval([this%temp(g,C)) for g
        this%my_invf%func => H_of_g
        call this%my_invf%compute(H(j), 0.0_r8, 1.0_r8, g, stat, errmsg)
        INSIST(stat == 0)
        C_liq(j) = C(n,j) / ((1-g)*this%part_coef(n) + g)
      else ! eutectic transformation
        C_liq(j) = Ceut
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C(:,j))
    end function

    real(r8) function H_of_g(x) result(H)
      real(r8), intent(in) :: x
      H = (1-x)*this%H_sol%eval([this%temp(x,C(:,j))]) + x*this%H_liq%eval([this%temp(x,C(:,j))])
    end function

  end subroutine compute_C_liq

  subroutine compute_C_sol(this, C, H, n, C_sol)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in) :: C(:,:), H(:)
    integer, intent(in) :: n
    real(r8), intent(out) :: C_sol(:)

    integer :: j, stat
    real(r8) :: Tmax, Hmax, Tmin, Hmin, g, geut, Heut, Ceut, T, dT
    character(:), allocatable :: errmsg

    do j = 1, size(C_sol)
      Tmax = this%temp(1.0_r8, C(:,j))
      Hmax = this%H_liq%eval([Tmax])
      if (H(j) >= Hmax) then ! pure liquid
        C_sol(j) = 0 ! actually undefined
        !g = 1
        cycle
      end if

      Tmin = this%temp(0.0_r8, C(:,j))
      if (Tmin < this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, C) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        dT = geut * this%eutectic_eps
        Tmin = this%Teut - dT
        Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
        Ceut = C(n,j) / ((1-geut)*this%part_coef(n) + geut)
      end if

      Hmin = this%H_sol%eval([Tmin])
      if (H(j) <= Hmin) then ! pure solid
        C_sol(j) = C(n,j)
        !g = 0
        cycle
      end if

      if (H(j) > merge(Heut, Hmin, Tmin < this%Teut)) then ! 2-phase region
        !solve H = (1-g)*this%H_sol%eval([this%temp(g,C)]) + g*this%H_liq%eval([this%temp(g,C)) for g
        this%my_invf%func => H_of_g
        call this%my_invf%compute(H(j), 0.0_r8, 1.0_r8, g, stat, errmsg)
        INSIST(stat == 0)
        C_sol(j) = this%part_coef(n)*C(n,j) / ((1-g)*this%part_coef(n) + g)
      else ! eutectic transformation
        !solve H = (1-g(T))*this%H_sol%eval([T]) + g(T)*this%H_liq%eval([T]) for T
        this%my_invf%func => H_of_T
        call this%my_invf%compute(H(j), Tmin, this%Teut, T, stat, errmsg)
        INSIST(stat == 0)
        g = (geut/dT)*(T - Tmin)
        C_sol(j) = (C(n,j) - g*Ceut) / (1 - g)
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C(:,j))
    end function

    real(r8) function H_of_g(x) result(H)
      real(r8), intent(in) :: x
      H = (1-x)*this%H_sol%eval([this%temp(x,C(:,j))]) + x*this%H_liq%eval([this%temp(x,C(:,j))])
    end function

    real(r8) function H_of_T(x) result(H)
      real(r8), intent(in) :: x
      real(r8) :: g
      g = (geut/dT)*(x - Tmin)
      H = (1-g)*this%H_sol%eval([x]) + g*this%H_liq%eval([x])
    end function

  end subroutine compute_C_sol

  !! Compute the residual of the algebraic liquid volume fraction relation
  !! given the liquid volume fraction g and enthalpy H.

  subroutine compute_f(this, C, g, H, T, fg, fH)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: C(:,:), g(:), H(:), T(:)
    real(r8), intent(out) :: fg(:), fH(:)

    integer :: j, stat
    character(:), allocatable :: errmsg
    real(r8) :: Tmax, Hmax, Tmin, Hmin, Heut, geut, TT, dT

    do j = 1, size(g)
      Tmax = this%temp(1.0_r8, C(:,j))
      Hmax = this%H_liq%eval([Tmax])
      if (H(j) >= Hmax) then ! pure liquid
        fg(j) = g(j) - 1
        cycle
      end if

      Tmin = this%temp(0.0_r8, C(:,j))
      if (Tmin < this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, C) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        dT = geut * this%eutectic_eps
        Tmin = this%Teut - dT
        Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
      end if

      Hmin = this%H_sol%eval([Tmin])
      if (H(j) <= Hmin) then ! pure solid
        fg(j) = g(j)
        cycle
      end if

      if (H(j) > merge(Heut, Hmin, Tmin < this%Teut)) then ! 2-phase region
        fg(j) = T(j) - this%temp(g(j),C(:,j))
      else ! eutectic transformation
        fg(j) = g(j) - (geut/dT)*(T(j) - Tmin)
      end if
    end do

    do j = 1, size(g)
      fH(j) = (1-g(j))*this%H_sol%eval([T(j)]) + g(j)*this%H_liq%eval([T(j)]) - H(j)
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C(:,j))
    end function

  end subroutine compute_f

  !! Compute the derivatives of the algebraic liquid volume fraction relation
  !! with respect to the liquid volume fraction g and enthalpy H.

  subroutine compute_f_jac(this, C, g, H, T, jac)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: C(:,:), g(:), H(:), T(:)
    type(multicomp_lever_jac), intent(out) :: jac(:)

    integer :: j, stat
    character(:), allocatable :: errmsg
    real(r8) :: Tmax, Hmax, Tmin, Hmin, Heut, geut, dT, TT

    do j = 1, size(g)

      Tmax = this%temp(1.0_r8, C(:,j))
      Hmax = this%H_liq%eval([Tmax])
      if (H(j) >= Hmax) then ! pure liquid
        jac(j)%dfgdg = 1
        jac(j)%dfgdH = 0
        jac(j)%dfgdT = 0
        cycle
      end if

      Tmin = this%temp(0.0_r8, C(:,j))
      if (Tmin < this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, C) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        dT = geut * this%eutectic_eps
        Tmin = this%Teut - dT
        Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
      end if

      Hmin = this%H_sol%eval([Tmin])
      if (H(j) <= Hmin) then ! pure solid
        jac(j)%dfgdg = 1
        jac(j)%dfgdH = 0
        jac(j)%dfgdT = 0
        cycle
      end if

      if (H(j) > merge(Heut, Hmin, Tmin < this%Teut)) then ! 2-phase region
        jac(j)%dfgdg = -temp_deriv(this, g(j), C(:,j))
        jac(j)%dfgdH = 0
        jac(j)%dfgdT = 1
      else ! eutectic transformation
        jac(j)%dfgdg = 1
        jac(j)%dfgdH = 0
        jac(j)%dfgdT = -geut/dT
      end if
    end do

    do j = 1, size(g)
      jac(j)%dfHdg = this%H_liq%eval([T(j)]) - this%H_sol%eval([T(j)])
      jac(j)%dfHdH = -1
      jac(j)%dfHdT = (1-g(j))*this%H_sol_deriv%eval([T(j)]) + g(j)*this%H_liq_deriv%eval([T(j)])
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C(:,j))
    end function

  end subroutine compute_f_jac

  subroutine lu_factor(this)
    class(multicomp_lever_jac), intent(inout) :: this
    this%dfHdg = this%dfHdg/this%dfgdg
    this%dfHdH = this%dfHdH - this%dfHdg*this%dfgdH
    this%dfHdT = this%dfHdT - this%dfHdg*this%dfgdT
  end subroutine

  subroutine lower_solve(this, fg, fH)
    class(multicomp_lever_jac), intent(in) :: this
    real(r8), intent(in) :: fg
    real(r8), intent(inout) :: fH
    fH = fH - this%dfHdg * fg
  end subroutine

  subroutine upper_solve(this, fg, fH, fT)
    class(multicomp_lever_jac), intent(in) :: this
    real(r8), intent(inout) :: fg, fH
    real(r8), intent(in) :: fT
    fH = (fH - this%dfHdT*fT) / this%dfHdH
    fg = (fg - this%dfgdH*fH - this%dfgdT*fT) / this%dfgdg
  end subroutine

end module
