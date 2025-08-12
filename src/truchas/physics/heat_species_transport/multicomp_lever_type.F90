#include "f90_assert.fpp"

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
    real(r8) :: T_f, Teut, dHdTmax
    integer :: num_comp
    class(scalar_func), allocatable :: h_sol, h_liq
    class(scalar_func), allocatable :: h_sol_deriv, h_liq_deriv
    type(invf) :: my_invf
  contains
    procedure :: init
    procedure :: solve
    procedure :: solve_for_H_g
    procedure :: compute_H_res, compute_H_jac
    procedure :: compute_g_res, compute_g_jac
    procedure :: compute_C_liq, compute_C_sol
    procedure, private :: temp
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

    call params%get('dhdtmax', this%dHdTmax, stat, errmsg, default = 100.0_r8)
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
      call this%my_invf%compute(H, T1, T2, T, stat, errmsg) !TODO: [T1,T2] bracket
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
      dT = geut * (this%H_liq%eval([this%Teut])-this%H_sol%eval([this%Teut])) / this%dHdTmax
      Tmin = this%Teut - dT
      Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
    end if

    Hmin = this%H_sol%eval([Tmin])
    if (H <= Hmin) then ! pure solid
      !solve H = this%H_sol([T]) for T
      this%my_invf%func => H_sol_of_T
      call this%my_invf%compute(H, T1, T2, T, stat, errmsg) !TODO: [T1,T2] bracket
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
      T = this%Teut - dT*(Heut - H)/(Heut - Hmin)
      g = (H - this%H_sol%eval([T])) / (this%H_liq%eval([T]) - this%H_sol%eval([T]))
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

  end subroutine solve

  !! Given the temperature T and solute concentrations C, this solves for
  !! the consistent enthalpy H and liquid volume fraction according to the
  !! lever rule.

  subroutine solve_for_H_g(this, C, T, H, g)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: C(:,:), T(:)
    real(r8), intent(out) :: H(:), g(:)

    integer  :: j, stat
    real(r8) :: Tmin, Tmax, dT, geut, Hmin, Heut
    character(:), allocatable :: errmsg

    ASSERT(size(H) == size(T))
    ASSERT(size(g) == size(T))

    do j = 1, size(T)
      Tmax = this%temp(1.0_r8, C(:,j))
      if (T(j) >= Tmax) then ! pure liquid
        H(j) = this%H_liq%eval([T(j)])
        g(j) = 1
        cycle
      end if

      Tmin = this%temp(0.0_r8, C(:,j))
      if (Tmin < this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, C) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        dT = geut * (this%H_liq%eval([this%Teut])-this%H_sol%eval([this%Teut])) / this%dHdTmax
        Tmin = this%Teut - dT
        Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
      end if

      if (T(j) <= Tmin) then ! pure solid
        H(j) = this%H_sol%eval([T(j)])
        g(j) = 0
        cycle
      end if

      if (Tmin < this%Teut) then
        Hmin = this%H_sol%eval([Tmin])
        Tmin = this%Teut
      end if

      if (T(j) > Tmin) then ! 2-phase region
        ! solve T = this%temp(g, C(:,j)) for g
        this%my_invf%func => T_of_g
        call this%my_invf%compute(T(j), 0.0_r8, 1.0_r8, g(j), stat, errmsg)
        INSIST(stat == 0)
        H(j) = (1-g(j))*this%H_sol%eval([T(j)]) + g(j)*this%H_liq%eval([T(j)])
      else ! smeared eutectic transformation
        H(j) = Heut + (T(j) - this%Teut) * ((Heut - Hmin)/dT)
        g(j) = (H(j) - this%H_sol%eval([T(j)])) / (this%H_liq%eval([T(j)]) - this%H_sol%eval([T(j)]))
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C(:,j))
    end function

  end subroutine solve_for_H_g

  !! Compute the residual of the algebraic enthalpy relation given the liquid
  !! volume fraction g, enthalpy H, and temperature T.

  subroutine compute_H_res(this, g, H, T, r)
    class(multicomp_lever), intent(in) :: this
    real(r8), intent(in)  :: g(:), H(:), T(:)
    real(r8), intent(out) :: r(:)
    integer :: j
    do j = 1, size(g)
      r(j) = (1-g(j))*this%H_sol%eval([T(j)]) + g(j)*this%H_liq%eval([T(j)]) - H(j)
    end do
  end subroutine

  !! Compute the derivatives of the algebraic enthalpy relation with respect
  !! to the liquid volume fraction g, enthalpy H, and temperature T.

  subroutine compute_H_jac(this, g, H, T, drdg, drdH, drdT)
    class(multicomp_lever), intent(in) :: this
    real(r8), intent(in)  :: g(:), H(:), T(:)
    real(r8), intent(out) :: drdg(:), drdH(:), drdT(:)
    integer :: j
    do j = 1, size(g)
      drdg(j) = this%H_liq%eval([T(j)]) - this%H_sol%eval([T(j)])
      drdH(j) = -1
      drdT(j) = (1-g(j))*this%H_sol_deriv%eval([T(j)]) + g(j)*this%H_liq_deriv%eval([T(j)])
    end do
  end subroutine

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
        dT = geut * (this%H_liq%eval([this%Teut])-this%H_sol%eval([this%Teut])) / this%dHdTmax
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
      else ! smeared eutectic transformation
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
        dT = geut * (this%H_liq%eval([this%Teut])-this%H_sol%eval([this%Teut])) / this%dHdTmax
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
      else ! smeared eutectic transformation
        T = this%Teut - dT*(Heut - H(j))/(Heut - Hmin)
        g = (H(j) - this%H_sol%eval([T])) / (this%H_liq%eval([T]) - this%H_sol%eval([T]))
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

  end subroutine compute_C_sol

  !! Compute the residual of the algebraic liquid volume fraction relation
  !! given the liquid volume fraction g and enthalpy H.

  subroutine compute_g_res(this, C, g, H, r)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: C(:,:), g(:), H(:)
    real(r8), intent(out) :: r(:)

    integer :: j, stat
    character(:), allocatable :: errmsg
    real(r8) :: Tmax, Hmax, Tmin, Hmin, Heut, geut, T, dT

    ASSERT(size(H) == size(r))
    ASSERT(size(g) == size(r))

    do j = 1, size(r)
      Tmax = this%temp(1.0_r8, C(:,j))
      Hmax = this%H_liq%eval([Tmax])
      if (H(j) >= Hmax) then ! pure liquid
        r(j) = g(j) - 1
        cycle
      end if

      Tmin = this%temp(0.0_r8, C(:,j))
      if (Tmin < this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, C) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        dT = geut * (this%H_liq%eval([this%Teut])-this%H_sol%eval([this%Teut])) / this%dHdTmax
        Tmin = this%Teut - dT
        Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
      end if

      Hmin = this%H_sol%eval([Tmin])
      if (H(j) <= Hmin) then ! pure solid
        r(j) = g(j)
        cycle
      end if

      if (H(j) > merge(Heut, Hmin, Tmin < this%Teut)) then ! 2-phase region
        T = this%temp(g(j),C(:,j))
        r(j) = (1-g(j))*this%H_sol%eval([T]) + g(j)*this%H_liq%eval([T]) - H(j)
      else ! smeared eutectic transformation
        T = this%Teut - dT*(Heut - H(j))/(Heut - Hmin)
        r(j) = (1-g(j))*this%H_sol%eval([T]) + g(j)*this%H_liq%eval([T]) - H(j)
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C(:,j))
    end function

  end subroutine compute_g_res

  !! Compute the derivatives of the algebraic liquid volume fraction relation
  !! with respect to the liquid volume fraction g and enthalpy H.

  subroutine compute_g_jac(this, C, g, H, drdg, drdH)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: C(:,:), g(:), H(:)
    real(r8), intent(out) :: drdg(:), drdH(:)

    integer :: j, stat
    character(:), allocatable :: errmsg
    real(r8) :: Tmax, Hmax, Tmin, Hmin, Heut, geut, dT, T

    ASSERT(size(H) == size(drdH))
    ASSERT(size(g) == size(drdH))

    do j = 1, size(drdH)

      Tmax = this%temp(1.0_r8, C(:,j))
      Hmax = this%H_liq%eval([Tmax])
      if (H(j) >= Hmax) then ! pure liquid
        drdg(j) = 1
        drdH(j) = 0
        cycle
      end if

      Tmin = this%temp(0.0_r8, C(:,j))
      if (Tmin < this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, C) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        dT = geut * (this%H_liq%eval([this%Teut])-this%H_sol%eval([this%Teut])) / this%dHdTmax
        Tmin = this%Teut - dT
        Heut = (1-geut)*this%H_sol%eval([this%Teut]) + geut*this%H_liq%eval([this%Teut])
      end if

      Hmin = this%H_sol%eval([Tmin])
      if (H(j) <= Hmin) then ! pure solid
        drdg(j) = 1
        drdH(j) = 0
        cycle
      end if

      if (H(j) > merge(Heut, Hmin, Tmin < this%Teut)) then ! 2-phase region
        T = this%temp(g(j),C(:,j))
        drdg(j) = ((1-g(j))*this%H_sol_deriv%eval([T]) + g(j)*this%H_liq_deriv%eval([T])) &
                * temp_deriv(this, g(j), C(:,j)) + (this%H_liq%eval([T]) - this%H_sol%eval([T]))
        drdH(j) = -1
      else ! smeared eutectic transformation
        T = this%Teut - dT*(Heut - H(j))/(Heut - Hmin)
        drdg(j) = this%H_liq%eval([T]) - this%H_sol%eval([T])
        drdH(j) = ((1-g(j))*this%H_sol_deriv%eval([T]) + g(j)*this%H_liq_deriv%eval([T])) * &
                  (dT/(Heut-Hmin)) - 1
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C(:,j))
    end function

  end subroutine

end module
