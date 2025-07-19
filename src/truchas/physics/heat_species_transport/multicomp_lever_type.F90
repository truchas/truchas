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
    real(r8), allocatable :: liq_slope(:), part_coef(:), C0(:)
    real(r8) :: T_f, Teut
    integer :: num_comp
    class(scalar_func), allocatable :: h_sol, h_liq
    class(scalar_func), allocatable :: h_sol_deriv, h_liq_deriv
    type(invf) :: my_invf
    ! temporary stuff for testing/development
    real(r8) :: solidus_temp, liquidus_temp
  contains
    procedure :: init
    procedure, private :: temp
    procedure :: solve
    procedure :: compute_H, compute_dHdT, compute_dHdg
    procedure :: compute_r, compute_drdH, compute_drdg
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

    !! Ridders solver to invert functions
    call this%my_invf%init(eps=1d-9,maxadj=10)

    !! Extract solidus and liquidus temperatures (TEMPORARY HACK)
    block
      use multiphase_matl_type
      select type (matl)
      type is (multiphase_matl)
        this%solidus_temp  = matl%pc_seq(1)%pc%solidus_temp()
        this%liquidus_temp = matl%pc_seq(1)%pc%liquidus_temp()
      end select
    end block

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

  subroutine solve(this, H, C, T1, T2, T, g, stat, errmsg)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: h, C(:), T1, T2
    real(r8), intent(out) :: g, T
    integer,  intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: Tmax, Hmax, Tmin, Hmin, hsol, hliq

    Tmax = this%temp(1.0_r8, C)
    Hmax = this%H_liq%eval([Tmax, C])
    if (H >= Hmax) then ! pure liquid
      g = 1
      !solve H = this%H_liq([T, C]) for T
      this%my_invf%func => H_liq_of_T
      call this%my_invf%compute(H, T1, T2, T, stat, errmsg) !TODO: [T1,T2] bracket
      if (stat /= 0) return
      return
    end if

    Tmin = max(this%temp(0.0_r8, C), this%Teut)
    Hmin = this%H_sol%eval([Tmin, C])
    if (H <= Hmin) then ! pure solid
      g = 0
      !solve H = this%H_sol([T,C]) for T
      this%my_invf%func => H_sol_of_T
      call this%my_invf%compute(H, T1, T2, T, stat, errmsg) !TODO: [T1,T2] bracket
      if (stat /= 0) return
      return
    end if

    if (Tmin == this%Teut) then ! eutectic transformation may be in play
      !solve this%Teut = this%temp(g, C) for g
      this%my_invf%func => T_of_g
      call this%my_invf%compute(Tmin, 0.0_r8, 1.0_r8, g, stat, errmsg)
      INSIST(stat == 0)
      Hmin = (1-g)*this%H_sol%eval([Tmin, C]) + g*this%H_liq%eval([Tmin, C])
    end if

    if (H >= Hmin) then ! 2-phase region
      !solve H = (1-g)*this%H_sol%eval([this%temp(g,C),C]) + g*this%H_liq%eval([this%temp(g,C),C]) for g
      this%my_invf%func => H_of_g
      call this%my_invf%compute(H, 0.0_r8, 1.0_r8, g, stat, errmsg)
      if (stat /= 0) return
      T = this%temp(g, C)
    else ! eutectic transformation at constant temperature
      g = (H - this%H_sol%eval([this%Teut,C])) / (this%H_liq%eval([this%Teut,C]) - this%H_sol%eval([this%Teut,C]))
      T = this%Teut
    end if

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, C)
    end function

    real(r8) function H_sol_of_T(x) result(H)
      real(r8), intent(in) :: x
      H = this%H_sol%eval([x, C])
    end function

    real(r8) function H_liq_of_T(x) result(H)
      real(r8), intent(in) :: x
      H = this%H_liq%eval([x, C])
    end function

    real(r8) function H_of_g(x) result(H)
      real(r8), intent(in) :: x
      real(r8) :: hsol, hliq
      H = (1-x)*this%H_sol%eval([this%temp(x,C),C]) + x*this%H_liq%eval([this%temp(x,C),C])
    end function

  end subroutine solve

  subroutine compute_H(this, T, g, H)
    class(multicomp_lever), intent(in) :: this
    real(r8), intent(in)  :: T(:), g(:)
    real(r8), intent(out) :: H(:)
    real(r8) :: arg(0:this%num_comp)
    integer :: j
    ASSERT(size(T) == size(H))
    ASSERT(size(g) == size(H))
    arg(1:) = this%C0
    do j = 1, size(H)
      arg(0) = T(j)
      H(j) = (1-g(j))*this%H_sol%eval(arg) + g(j)*this%H_liq%eval(arg)
    end do
  end subroutine

  subroutine compute_dHdT(this, T, g, dHdT)
    class(multicomp_lever), intent(in) :: this
    real(r8), intent(in)  :: T(:), g(:)
    real(r8), intent(out) :: dHdT(:)
    integer :: j
    real(r8) :: arg(0:this%num_comp)
    ASSERT(size(T) == size(dHdT))
    ASSERT(size(g) == size(dHdT))
    arg(1:) = this%C0
    do j = 1, size(dHdT)
      arg(0) = T(j)
      dHdT(j) = (1-g(j))*this%H_sol_deriv%eval(arg) + g(j)*this%H_liq_deriv%eval(arg)
    end do
  end subroutine

  subroutine compute_dHdg(this, T, g, dHdg)
    class(multicomp_lever), intent(in) :: this
    real(r8), intent(in)  :: T(:), g(:)
    real(r8), intent(out) :: dHdg(:)
    integer :: j
    real(r8) :: arg(0:this%num_comp)
    ASSERT(size(T) == size(dHdg))
    ASSERT(size(g) == size(dHdg))
    arg(1:) = this%C0
    do j = 1, size(dHdg)
      arg(0) = T(j)
      dHdg(j) = this%H_liq%eval(arg) - this%H_sol%eval(arg)
    end do
  end subroutine

  subroutine compute_r(this, H, g, r)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: H(:), g(:)
    real(r8), intent(out) :: r(:)

    integer :: j, stat
    character(:), allocatable :: errmsg
    real(r8) :: Tmax, Hmax, Tmin, Hmin, geut, arg(0:this%num_comp)

    ASSERT(size(H) == size(r))
    ASSERT(size(g) == size(r))

    arg(1:) = this%C0

    do j = 1, size(r)
      Tmax = this%temp(1.0_r8, this%C0)
      Hmax = this%H_liq%eval([Tmax, this%C0])
      if (H(j) >= Hmax) then ! pure liquid
        r(j) = g(j) - 1
        cycle
      end if

      Tmin = max(this%temp(0.0_r8, this%C0), this%Teut)
      Hmin = this%H_sol%eval([Tmin, this%C0])
      if (H(j) <= Hmin) then ! pure solid
        r(j) = g(j)
        cycle
      end if

      if (Tmin == this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, this%C0) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        arg(0) = this%Teut
        Hmin = (1-geut)*this%H_sol%eval(arg) + geut*this%H_liq%eval(arg)
      end if

      if (H(j) >= Hmin) then ! 2-phase region
        arg(0) = this%temp(g(j),this%C0)
        r(j) = (1-g(j))*this%H_sol%eval(arg) + g(j)*this%H_liq%eval(arg) - H(j)
      else ! eutectic transformation at constant temperature
        arg(0) = this%Teut
        r(j) = g(j) - (H(j) - this%H_sol%eval(arg)) / (this%H_liq%eval(arg) - this%H_sol%eval(arg))
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, this%C0)
    end function

  end subroutine

  subroutine compute_drdH(this, H, g, drdH)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: H(:), g(:)
    real(r8), intent(out) :: drdH(:)

    integer :: j, stat
    character(:), allocatable :: errmsg
    real(r8) :: Tmax, Hmax, Tmin, Hmin, geut, arg(1+this%num_comp)

    ASSERT(size(H) == size(drdH))
    ASSERT(size(g) == size(drdH))

    do j = 1, size(drdH)

      drdH(j) = 0 ! default

      Tmax = this%temp(1.0_r8, this%C0)
      Hmax = this%H_liq%eval([Tmax, this%C0])
      if (H(j) >= Hmax) cycle ! pure liquid

      Tmin = max(this%temp(0.0_r8, this%C0), this%Teut)
      Hmin = this%H_sol%eval([Tmin, this%C0])
      if (H(j) <= Hmin) cycle ! pure solid

      if (Tmin == this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, this%C0) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        arg = [this%Teut, this%C0]
        Hmin = (1-geut)*this%H_sol%eval(arg) + geut*this%H_liq%eval(arg)
      end if

      if (H(j) >= Hmin) then ! 2-phase region
        drdH(j) = -1
      else ! eutectic transformation at constant temperature
        arg = [this%Teut, this%C0]
        drdH(j) = - 1.0_r8 / (this%H_liq%eval(arg) - this%H_sol%eval(arg))
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, this%C0)
    end function

  end subroutine

  subroutine compute_drdg(this, H, g, drdg)

    class(multicomp_lever), intent(inout) :: this
    real(r8), intent(in)  :: H(:), g(:)
    real(r8), intent(out) :: drdg(:)

    integer :: j, stat
    character(:), allocatable :: errmsg
    real(r8) :: Tmax, Hmax, Tmin, Hmin, geut, arg(1+this%num_comp)

    ASSERT(size(H) == size(drdg))
    ASSERT(size(g) == size(drdg))

    do j = 1, size(drdg)
      Tmax = this%temp(1.0_r8, this%C0)
      Hmax = this%H_liq%eval([Tmax, this%C0])
      if (H(j) >= Hmax) then ! pure liquid
        drdg(j) = 1
        cycle
      end if

      Tmin = max(this%temp(0.0_r8, this%C0), this%Teut)
      Hmin = this%H_sol%eval([Tmin, this%C0])
      if (H(j) <= Hmin) then ! pure solid
        drdg(j) = 1
        cycle
      end if

      if (Tmin == this%Teut) then ! eutectic transformation may be in play
        !solve this%Teut = this%temp(geut, this%C0) for geut
        this%my_invf%func => T_of_g
        call this%my_invf%compute(this%Teut, 0.0_r8, 1.0_r8, geut, stat, errmsg)
        INSIST(stat == 0)
        arg = [this%Teut, this%C0]
        Hmin = (1-geut)*this%H_sol%eval(arg) + geut*this%H_liq%eval(arg)
      end if

      if (H(j) >= Hmin) then ! 2-phase region
        arg = [this%temp(g(j),this%C0), this%C0]
        drdg(j) = ((1-g(j))*this%H_sol_deriv%eval(arg) + g(j)*this%H_liq_deriv%eval(arg)) &
                * temp_deriv(this, g(j), this%C0) + (this%H_liq%eval(arg) - this%H_sol%eval(arg))
      else ! eutectic transformation at constant temperature
        drdg(j) = 1
      end if
    end do

  contains

    real(r8) function T_of_g(x) result(T)
      real(r8), intent(in) :: x
      T = this%temp(x, this%C0)
    end function

  end subroutine

end module
