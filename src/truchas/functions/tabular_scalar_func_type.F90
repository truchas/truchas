!!
!! TABULAR_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class SCALAR_FUNC.  This
!! implementation defines a tabular function given by user-specified data
!! points with intervening interpolation.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! By default, linear interpolation is used between data points. If the function
!! is instantiated using the SMOOTH=.TRUE. option, then function evaluation will
!! use smooth C1 interpolation using Hermite cubics on each data interval.  The
!! slopes at the points are given by a modified form of the Akima algorithm [1],
!! used by Matlabs tablelookup function.  The interpolation tries to preserve
!! the slope and avoid undulations where the data suggests a flat region.
!!
!! [1] Akima, Hiroshi. "A new method of interpolation and smooth curve fitting
!! based on local procedures." Journal of the ACM (JACM) , vol 17, no 4, 1970.
!!
!! IMPLEMENTATION NOTES
!!
!! No tabular_ad_scalar_func constructor is defined for tabular_scalar_func
!! objects created with the smooth option (there is an assertion check). It
!! is feasible, though complex, to define one.  However, the anti-derivative
!! of a non-smooth tabular function is already C1, so it does not seem worth
!! the effort at this time.
!!

#include "f90_assert.fpp"

module tabular_scalar_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  implicit none
  private

  type, extends(scalar_func), public :: tabular_scalar_func
    !private  ! scalar_func_tools needs access
    real(r8), allocatable :: x(:), y(:), t(:)
    integer :: dim = 1
    integer :: extrap
  contains
    procedure :: eval
  end type tabular_scalar_func

  !! Defined constructor
  interface tabular_scalar_func
    procedure tabular_scalar_func_value
  end interface

  type, extends(scalar_func), public :: tabular_ad_scalar_func
    !private  ! scalar_func_tools needs access
    real(r8), allocatable :: x(:), y(:), c(:)
    integer :: dim = 1
    integer :: extrap
  contains
    procedure :: eval => eval_ad
  end type tabular_ad_scalar_func

  !! Defined constructor
  interface tabular_ad_scalar_func
    procedure tabular_ad_scalar_func_deriv
  end interface

  integer, parameter :: EXTRAP_NEAREST = 1, EXTRAP_LINEAR = 2

contains

  !! Constructor for TABULAR_SCALAR_FUNC objects
  function tabular_scalar_func_value(x, y, dim, smooth, extrap) result(f)

    real(r8), intent(in) :: x(:), y(:)
    integer, intent(in), optional :: dim
    logical, intent(in), optional :: smooth
    character(*), intent(in), optional :: extrap
    type(tabular_scalar_func) :: f

    integer :: j, n
    real(r8) :: wl, wr, sw
    real(r8), allocatable :: m(:)

    ASSERT(size(x) > 1)
    ASSERT(size(y) == size(x))

    f%x = x
    f%y = y
    if (present(dim)) f%dim = dim
    if (present(smooth)) then
      if (smooth) allocate(f%t(size(x)))
    end if

    f%extrap = EXTRAP_NEAREST
    if (present(extrap)) then
      select case (extrap)
      case ('nearest')
        f%extrap = EXTRAP_NEAREST
      case ('linear')
        f%extrap = EXTRAP_LINEAR
      case default
        INSIST(.false.)
      end select
    end if

    !! Tangent slope at nodes for Akima interpolation algorithm.
    if (allocated(f%t)) then
      if (size(x) == 2) then
        f%t = (y(2)-y(1))/(x(2)-x(1))
      else
        !! Divided differences of table data.
        allocate(m(0:size(x)+3))
        do j = 2, size(x)
          m(j) = (y(j)-y(j-1))/(x(j)-x(j-1))
        end do
        !! Add 2 data points at either end using quadratic extrapolation.
        m(1) = 2*m(2) - m(3)
        m(0) = 2*m(1) - m(2)
        n = size(x)
        m(n+1) = 2*m(n) - m(n-1)
        m(n+2) = 2*m(n+1) - m(n)
        !! Interpolant tangent slope at nodes.
        do j = 1, size(x)
          wl = abs(m(j+2) - m(j+1))
          wr = abs(m(j) - m(j-1))
          sw = wl + wr
          if (sw == 0.0_r8) then
            !! Akima would set wl = wr = 0.5, producing undulations on both
            !! sides.  This choice employed by Matlab preserves flatness of
            !! the more horizontal side.
            if (abs(m(j)) < abs(m(j+1))) then
              wl = 1.0_r8
              wr = 0.0_r8
            else
              wl = 0.0_r8
              wr = 1.0_r8
            end if
            sw = 1.0_r8
          end if
          f%t(j) = (wl*m(j) + wr*m(j+1))/sw
        end do
      end if
    end if

  end function tabular_scalar_func_value

  function eval(this, x) result(fx)
    class(tabular_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)  ! only x(this%dim) is used
    real(r8) :: fx, xdim, t
    integer :: n, j, j1, j2
    n = size(this%x)
    xdim = x(this%dim)
    if (xdim <= this%x(1)) then
      select case (this%extrap)
      case (EXTRAP_NEAREST)
        fx = this%y(1)
      case (EXTRAP_LINEAR)
        if (allocated(this%t)) then
          fx = this%y(1) + this%t(1)*(xdim-this%x(1))
        else
          fx = this%y(1) + ((this%y(2)-this%y(1))/(this%x(2)-this%x(1)))*(xdim-this%x(1))
        end if
      end select
    else if (xdim >= this%x(n)) then
      select case (this%extrap)
      case (EXTRAP_NEAREST)
        fx = this%y(n)
      case (EXTRAP_LINEAR)
        if (allocated(this%t)) then
          fx = this%y(n) + this%t(n)*(xdim-this%x(n))
        else
          fx = this%y(n) + ((this%y(n)-this%y(n-1))/(this%x(n)-this%x(n-1)))*(xdim-this%x(n))
        end if
      end select
    else
      !! Binary search to find the interval x(j1) < x <= x(j2), j2 = j1+1.
      j1 = 1; j2 = n
      do while (j2 - j1 > 1)
        j = (j1 + j2) / 2
        if (xdim > this%x(j)) then
          j1 = j
        else
          j2 = j
        end if
      end do
      if (allocated(this%t)) then
        !! Hermite cubit interpolation over the interval.
        t = (this%x(j2)-xdim)/(this%x(j2)-this%x(j1))
        fx = this%y(j1)*(3-2*t)*t**2 + this%y(j2)*(1+2*t)*(1-t)**2 &
              + (this%x(j2)-this%x(j1))*t*(1-t)*(this%t(j1)*t-this%t(j2)*(1-t))
      else
        !! Linear interpolation over the interval.
        fx = ((this%x(j2)-xdim)*this%y(j1) + (xdim-this%x(j1))*this%y(j2))/(this%x(j2)-this%x(j1))
      end if
    end if
  end function eval

  !! Constructor for TABULAR_AD_SCALAR_FUNC objects
  function tabular_ad_scalar_func_deriv(f, x0, y0) result(adf)
    type(tabular_scalar_func), intent(in) :: f
    real(r8), intent(in) :: x0, y0
    type(tabular_ad_scalar_func) :: adf
    integer :: j, n
    INSIST(.not.allocated(f%t))
    n = size(f%x)
    allocate(adf%x(n), adf%y(n), adf%c(0:n))
    adf%x = f%x
    adf%y(1) = 0.0_r8
    do j = 1, n-1
      adf%y(j+1) = adf%y(j) + 0.5_r8*(f%y(j+1)+f%y(j))*(f%x(j+1)-f%x(j))
      adf%c(j) = -0.5_r8*(f%y(j+1)-f%y(j))*(f%x(j+1)-f%x(j))
    end do
    adf%c(0) = f%y(1)
    adf%c(n) = f%y(n)
    adf%extrap = f%extrap
    adf%y = adf%y - adf%eval([x0]) + y0 ! shift to satisfy adf(x0)==y0
    adf%dim = f%dim
  end function tabular_ad_scalar_func_deriv

  function eval_ad(this, x) result(fx)
    class(tabular_ad_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)  ! only x(this%dim) is used
    real(r8) :: fx, xdim, a1, a2
    integer :: n, j, j1, j2
    n = size(this%x)
    xdim = x(this%dim)
    if (xdim <= this%x(1)) then
      select case (this%extrap)
      case (EXTRAP_NEAREST)
        fx = this%y(1) + this%c(0)*(xdim - this%x(1))
      case (EXTRAP_LINEAR)
        fx = this%y(1) + this%c(0)*(xdim - this%x(1)) &
           - this%c(1)*((xdim-this%x(1))/(this%x(2)-this%x(1)))**2
      end select
    else if (xdim >= this%x(n)) then
      select case (this%extrap)
      case (EXTRAP_NEAREST)
        fx = this%y(n) + this%c(n)*(xdim - this%x(n))
      case (EXTRAP_LINEAR)
        fx = this%y(n) + this%c(n)*(xdim - this%x(n)) &
           - this%c(n-1)*((xdim - this%x(n))/(this%x(n)-this%x(n-1)))**2
      end select
    else
      !! Binary search to find the interval x(j1) < x <= x(j2), j2 = j1+1.
      j1 = 1; j2 = n
      do while (j2 - j1 > 1)
        j = (j1 + j2) / 2
        if (xdim > this%x(j)) then
          j1 = j
        else
          j2 = j
        end if
      end do
      !! Quadratic interpolation within the interval
      a1 = (this%x(j2)-xdim)/(this%x(j2)-this%x(j1))
      a2 = (xdim-this%x(j1))/(this%x(j2)-this%x(j1))
      fx = a1*this%y(j1) + a2*this%y(j2) + a1*a2*this%c(j1)
    end if
  end function eval_ad

end module tabular_scalar_func_type
