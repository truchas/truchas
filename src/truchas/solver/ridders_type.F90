#include "f90_assert.fpp"

module ridders_type

  use kinds, only: r8
  implicit none
  private

  type, public :: ridders
    real(r8) :: eps     ! convergence tolerance
    real(r8) :: maxitr  ! maximum number of iterations allowed
    real(r8) :: error = 0.0_r8  ! estimate of the error in the root
    integer  :: numitr = 0      ! number of iterations taken
  end type ridders
  
  public :: ridders_find_root
  
contains

  subroutine ridders_find_root (this, f, xmin, xmax, root, stat)
  
    type(ridders), intent(inout) :: this
    real(r8), intent(in) :: xmin, xmax
    real(r8), intent(out) :: root
    integer, intent(out) :: stat
    
    interface
      function f(x) result(fx)
        use kinds, only: r8
        real(r8), intent(in) :: x
        real(r8) :: fx
      end function
    end interface
    
    real(r8) :: a, b, c, m, fa, fb, fc, fm
    
    ASSERT(xmin <= xmax)
    a = xmin
    b = xmax
    
    this%error = 0.0_r8
    this%numitr = 0
    stat = 0
    
    fa = f(a)
    if (fa == 0.0_r8) then
      root = a
      return
    end if
    
    fb = f(b)
    if (fb == 0.0_r8) then
      root = b
      return
    end if
    
    if (fa == sign(fa,fb)) then ! interval doesn't bracket a root
      stat = -1
      return
    end if
    
    do  ! until converged
      if (this%numitr == this%maxitr) then
        stat = 1
        return
      end if
      this%numitr = this%numitr + 1
      !! Next approximate root C.
      m = 0.5d0*(a + b)
      fm = f(m)
      c = m + (m-a) * fm * sign(1.0d0/sqrt(fm*fm - fa*fb), fa)
      !! First convergence check.
      if (this%numitr > 1) then
        this%error = abs(c - root)
        if (this%error <= this%eps) then
          root = c
          return
        end if
      end if
      root = c
      !! Update the interval bracketing the root.
      fc = f(c)
      if (fc == 0.0_r8) then
        this%error = 0.0_r8
        return
      else if (fm /= sign(fm,fa)) then  ! [a,m]
        if (fm /= sign(fm,fc)) then ! [c,m]
          a = c; fa = fc
          b = m; fb = fm
        else ! [a,c]
          b = c; fb = fc
        end if
      else  ! [m,b]
        if (fm /= sign(fm,fc)) then ! [m,c]
          a = m; fa = fm
          b = c; fb = fc
        else  ! [c,b]
          a = c; fa = fc
        end if
      endif
      ASSERT(a < b)
      ASSERT(fa*fb < 0.0_r8)
      !! Second convergence check: [a,b] brackets the root.
      this%error = b - a
      if (this%error <= this%eps) return
    end do
    
  end subroutine ridders_find_root

end module ridders_type
