!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module root_finder_type

  use kinds, only: r8

  type, public :: root_finder
    integer :: maxitr
    integer :: numitr
    integer :: max_try
    real(r8) :: eps
    real(r8) :: error
    real(r8) :: root
    logical :: is_bracketed
    logical :: found_root
  contains
    procedure :: init => init_root_finder
    procedure :: bracket_root
    procedure :: compute_root => compute_root_via_ridders
  end type root_finder
  
  contains

    subroutine init_root_finder(this, maxitr, numitr, max_try, eps)

      class(root_finder), intent(inout) :: this
      integer, intent(in) :: maxitr
      integer, intent(in) :: numitr
      integer, intent(in) :: max_try
      real(r8), intent(in) :: eps

      this%maxitr = maxitr
      this%numitr = numitr
      this%max_try = max_try
      this%eps = eps
      this%is_bracketed = .false.
      this%found_root = .false.
      this%root = -1


    end subroutine init_root_finder


    subroutine bracket_root(this, f, params, a, b)

      class(root_finder), intent(inout) :: this
      real(r8), intent(in) :: params(:)
      real(r8), intent(out) :: a, b

      interface
        function f(x, params) 
          use kinds, only: r8
          real(r8), intent(in) :: x
          real(r8), intent(in) :: params(:)
          real(r8) :: f
        end function
      end interface

      integer :: j
      real(r8) :: d, fa, fb
 
      a = 0.0d0
      d = 1.0d0

      ! First check if the root is in [0,1]
      b = d
      fa = f(a, params)
      fb = f(b, params)

      if ( fa < 0 .and. fb > 0 ) then
        this%is_bracketed = .true.
        return 
      else
        b = 0
      end if

      ! Not in [0,1]---starting checking wider and wider intervals
      j = 0
      do while (fb < 0 .and. j < this%max_try )
        a = b
        b = b + d
        fa = f(a, params)
        fb = f(b, params)
        d = d*10
        j = j+1
      end do


      fa = f(a, params)
      fb = f(b, params)
      if ( fa < 0 .and. fb > 0 ) then
        this%is_bracketed = .true.
      else
        this%is_bracketed = .false.
      end if


    end subroutine bracket_root



    subroutine compute_root_via_ridders(this, f, params, xmin, xmax)

      class(root_finder), intent(inout) :: this
      real(r8), intent(in) :: xmin, xmax, params(:)
    
      interface
        function f(x, params) 
          use kinds, only: r8
          real(r8), intent(in) :: x, params(:)
          real(r8) :: f
        end function
      end interface
    
      integer :: numitr, maxitr
      real(r8) :: a, b, c, m, fa, fb, fc, fm, error, &
                  eps
    

      numitr = this%numitr
      maxitr = this%maxitr
      eps = this%eps

      a = xmin
      b = xmax
    
      error = 0.0d0
      numitr = 0
      this%found_root = .false.
    
      fa = f(a, params) 
      if (fa == 0.0d0) then
        this%found_root = .true.
        this%root = a
        this%error = 0
        return
      end if
    
      fb = f(b, params) 
      if (fb == 0.0d0) then
        this%found_root = .true.
        this%root = b
        this%error = 0
        return
      end if
    
      if (fa == sign(fa,fb)) then ! interval doesn't bracket a root
        return
      end if
    
      do  ! until converged

        if (numitr == maxitr) then
          this%found_root = .false.
          return
        end if
        numitr = numitr + 1
        !! Next approximate root C.
        m = 0.5d0*(a + b)
        fm = f(m, params) 
        c = m + (m-a) * fm * sign(1.0d0/sqrt(fm*fm - fa*fb), fa)
        !! First convergence check.
        if (numitr > 1) then
          this%error = abs(c - this%root)
          if (this%error <= eps) then
            this%root = c
            this%found_root = .true.
            return
          end if
        end if
        this%root = c
        !! Update the interval bracketing the root.
        fc = f(c, params)  
        if (fc == 0.0d0) then
          this%error = 0.0d0
          this%found_root = .true.
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

      !! Second convergence check: [a,b] brackets the root.
        this%error = b - a
        if (this%error <= eps) then
          this%found_root = .true.
          return
        end if
      end do
    
    end subroutine compute_root_via_ridders



end module root_finder_type


