!!
!! MPOLY_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class SCALAR_FUNC.
!! This implementation defines a general multivariable polynomial with
!! user-specified coefficients and integral exponents.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!

#include "f90_assert.fpp"

module mpoly_scalar_func_type

  use kinds, only: r8
  use scalar_func_class
  implicit none
  private

  type, extends(scalar_func), public :: mpoly_scalar_func
    !private  ! scalar_func_tools needs access
    real(r8), allocatable :: x0(:)      ! reference point
    integer,  allocatable :: expon(:,:) ! array of exponents
    real(r8), allocatable :: coef(:)    ! array of coefficients
  contains
    procedure :: eval
  end type mpoly_scalar_func

  !! Defined constructor
  interface mpoly_scalar_func
    procedure mpoly_scalar_func_value
  end interface

contains

  !! Constructor for MPOLY_SCALAR_FUNC objects.
  function mpoly_scalar_func_value (c, e, x0) result (f)

    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:,:)
    real(r8), intent(in), optional :: x0(:)
    type(mpoly_scalar_func) :: f

    INSIST(size(c) > 0 .and. size(c) == size(e,dim=2))
    INSIST(size(e,dim=1) > 1)

    allocate(f%x0(size(e,dim=1)))
    if (present(x0)) then
      INSIST(size(x0) == size(f%x0))
      f%x0 = x0
    else
      f%x0 = 0.0_r8
    end if

    f%expon = e
    f%coef = c

  end function mpoly_scalar_func_value

  function eval (this, x) result(fx)

    class(mpoly_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx

    integer  :: i, j
    real(r8) :: t

    ASSERT(size(x) >= size(this%x0))

    fx = 0.0_r8
    do j = 1, size(this%coef)
      t = this%coef(j)
      do i = 1, size(this%x0)
        if (this%expon(i,j) /= 0) then
          t = t * (x(i)-this%x0(i))**this%expon(i,j)
        end if
      end do
      fx = fx + t
    end do

  end function eval

end module mpoly_scalar_func_type
