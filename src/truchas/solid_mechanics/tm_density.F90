!!
!! TM_DENSITY
!!
!! This module provides a procedure for defining a SCALAR_FUNC object that
!! describes the material density as a function of temperature for use by
!! the thermo-mechanics physics kernel.  The function has the following form
!!
!!    \[ \rho(T) = \rho_0 \exp(-3 \int_{T_0}^{T} \alpha(T) dT) \]
!!
!! where $\rho_0$ and $\T_0$ are the reference density and temperature, and
!! $\alpha(T)$ is the linear coefficient of thermal expansion.
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
!! This module defines the extension TM_DENSITY_FUNC of the abstract base class
!! SCALAR_FUNC.  However, the type is private; only the following procedure
!! for instantiating a polymorphic SCALAR_FUNC variable of that type is needed;
!! application code should only need to use the type-bound class method (EVAL).
!!
!!  CALL ALLOC_TM_DENSITY_FUNC (RHO, REF_DENS, REF_TEMP, CTE, STAT, ERRMSG)
!!    allocates and defines the allocatable SCALAR_FUNC class argument RHO.
!!    REF_DENS and REF_TEMP are the reference density and temperature, and
!!    CTE is a SCALAR_FUNC class object that gives the linear coefficient of
!!    thermal expansion. CTE must be of a form that ALLOC_SCALAR_FUNC_ANTIDERIV
!!    can be used with, namely a constant function or a single-variable
!!    polynomial.  STAT returns a non-zero value if the antiderivative of CTE
!!    cannot be created, and an explanatory error message is returned in the
!!    deferred-length allocatable character argument ERRMSG.
!!

#include "f90_assert.fpp"

module tm_density

  use kinds, only: r8
  use scalar_func_class
  implicit none
  private

  public :: alloc_tm_density_func

  type, extends(scalar_func), private :: tm_density_func
    private
    real(r8) :: d0 = 0.0_r8 ! reference density
    class(scalar_func), allocatable :: int_cte  ! indefinite integral of the linear CTE
  contains
    procedure :: eval
  end type tm_density_func

contains

  function eval (this, x) result (fx)
    class(tm_density_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    fx =  this%d0 * exp(-3*this%int_cte%eval(x))
  end function eval

  subroutine alloc_tm_density_func (rho, ref_dens, ref_temp, cte, stat, errmsg)

    use scalar_func_tools, only: alloc_scalar_func_antideriv

    class(scalar_func), allocatable, intent(out) :: rho
    real(r8), intent(in) :: ref_dens  ! reference density
    real(r8), intent(in) :: ref_temp  ! reference temperature
    class(scalar_func),intent(in) :: cte ! linear coefficient of thermal expansion
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(tm_density_func), allocatable :: this

    allocate(this)
    this%d0 = ref_dens
    call alloc_scalar_func_antideriv (cte, ref_temp, 0.0_r8, this%int_cte, stat, errmsg)
    call move_alloc (this, rho) ! move it to the polymorphic return variable

  end subroutine alloc_tm_density_func

end module tm_density
