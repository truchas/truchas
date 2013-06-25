! -*- Mode: f90 -*-

!*******************************************************************
! $Id: Driver.F,v 1.14 2002/08/30 20:01:41 turner Exp $
!
!*******************************************************************
module Driver_module
#if !defined(absoft) && !defined(sgi) && !defined(sun)
  use Ubik_module
#endif
  use Iterate_module
  use MatVec_module
  use Precond_module
  use SetUpSaad_module
  use test_data_module
#ifdef use_PGSLib
  use PGSLib_module
#endif
#ifdef NAG_COMPILER
  use f90_unix_proc
#endif
  implicit none

  interface Driver
     module procedure Driver_Full
     module procedure Driver_ELL
  end interface

  private
  public :: Driver

contains

  !*******************************************************************

  subroutine Driver_Full (coeff_save, &
       nsizeb, nrows, dvalue, scale, plot, &
       bsave, ref, x, b, control, proc_time, elapsed_time)

#if defined(absoft) || defined(sgi) || defined(sun)
    use Ubik_module
#endif
    use Full_data_module, coeff => coeff_full, &
                          precond => precond_full

    real(Ubik_real_type), intent(inout), dimension(:,:) :: coeff_save

#include "Driver-guts.F"

  end subroutine Driver_Full

  !*******************************************************************

  subroutine Driver_ELL (coeff_save, &
       nsizeb, nrows, dvalue, scale, plot, &
       bsave, ref, x, b, control, proc_time, elapsed_time)

#if defined(absoft) || defined(sgi) || defined(sun)
    use Ubik_module
#endif
    use ELL_data_module, coeff => coeff_ell, &
                         precond => precond_ell

    type(Ubik_ELL_matrix), intent(inout) :: coeff_save

#include "Driver-guts.F"

  end subroutine Driver_ELL

  !*******************************************************************

end module Driver_module
