! -*- Mode: f90 -*-

module ELL_data_module
  use UbikMatrix_module

  implicit none
  type(Ubik_ELL_matrix), save :: coeff_ell, precond_ell

  private
  public :: Ubik_ELL_matrix, &
            coeff_ell, &
            precond_ell

end module ELL_data_module
