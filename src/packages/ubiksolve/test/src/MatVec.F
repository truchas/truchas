! -*- Mode: f90 -*-

module MatVec_module
  use UbikMatrix_module
  implicit none

  private
  public :: MatVec

contains

  subroutine MatVec_Full (coeff, x_vec, y, status)
    real(Ubik_real_type), intent(in), dimension(:,:) :: coeff

#include "MatVec-guts.F"

    return
  end subroutine MatVec_Full

  subroutine MatVec_ELL (coeff, x_vec, y, status)
    type(Ubik_ELL_matrix), intent(in) :: coeff

#include "MatVec-guts.F"

    return
  end subroutine MatVec_ELL

  subroutine MatVec (x_vec, y, status)
    use Full_data_module
    use ELL_data_module
    use test_data_module

    type(Ubik_vector_type), intent(inout) :: x_vec
    real(Ubik_real_type), intent(inout), dimension(:), target :: y
    integer(Ubik_int_type), intent(out) :: status

    select case (storage_format)
    case (Ubik_Full)
       call MatVec_Full (coeff_full, x_vec, y, status)
    case (Ubik_ELL)
       call MatVec_ELL (coeff_ell, x_vec, y, status)
    end select

    return
  end subroutine MatVec

end module MatVec_module
