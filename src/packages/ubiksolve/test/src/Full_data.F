! -*- Mode: f90 -*-

module Full_data_module
  use UbikMatrix_module, only: Ubik_real_type

  implicit none
  real(Ubik_real_type), allocatable, dimension(:,:), save :: coeff_full, precond_full

  private
  public :: coeff_full, precond_full

end module Full_data_module
