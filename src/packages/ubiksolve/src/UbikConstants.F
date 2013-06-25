! -*- Mode: f90 -*-

module ubik_constants_module
  use ubik_types_module, only: int_type

  private :: int_type

  ! Generally useful constants.
  integer(int_type), parameter :: no = 0
  integer(int_type), parameter :: yes = 1

  integer(int_type), parameter :: false = 0
  integer(int_type), parameter :: true = 1

  ! Constants denoting storage types.
  integer(int_type), parameter :: Full = 1
  integer(int_type), parameter :: ELL = 8

  ! Constants denoting methods.
  integer(int_type), parameter :: method_min = -2
  integer(int_type), parameter :: precond_none = -2
  integer(int_type), parameter :: method_LU = -1
  integer(int_type), parameter :: method_CG = 0
  integer(int_type), parameter :: method_GMRES_Left_Radicati = 1
  integer(int_type), parameter :: method_BCGS = 2
  integer(int_type), parameter :: method_TFQMR = 3
  integer(int_type), parameter :: method_BCGSTAB = 4
  integer(int_type), parameter :: method_SOR = 5
  integer(int_type), parameter :: method_SSOR = 6
  integer(int_type), parameter :: method_Jacobi = 7
  integer(int_type), parameter :: method_IC = 8
  integer(int_type), parameter :: method_ILU = 9
  integer(int_type), parameter :: method_FGMRES = 10
  integer(int_type), parameter :: method_GMRES_Left_Saad = 11
  integer(int_type), parameter :: method_ApproxInv = 12
  integer(int_type), parameter :: method_GMRES_Right_Saad = 13
  integer(int_type), parameter :: method_LAMG = 14
  integer(int_type), parameter :: method_max = 14
  integer(int_type), parameter :: method_GMRES = method_GMRES_Left_Radicati

  ! Constants denoting return status values.
  integer(int_type), parameter :: success__ = 0
  integer(int_type), parameter :: converged__ = 0
  integer(int_type), parameter :: invalid_input__ = -1
  integer(int_type), parameter :: alloc_failure__ = -2
  integer(int_type), parameter :: internal_error__ = -3
  integer(int_type), parameter :: itmax_exceeded__ = -4
  integer(int_type), parameter :: breakdown__ = -5
  integer(int_type), parameter :: matmul_error__ = -6
  integer(int_type), parameter :: precond_error__ = -7
  integer(int_type), parameter :: next_iter__ = 1

end module ubik_constants_module
