! -*- Mode: f90 -*-

module test_data_module
  use UbikBase_module

  ! Constants denoting problem type.
  integer(Ubik_int_type), parameter :: Read_System=0
  integer(Ubik_int_type), parameter :: Saad_Full=1
  integer(Ubik_int_type), parameter :: Saad_ELL=2

  ! prob is compared with Read_System, Saad_Full and Saad_ELL
  integer(Ubik_int_type) :: prob=0

  ! storage_format is compared with Ubik_Full and Ubik_ELL
  integer(Ubik_int_type) :: storage_format

end module test_data_module
