module truchas_danu_output_data
  use,intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  type(c_ptr), save, public :: fid = C_NULL_PTR ! h5 file id
  type(c_ptr), save, public :: sid = C_NULL_PTR ! Danu simulation id
end module truchas_danu_output_data
