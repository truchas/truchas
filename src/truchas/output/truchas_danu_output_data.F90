module truchas_danu_output_data
  use,intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
#ifdef PATHSCALE_COMPILER_WORKAROUND
  type(c_ptr), save, public :: fid ! h5 file id
  type(c_ptr), save, public :: sid ! Danu simulation id
#else
  type(c_ptr), save, public :: fid = C_NULL_PTR ! h5 file id
  type(c_ptr), save, public :: sid = C_NULL_PTR ! Danu simulation id
#endif
end module truchas_danu_output_data
