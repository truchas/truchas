target_sources(truchas PRIVATE
  netcdf_c_binding.F90
  netcdf_file_type.F90
  rad_encl_type.F90
  rad_encl_func_type.F90
  rad_solver_type.F90
  rad_problem_type.F90
  rad_encl_gmv.F90
  rad_problem_gmv.F90
  rad_encl_file_type.F90
  enclosure_radiation_namelist.F90
  rad_system_type.F90
  encl_vf_class.F90
  static_vf_class.F90
  facet_vf_type.F90
  patch_vf_type.F90
  moving_vf_type.F90
  vf_data_type.F90
)

add_subdirectory(Test)
