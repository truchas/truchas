target_sources(truchas PRIVATE
  func_table.F90
  function_namelist.F90
  vfunction_namelist.F90
  complex_function_namelist.F90
  complex_vfunction_namelist.F90
  scalar_func_class.F90
  const_scalar_func_type.F90
  poly_scalar_func_type.F90
  mpoly_scalar_func_type.F90
  tabular_scalar_func_type.F90
  fptr_scalar_func_type.F90
  scalar_func_factories.F90
  scalar_func_containers.F90
  vector_func_containers.F90
  scalar_func_tools.F90
  scalar_func_map_type.F90
  vector_func_map_type.F90
  vector_func_class.F90
  const_vector_func_type.F90
  tabular_vector_func_type.F90
  fptr_vector_func_type.F90
  vector_func_factories.F90
  inverse_func_class.F90
  complex_scalar_func_class.F90
  complex_scalar_func_containers.F90
  complex_scalar_func_factories.F90
  const_complex_scalar_func_type.F90
  complex_vector_func_class.F90
  complex_vector_func_containers.F90
  complex_vector_func_factories.F90
  const_complex_vector_func_type.F90
  fptr_complex_vector_func_type.F90
  const_complex_scalar_func_type.F90
  dl_complex_scalar_func_type.F90
  dl_complex_vector_func_type.F90
)

if(ENABLE_DYNAMIC_LOADING)
  target_sources(truchas PRIVATE
    dl_scalar_func_type.F90
    dl_vector_func_type.F90
  )
endif()

add_subdirectory(Test)
