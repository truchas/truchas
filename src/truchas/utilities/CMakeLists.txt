target_sources(truchas PRIVATE
  f90_assert.F90
  file_utility.F90
  graph_type.F90
  input_utilities.F90
  integer_set_type.F90
  integer_set_type_wavl.F90
  integer_map_type.F90
  math_constants.F90
  permutations.F90
  process_info_module.F90
  signal_handler.F90
  sort_module.F90
  sort_utilities.F90
  string_set_type.F90
  string_utilities.F90
  truchas_env.F90
  truchas_logging_services.F90
  truchas_timers.F90
  utilities_module.F90
  gmv/gmvwrite_c_binding.F90
  get_process_size.c
  make_directory.c
  gmv/gmvwrite.c
  serialization_tools.F90
)

add_subdirectory(Test)
