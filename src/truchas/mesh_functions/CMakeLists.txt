target_sources(truchas PRIVATE
  bndry_vfunc_class.F90
  bndry_face_vfunc_type.F90
  bndry_func1_class.F90
  bndry_func2_class.F90
  bndry_func3_class.F90
  intfc_func2_class.F90
  intfc_func3_class.F90
  bndry_face_func_type.F90
  bndry_face_group_builder_type.F90
  intfc_link_group_builder_type.F90
  cell_group_builder_type.F90
  vflux_bndry_func_type.F90
  htc_bndry_func_type.F90
  mtc_bndry_func_type.F90
  rad_bndry_func_type.F90
  htc_intfc_func_type.F90
  mtc_intfc_func_type.F90
  rad_intfc_func_type.F90
  scalar_mesh_func_class.F90
  scalar_mesh_func2_class.F90
  scalar_cell_func_type.F90
  scalar_cell_func1_type.F90
  scalar_cell_func2_type.F90
  scalar_mesh_multifunc_type.F90
  bndry_edge_group_builder_type.F90
  bndry_cfunc1_class.F90
  bndry_face_cfunc_type.F90
)

#TODO?
#if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#  list(APPEND mesh_func_source_flags "-standard-semantics -assume nostd_mod_proc_name")
#endif()
