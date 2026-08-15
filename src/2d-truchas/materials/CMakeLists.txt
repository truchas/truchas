target_sources(truchas-extra PRIVATE
  new_mesh_func_class.F90
  cell_matl_prop_func_type.F90
  material_composition_type.F90
)

add_subdirectory(test)
