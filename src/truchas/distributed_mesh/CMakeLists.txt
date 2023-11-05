target_sources(truchas PRIVATE
  bitfield_type.F90
  base_mesh_class.F90
  unstr_base_mesh_class.F90
  facet_hash_type.F90
  facet_table_type.F90
  face_neighbor_table_type.F90
  mesh_manager.F90
  parallel_permutations.F90
  cell_topology.F90
  cell_topology_2d.F90
  cell_geometry.F90
  unstr_mesh_type.F90
  unstr_mesh_factory.F90
  unstr_mesh_tools.F90
  unstr_mesh_gmv.F90
  unstr_2d_mesh_type.F90
  unstr_2d_mesh_factory.F90
  unstr_2d_mesh_tools.F90
  ext_exodus_mesh_type.F90
  exodus_mesh_tools.F90
  simplex_topology.F90
  simplex_geometry.F90
  simpl_mesh_type.F90
  simpl_mesh_factory.F90
  simpl_mesh_tools.F90
  simpl_mesh_gmv.F90
  xdmf_file_type.F90
)

# not needed?
#target_include_directories(truchas PRIVATE ${CMAKE_SOURCE_DIR})
