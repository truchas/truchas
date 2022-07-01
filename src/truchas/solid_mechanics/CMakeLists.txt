target_sources(truchas PRIVATE
  solid_mechanics_driver.F90
  solid_mechanics_type.F90
  solid_mechanics_namelist.F90
  solid_mechanics_bc_namelist.F90
  viscoplastic_model_namelist.F90
  integration_geometry_type.F90
  integration_cell_type.F90
  bndry_ip_func_type.F90

  sm_bc_utilities.F90
  sm_normal_traction_bc_type.F90
  sm_bc_class.F90
  sm_bc_manager_type.F90
  sm_bc_c0d1_type.F90
  sm_bc_c0d2_type.F90
  sm_bc_c0d3_type.F90
  sm_bc_c1d0_type.F90
  sm_bc_c1d1_type.F90
  sm_bc_c1d2_type.F90
  sm_bc_c2d0_type.F90
  sm_bc_c2d1_type.F90
  sm_bc_c3d0_type.F90
  sm_bc_face_list_type.F90
  sm_bc_node_list_type.F90
  sm_bc_list_type.F90

  sm_material_model_type.F90
  viscoplastic_model_type.F90
  viscoplastic_solver_type.F90
  viscoplastic_material_model_types.F90
  viscoplastic_jacob_idaesol_model_type.F90
  viscoplastic_jfree_idaesol_model_type.F90
  sm_model_type.F90
  sm_ds_precon_type.F90
  sm_nlsol_model_type.F90
  tm_density.F90
)
