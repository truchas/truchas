target_sources(truchas PRIVATE
  ih_driver.F90
  ih_legacy_bc.F90
  geometry_model_type.F90
  em_properties.F90
  EM_graphics_output.F90
  data_explorer.F90
  field_probes.F90
  mimetic_discretization.F90
  em_mesh_namelist.F90
  electromagnetics_namelist.F90
  electromagnetic_bc_namelist.F90
  induction_source_field_namelist.F90
  ih_source_factory_type.F90
  ih_source_func_type.F90
  solenoid_fields.F90
  elliptic_integrals.F90
  em_bc_factory_type.F90
  nxh_bndry_func_type.F90
  pec_bndry_func_type.F90
  tdme_model_type.F90
  tdme_solver_type.F90
  tdme_cg_solver_type.F90
  tdme_joule_heat_sim_type.F90
)
