target_sources(truchas PRIVATE
  em_properties.F90
  em_heat_driver.F90
  induction_heat_solver_type.F90
  ih_hfield_factory_type.F90
  ih_hfield_func_type.F90
  solenoid_fields.F90
  elliptic_integrals.F90
  microwave_heat_solver_type.F90
  port_feed_func_factory.F90
  tdme_joule_heat_sim_type.F90
  tdme_model_type.F90
  tdme_solver_type.F90
  tdme_cg_solver_type.F90
  fdme_solver_type.F90
  fdme_model_type.F90
  fdme_minres_solver_type.F90
  fdme_mixed_minres_solver_type.F90
  fdme_mixed_zvector_type.F90
  fdme_vtk_graphics_proc.F90
  em_bc_factory_type.F90
  pec_bndry_func_type.F90
  nxh_bndry_func_type.F90
  fd_robin_bndry_func_type.F90
  wg_port_bc_plist_factory_type.F90
  mimetic_discretization.F90
  index_corrector_type.F90
  electromagnetics_namelist.F90
  electromagnetic_bc_namelist.F90
  induction_source_field_namelist.F90
  em_mesh_namelist.F90
)

if(USE_MUMPS)
  target_sources(truchas PRIVATE fdme_mumps_solver_type.F90)
endif()
