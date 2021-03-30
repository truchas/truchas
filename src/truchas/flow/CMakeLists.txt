target_sources(truchas PRIVATE
  flow_domain_types.F90
  flow_input_utils.F90
  fischer_guess_type.F90
  flow_bc_factory_type.F90
  flow_bc_type.F90
  flow_driver.F90
  flow_operators.F90
  turbulence_model_class.F90
  default_turb_model_type.F90
  algebraic_turb_model_type.F90
  turbulence_models.F90
  flow_projection_type.F90
  flow_prediction_type.F90
  flow_props_type.F90
  flow_type.F90
  flow_namelist.F90
  flow_bc_namelist.F90
  flow_surface_tension_bc_type.F90
  flow_solver_namelists.F90
  turbulence_namelist.F90
)

add_subdirectory(Test)
