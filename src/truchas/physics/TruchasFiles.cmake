# Truchas files in directory
#   physics

# - enclosure_radiation
set(PHYSICS_SOURCE_FILES
           physics/enclosure_radiation/netcdf_c_binding.F90
           physics/enclosure_radiation/netcdf_file_type.F90
           physics/enclosure_radiation/rad_encl_type.F90
           physics/enclosure_radiation/rad_encl_func_type.F90
           physics/enclosure_radiation/rad_solver_type.F90
           physics/enclosure_radiation/rad_problem_type.F90
           physics/enclosure_radiation/rad_encl_gmv.F90
           physics/enclosure_radiation/rad_problem_gmv.F90
           physics/enclosure_radiation/rad_encl_file_type.F90
           physics/enclosure_radiation/enclosure_radiation_namelist.F90
           physics/enclosure_radiation/rad_system_type.F90
           physics/enclosure_radiation/encl_vf_class.F90
           physics/enclosure_radiation/static_vf_class.F90
           physics/enclosure_radiation/facet_vf_type.F90
           physics/enclosure_radiation/patch_vf_type.F90
           physics/enclosure_radiation/moving_vf_type.F90
           physics/enclosure_radiation/vf_data_type.F90)


# - volume_tracking
list(APPEND PHYSICS_SOURCE_FILES
          physics/volume_tracking/advection_velocity_namelist.F90
	  physics/volume_tracking/brent_root_class.F90
	  physics/volume_tracking/locate_plane_nd_module.F90
	  physics/volume_tracking/multimat_cell_type.F90
	  physics/volume_tracking/near_zero_function.F90
	  physics/volume_tracking/plane_type.F90
	  physics/volume_tracking/polygon_type.F90
	  physics/volume_tracking/polyhedron_type.F90
	  physics/volume_tracking/pure_polyhedron_type.F90
	  physics/volume_tracking/volume_tracker_class.F90
          physics/volume_tracking/geometric_volume_tracker_type.F90
          physics/volume_tracking/simple_volume_tracker_type.F90
	  physics/volume_tracking/vtrack_driver.F90
	  physics/volume_tracking/cell_geometry_type.F90
	  physics/volume_tracking/locate_plane_os_function.F90
	  physics/volume_tracking/truncation_volume_type.F90
	  physics/volume_tracking/wisp_redistribution_type.F90)

# - heat_species_transport
list(APPEND PHYSICS_SOURCE_FILES
           physics/heat_species_transport/thermal_bc_namelist.F90
           physics/heat_species_transport/thermal_bc_factory_class.F90
           physics/heat_species_transport/thermal_bc_factory1_type.F90
           physics/heat_species_transport/species_bc_namelist.F90
           physics/heat_species_transport/species_bc_factory_class.F90
           physics/heat_species_transport/species_bc_factory1_type.F90
           physics/heat_species_transport/FHT_model_factory.F90
           physics/heat_species_transport/FHT_model_type.F90
           physics/heat_species_transport/FHT_norm_type.F90
           physics/heat_species_transport/FHT_precon_type.F90
           physics/heat_species_transport/FHT_solver_factory.F90
           physics/heat_species_transport/FHT_solver_type.F90
           physics/heat_species_transport/HTSD_idaesol_model_type.F90
           physics/heat_species_transport/HTSD_model_factory.F90
           physics/heat_species_transport/HTSD_model_type.F90
           physics/heat_species_transport/HTSD_norm_type.F90
           physics/heat_species_transport/HTSD_precon_type.F90
           physics/heat_species_transport/HTSD_solver_factory.F90
           physics/heat_species_transport/HTSD_solver_type.F90
           physics/heat_species_transport/HTSD_init_cond_type.F90
           physics/heat_species_transport/TofH_type.F90
           physics/heat_species_transport/data_layout_type.F90
           physics/heat_species_transport/diff_precon_type.F90
           physics/heat_species_transport/diffusion_matrix.F90
           physics/heat_species_transport/diffusion_solver.F90
           physics/heat_species_transport/diffusion_solver_data.F90
           physics/heat_species_transport/ds_source_input.F90
           physics/heat_species_transport/mesh_interop.F90
           physics/heat_species_transport/mfd_disc_type.F90
           physics/heat_species_transport/prop_mesh_func_type.F90
           physics/heat_species_transport/source_mesh_function.F90
           physics/heat_species_transport/upper_packed_matrix.F90
           physics/heat_species_transport/enthalpy_advector_class.F90
           physics/heat_species_transport/enthalpy_advector1_type.F90
           physics/heat_species_transport/evap_heat_flux_type.F90
           physics/heat_species_transport/evaporation_namelist.F90
           physics/heat_species_transport/thermal_source_namelist.F90
           physics/heat_species_transport/thermal_source_factory_type.F90
           physics/heat_species_transport/enthalpy_redistributor_type.F90
)

# - induction_heating
list(APPEND PHYSICS_SOURCE_FILES
           physics/induction_heating/altmesh_namelist.F90
           physics/induction_heating/CGSolver.F90
           physics/induction_heating/EM.F90
           physics/induction_heating/EM_boundary_data.F90
           physics/induction_heating/EM_data_proxy.F90
           physics/induction_heating/EM_graphics_output.F90
           physics/induction_heating/EM_properties.F90
           physics/induction_heating/GeometricModeler.F90
           physics/induction_heating/MaxwellBoundaryData.F90
           physics/induction_heating/MaxwellEddy.F90
           physics/induction_heating/data_explorer.F90
           physics/induction_heating/debug_EM.F90
           physics/induction_heating/elliptic_integrals.F90
           physics/induction_heating/field_probes.F90
           physics/induction_heating/mimetic_discretization.F90
           physics/induction_heating/solenoid_fields.F90
           physics/induction_heating/sparse_matrix.F90)

# - solid_mechanics
list(APPEND PHYSICS_SOURCE_FILES
           physics/solid_mechanics/smold_nlsol_model_type.F90
           physics/solid_mechanics/mech_bc_data_module.F90
           physics/solid_mechanics/node_op_setup_module.F90
           physics/solid_mechanics/node_operator_module.F90
           physics/solid_mechanics/solid_mech_constraints.F90
           physics/solid_mechanics/solid_mechanics_data.F90
           physics/solid_mechanics/solid_mechanics_input.F90
           physics/solid_mechanics/solid_mechanics_mesh.F90
           physics/solid_mechanics/solid_mechanics_module.F90
           physics/solid_mechanics/solid_mechanics_output.F90
           physics/solid_mechanics/VP_model_class.F90
           physics/solid_mechanics/MTS_VP_model_type.F90
           physics/solid_mechanics/power_law_VP_model_type.F90
           physics/solid_mechanics/viscoplastic_model_namelist.F90
           physics/solid_mechanics/viscoplasticity.F90
           physics/solid_mechanics/tm_density.F90)

# - additive manufacturing
list(APPEND PHYSICS_SOURCE_FILES
           physics/additive_manufacturing/ded_head_driver.F90
           physics/additive_manufacturing/ded_head_namelist.F90
           physics/additive_manufacturing/ded_head_type.F90
           physics/additive_manufacturing/laser_irrad_class.F90
           physics/additive_manufacturing/laser_irrad_factory.F90
           physics/additive_manufacturing/beam_laser_irrad_type.F90
           physics/additive_manufacturing/gauss_laser_irrad_type.F90)


# Set compile flags
include(BuildWhitespaceString)
set(fc_flags
  -I${PGSLib_MODULE_DIR}
  -I${Truchas_utilities_dir}
  -I${Truchas_solver_dir})

build_whitespace_string(PHYSICS_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${PHYSICS_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${PHYSICS_COMPILE_FLAGS})

set_source_files_properties(${TruchasExe_BINARY_DIR}/ER_file.f90 PROPERTIES
  COMPILE_FLAGS ${PHYSICS_COMPILE_FLAGS})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${PHYSICS_SOURCE_FILES})
