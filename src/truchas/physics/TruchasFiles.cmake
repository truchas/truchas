# Truchas files in directory
#   physics

# List of files to  process
set(PHYSICS_FILES)

# List of files to add to the Truchas library
set(PHYSICS_SOURCE_FILES)

# Process target name
set(PHYSICS_TARGET_NAME ProcessTruchasPhysicsFiles)

# - enclosure_radiation
list(APPEND PHYSICS_FILES
           physics/enclosure_radiation/ER_dist_encl.F90
           physics/enclosure_radiation/ER_driver.F90
           physics/enclosure_radiation/ER_driver_gmv.F90
           physics/enclosure_radiation/ER_encl_func.F90
           physics/enclosure_radiation/ER_file.F90
           physics/enclosure_radiation/ER_input.F90
           physics/enclosure_radiation/ER_solver.F90
           physics/enclosure_radiation/ER_solver_gmv.F90
           physics/enclosure_radiation/ER_system.F90)

# - fluid_flow
list(APPEND PHYSICS_FILES
           physics/fluid_flow/flow_phase_change.F90
           physics/fluid_flow/fluid_data_module.F90
           physics/fluid_flow/fluid_flow_module.F90
           physics/fluid_flow/fluid_type_module.F90
           physics/fluid_flow/fluid_utilities_module.F90)

# - fluid_flow/advection	 
list(APPEND PHYSICS_FILES
           physics/fluid_flow/advection/advect_volume_module.F90
           physics/fluid_flow/advection/advection_data.F90
           physics/fluid_flow/advection/advection_module.F90
           physics/fluid_flow/advection/hoadvection.F90
           physics/fluid_flow/advection/limiter.F90
           physics/fluid_flow/advection/limiter_data.F90
           physics/fluid_flow/advection/limiter_module.F90
           physics/fluid_flow/advection/riemann_module.F90)

# - fluid_flow/body_force	 
list(APPEND PHYSICS_FILES
           physics/fluid_flow/body_force/body_data_module.F90
           physics/fluid_flow/body_force/body_force_module.F90)

# - fluid_flow/porous_drag
list(APPEND PHYSICS_FILES
           physics/fluid_flow/porous_drag/porous_drag_data.F90
           physics/fluid_flow/porous_drag/porous_drag_module.F90)

# - fluid_flow/predictor
list(APPEND PHYSICS_FILES
           physics/fluid_flow/predictor/predictor_module.F90
           physics/fluid_flow/predictor/y_eq_Ax_vel.F90)

# - fluid_flow/projection	 
list(APPEND PHYSICS_FILES
           physics/fluid_flow/projection/coordinates_module.F90
           physics/fluid_flow/projection/fischer_module.F90
           physics/fluid_flow/projection/projection_data_module.F90
           physics/fluid_flow/projection/projection_module.F90
           physics/fluid_flow/projection/y_eq_Ax_prs.F90)

# - fluid_flow/surface_tension
list(APPEND PHYSICS_FILES
           physics/fluid_flow/surface_tension/kernel_interpolation_module.F90
           physics/fluid_flow/surface_tension/surface_tension_module.F90)

# - fluid_flow/viscous
list(APPEND PHYSICS_FILES
           physics/fluid_flow/viscous/turbulence_module.F90
           physics/fluid_flow/viscous/viscous_data_module.F90
           physics/fluid_flow/viscous/viscous_module.F90)

# - fluid_flow/vof
list(APPEND PHYSICS_FILES
           physics/fluid_flow/vof/flux_volume_module.F90
           physics/fluid_flow/vof/interface_module.F90
           physics/fluid_flow/vof/interface_triangle_module.F90
           physics/fluid_flow/vof/locate_plane_module.F90
           physics/fluid_flow/vof/mollify.F90
           physics/fluid_flow/vof/truncate_volume_module.F90
           physics/fluid_flow/vof/vof_data_module.F90
           physics/fluid_flow/vof/volume_track_module.F90)

# - heat_species_transport
list(APPEND PHYSICS_FILES
           physics/heat_species_transport/FHT_model_factory.F90
           physics/heat_species_transport/FHT_model_type.F90
           physics/heat_species_transport/FHT_norm_type.F90
           physics/heat_species_transport/FHT_precon_type.F90
           physics/heat_species_transport/FHT_solver_factory.F90
           physics/heat_species_transport/FHT_solver_type.F90
           physics/heat_species_transport/HTSD_BDF2_model.F90
           physics/heat_species_transport/HTSD_model_factory.F90
           physics/heat_species_transport/HTSD_model_type.F90
           physics/heat_species_transport/HTSD_norm_type.F90
           physics/heat_species_transport/HTSD_precon_type.F90
           physics/heat_species_transport/HTSD_solver_factory.F90
           physics/heat_species_transport/HTSD_solver_type.F90
           physics/heat_species_transport/TofH_callback.F90
           physics/heat_species_transport/TofH_type.F90
           physics/heat_species_transport/boomer_amg_precon_type.F90
           physics/heat_species_transport/boundary_data.F90
           physics/heat_species_transport/data_layout_type.F90
           physics/heat_species_transport/diff_precon_type.F90
           physics/heat_species_transport/diffusion_matrix.F90
           physics/heat_species_transport/diffusion_solver.F90
           physics/heat_species_transport/diffusion_solver_data.F90
           physics/heat_species_transport/ds_boundary_condition_input.F90
           physics/heat_species_transport/ds_interface_condition_input.F90
           physics/heat_species_transport/ds_source_input.F90
           physics/heat_species_transport/ds_utilities.F90
           physics/heat_species_transport/interface_data.F90
           physics/heat_species_transport/mesh_interop.F90
           physics/heat_species_transport/mfd_disc_type.F90
           physics/heat_species_transport/parallel_csr_matrix.F90
           physics/heat_species_transport/property_mesh_function.F90
           physics/heat_species_transport/source_mesh_function.F90
           physics/heat_species_transport/ssor_precon_type.F90
           physics/heat_species_transport/upper_packed_matrix.F90)

# - induction_heating
list(APPEND PHYSICS_FILES
           physics/induction_heating/CGSolver.F90
           physics/induction_heating/EM.F90
           physics/induction_heating/EM_boundary_data.F90
           physics/induction_heating/EM_data_proxy.F90
           physics/induction_heating/EM_graphics_output.F90
           physics/induction_heating/EM_hex_tet_mapping.F90
           physics/induction_heating/EM_utilities.F90
           physics/induction_heating/GeometricModeler.F90
           physics/induction_heating/MaxwellBoundaryData.F90
           physics/induction_heating/MaxwellEddy.F90
           physics/induction_heating/data_explorer.F90
           physics/induction_heating/debug_EM.F90
           physics/induction_heating/elliptic_integrals.F90
           physics/induction_heating/field_probes.F90
           physics/induction_heating/joule_xml_utilities.F90
           physics/induction_heating/mimetic_discretization.F90
           physics/induction_heating/solenoid_fields.F90
           physics/induction_heating/sparse_matrix.F90
           physics/induction_heating/system_io.F90)

# - properties	 
list(APPEND PHYSICS_FILES
           physics/properties/physical_constants.F90
           physics/properties/property_data_module.F90
           physics/properties/property_module.F90)

# - solid_mechanics 	 
list(APPEND PHYSICS_FILES
           physics/solid_mechanics/mech_bc_data_module.F90
           physics/solid_mechanics/node_op_setup_module.F90
           physics/solid_mechanics/node_operator_module.F90
           physics/solid_mechanics/solid_mech_constraints.F90
           physics/solid_mechanics/solid_mechanics_data.F90
           physics/solid_mechanics/solid_mechanics_module.F90
           physics/solid_mechanics/viscoplasticity.F90)

# Preprocess flags	 
set(PHYSICS_FPP_FLAGS
        -I${TruchasExe_SOURCE_DIR}/utilities
        -I${TruchasExe_SOURCE_DIR}/solver
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(PHYSICS_SOURCE_FILES
                         FILES ${PHYSICS_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${PHYSICS_FPP_FLAGS}
			 PROCESS_TARGET ProcessTruchasPhysicsFiles)

# Set compile flags
include(BuildWhitespaceString)
if (ENABLE_PGSLib)
  set(fc_flags -I${PGSLib_MODULE_DIR})
endif()
if (ENABLE_UbikSolve)
  list(APPEND fc_flags -I${UbikSolve_MODULE_DIR})
endif()
if(ENABLE_Danu)
  list(APPEND fc_flags -I${Danu_Fortran_MODULE_DIR})
endif()
build_whitespace_string(PHYSICS_COMPILE_FLAGS ${fc_flags} -I${NetCDF_INCLUDE_DIR})
set_source_files_properties(${PHYSICS_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${PHYSICS_COMPILE_FLAGS})

set_source_files_properties(${TruchasExe_BINARY_DIR}/ER_file.f90 PROPERTIES
                            COMPILE_FLAGS "${PHYSICS_COMPILE_FLAGS} -I${NetCDF_INCLUDE_DIR}")

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${PHYSICS_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ProcessTruchasPhysicsFiles)





            


