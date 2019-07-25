# Truchas files in directory
#   discrete_operators

# List of files to  process
set(OUTPUT_FILES)

# List of files to add to the Truchas library
set(OUTPUT_SOURCE_FILES)

# Process target name
set(OUTPUT_TARGET_NAME ProcessTruchasOutputFiles)


set(OUTPUT_FILES
          output/cycle_output_module.F90
          output/diagnostics_module.F90
          output/edit_module.F90
          output/gap_output.F90
          output/interface_output_module.F90
          output/output_control.F90
          output/output_utilities.F90
          output/probes_type.F90
          output/probes_driver.F90
          output/probe_namelist.F90
          output/truchas_probe_field_factory_type.F90
          output/truchas_danu_output.F90
          output/truchas_danu_output_data.F90
          output/truchas_danu_output_tools.F90)

set(OUTPUT_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
        -I${TruchasExe_SOURCE_DIR}/output
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(OUTPUT_SOURCE_FILES
                         FILES ${OUTPUT_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${OUTPUT_FPP_FLAGS}
			 PROCESS_TARGET ${OUTPUT_TARGET_NAME})
set(fc_flags -I${UbikSolve_MODULE_DIR})
build_whitespace_string(OUTPUT_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${OUTPUT_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${OUTPUT_COMPILE_FLAGS})


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${OUTPUT_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS ${OUTPUT_TARGET_NAME})
