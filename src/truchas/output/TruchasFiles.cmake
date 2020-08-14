# Truchas files in directory
#   discrete_operators

set(OUTPUT_SOURCE_FILES
          output/cycle_output_module.F90
          output/diagnostics_module.F90
          output/edit_module.F90
          output/gap_output.F90
          output/output_control.F90
          output/output_utilities.F90
          output/probes_type.F90
          output/probes_driver.F90
          output/probe_namelist.F90
          output/truchas_probe_field_factory_type.F90
          output/truchas_danu_output.F90
          output/truchas_danu_output_data.F90
          output/truchas_danu_output_tools.F90)


set(fc_flags -I${Truchas_utilities_dir})
build_whitespace_string(OUTPUT_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${OUTPUT_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${OUTPUT_COMPILE_FLAGS})


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${OUTPUT_SOURCE_FILES})
