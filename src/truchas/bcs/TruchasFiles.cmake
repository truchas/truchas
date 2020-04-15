# Truchas files in directory
#   bcs

set(BCS_SOURCE_FILES
        bcs/Update_BCS.F90
        bcs/bc_atlases.F90
        bcs/bc_atlases_data_types.F90
        bcs/bc_atlases_internal.F90
        bcs/bc_atlases_util.F90
        bcs/bc_charts.F90
        bcs/bc_charts_atlases.F90
        bcs/bc_data_types.F90
        bcs/bc_displacement_init.F90
        bcs/bc_enum_types.F90
        bcs/bc_initialize.F90
        bcs/bc_operations.F90
        bcs/bc_operators.F90
        bcs/bc_pressure_init.F90
        bcs/bc_regions.F90
        bcs/bc_specifications.F90)

set_source_files_properties(${BCS_SOURCE_FILES} PROPERTIES
  COMPILE_FLAGS "-I${PGSLib_MODULE_DIR} -I${Truchas_utilities_dir}")


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${BCS_SOURCE_FILES})
