# Truchas files in directory
#   material

set(MAT_SOURCE_FILES
       materials/material_class.F90
       materials/material_database_type.F90
       materials/material_factory.F90
       materials/material_namelist.F90
       materials/single_phase_matl_type.F90
       materials/multi_phase_matl_type.F90
       materials/phase_change_class.F90
       materials/smooth_phase_change_type.F90
       materials/tabular_phase_change_type.F90
       materials/phase_change_factory.F90
       materials/material_model_type.F90
       materials/matl_mesh_func_type.F90
       materials/matl_prop_class.F90
       materials/avg_matl_prop_type.F90
       materials/avg_phase_prop_type.F90
       materials/equil_temp_type.F90
       materials/material_model_driver.F90
       materials/material_utilities.F90
)

# Define compile flags
include(BuildWhitespaceString)
build_whitespace_string(mat_source_flags_string
  -I${PGSLib_MODULE_DIR}
  -I${PETACA_MODULE_DIR}
  -I${Truchas_utilities_dir})
set_source_files_properties(${MAT_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${mat_source_flags_string})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${MAT_SOURCE_FILES})
