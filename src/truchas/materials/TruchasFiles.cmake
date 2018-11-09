# Truchas files in directory
#   material

# List of files to  process
set(MAT_FILES)

# List of files to add to the Truchas library
set(MAT_SOURCE_FILES)

# Process target name
set(MAT_TARGET_NAME ProcessTruchasMaterialFiles)


set(MAT_FILES
       materials/material_interop.F90
       materials/matl_mesh_func_type.F90
       materials/material_property.F90
       materials/material_system.F90
       materials/material_system_namelist.F90
       materials/material_table.F90
       materials/material_utilities.F90
       materials/phase_namelist.F90
       materials/phase_property_table.F90
       materials/legacy_matl_api.F90
       materials/legacy_matl_adapter_class.F90
       materials/legacy_matl_adapter0_type.F90
       materials/legacy_matl_adapter1_type.F90
       materials/legacy_matl_adapter2_type.F90
)

set(MAT_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
        ${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(MAT_SOURCE_FILES
                         FILES ${MAT_FILES}
                         FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
                         FPP_FLAGS ${MAT_FPP_FLAGS}
                         PROCESS_TARGET ${MAT_TARGET_NAME})

                       
# Define compile flags                       
include(BuildWhitespaceString)
build_whitespace_string(mat_source_flags_string 
                        -I${PGSLib_MODULE_DIR} -I${Danu_Fortran_MODULE_DIR} -I${PETACA_MODULE_DIR})
set_source_files_properties(${MAT_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${mat_source_flags_string})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${MAT_SOURCE_FILES})       
list(APPEND Truchas_PROCESS_TARGETS ${MAT_TARGET_NAME})


