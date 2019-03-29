# Truchas files in directory
#   bcs

# List of files to  process
set(BCS_FILES)

# List of files to add to the Truchas library
set(BCS_SOURCE_FILES)

set(BCS_FILES
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

set(BCS_FPP_FLAGS 
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(BCS_SOURCE_FILES
                         FILES ${BCS_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${OUTPUT_FPP_FLAGS}
			 PROCESS_TARGET ProcessTruchasBcsFiles)
set_source_files_properties(${BCS_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS "-I${PGSLib_MODULE_DIR}")


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${BCS_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS ProcessTruchasBcsFiles)
