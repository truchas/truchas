# Truchas files in directory
#   input

# List of files to  process
set(INPUT_FILES)

# List of files to add to the Truchas library
set(INPUT_SOURCE_FILES)

# Process target name
set(INPUT_TARGET_NAME ProcessTruchasInputFiles)


set(INPUT_FILES
         input/EM_input.F90
         input/bc_input_module.F90
         input/body_input_module.F90
         input/input_driver.F90
         input/interfaces_input_module.F90
         input/lin_solver_input.F90
         input/material_input_module.F90
         input/nonlin_solver_input.F90
         input/numerics_input_module.F90
         input/outputs_input_module.F90
         input/physics_input_module.F90
         input/probe_data_module.F90
         input/probe_input_module.F90
         input/region_data.F90
         input/region_input_module.F90
         input/solid_mechanics_namelist.F90
         input/legacy_flow_namelist.F90
)

set(INPUT_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(INPUT_SOURCE_FILES
                         FILES ${INPUT_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${INPUT_FPP_FLAGS}
			 PROCESS_TARGET ${INPUT_TARGET_NAME})
set(fc_flags -I${Danu_Fortran_MODULE_DIR})
list(APPEND fc_flags -I${UbikSolve_MODULE_DIR})
list(APPEND fc_flags -I${PGSLib_MODULE_DIR})
build_whitespace_string(INPUT_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${INPUT_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${INPUT_COMPILE_FLAGS})


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${INPUT_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${INPUT_TARGET_NAME})


