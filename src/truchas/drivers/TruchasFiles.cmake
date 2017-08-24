# Truchas files in directory
#   drivers

# List of files to  process
set(DRIVERS_FILES)

# List of files to add to the Truchas library
set(DRIVERS_SOURCE_FILES)

# Process target name
set(DRIVERS_TARGET_NAME ProcessTruchasDriversFiles)

set(DRIVERS_FILES
           drivers/drivers.F90
           drivers/physics_module.F90
           drivers/time_step_module.F90
           drivers/time_step_sync_type.F90
           drivers/simulation_event_queue.F90
           drivers/hijack_truchas.F90)

set(DRIVERS_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(DRIVERS_SOURCE_FILES
                         FILES ${DRIVERS_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${DRIVERS_FPP_FLAGS}
			 PROCESS_TARGET ${DRIVERS_TARGET_NAME})

# Set compile flags		       
include(BuildWhitespaceString)
set(fc_flags -I${Danu_Fortran_MODULE_DIR})
list(APPEND fc_flags -I${PGSLib_MODULE_DIR})
build_whitespace_string(DRIVERS_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${DRIVERS_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${DRIVERS_COMPILE_FLAGS})

# drivers.F90 requires extra flags
list(APPEND fc_flags -I${Danu_Fortran_MODULE_DIR})
list(APPEND fc_flags -I${UbikSolve_MODULE_DIR})
build_whitespace_string(DRIVERS_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${TruchasExe_BINARY_DIR}/drivers.f90
                            COMPILE_FLAGS ${DRIVERS_COMPILE_FLAGS})


# Add the C source files
set(DRIVERS_C_SOURCE_FILES drivers/runinfo.c)
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DRIVERS_C_SOURCE_FILES})
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DRIVERS_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${DRIVERS_TARGET_NAME})



