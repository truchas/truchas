# Truchas files in directory
#   functions

# List of files to  process
set(FUNC_FILES)

# List of files to add to the Truchas library
set(FUNC_SOURCE_FILES)

# Process target name
set(FUNC_TARGET_NAME ProcessTruchasFunctionFiles)


set(FUNC_FILES
        functions/call_dll_scafun.F90
        functions/function_namelist.F90
        functions/function_table.F90
        functions/scalar_functions.F90
        functions/user_scafun_stub.F90)

set(FUNC_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(FUNC_SOURCE_FILES
                         FILES ${FUNC_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${FUNC_FPP_FLAGS}
			 PROCESS_TARGET ${FUNC_TARGET_NAME})
set_source_files_properties(${FUNC_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS -I${PGSLib_MODULE_DIR})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${FUNC_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${FUNC_TARGET_NAME})


