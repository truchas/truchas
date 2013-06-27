# Truchas files in directory
#   nka

# List of files to  process
set(NKA_FILES)

# List of files to add to the Truchas library
set(NKA_SOURCE_FILES)

# Process target name
set(NKA_TARGET_NAME ProcessTruchasNkaFiles)


set(NKA_FILES nka/nka_type.F90)

set(NKA_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(NKA_SOURCE_FILES
                         FILES ${NKA_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${NKA_FPP_FLAGS}
			 PROCESS_TARGET ${NKA_TARGET_NAME})
set_source_files_properties(${NKA_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS -I${PGSLib_MODULE_DIR})


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${NKA_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${NKA_TARGET_NAME})


