# Truchas files in directory
#   fpa

# List of files to  process
set(FPA_FILES)

# List of files to add to the Truchas library
set(FPA_SOURCE_FILES)

# Process target name
set(FPA_TARGET_NAME ProcessTruchasFpaFiles)


set(FPA_FILES
       fpa/fixed_point_accelerator.F90
       fpa/fpa_kinds.F90)

set(FPA_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(FPA_SOURCE_FILES
                         FILES ${FPA_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${FPA_FPP_FLAGS}
			 PROCESS_TARGET ${FPA_TARGET_NAME})
set_source_files_properties(${FPA_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS -I${PGSLib_MODULE_DIR})


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${FPA_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${FPA_TARGET_NAME})


