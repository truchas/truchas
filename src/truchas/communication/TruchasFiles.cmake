# Truchas files in directory
#   communication

# List of files to  process
set(COMM_FILES)

# List of files to add to the Truchas library
set(COMM_SOURCE_FILES)

set(COMM_FILES
         communication/ee_gather_module.F90
         communication/en_gather_module.F90
         communication/en_scatter_module.F90
         communication/gather_module.F90
         communication/gs_info_module.F90
         communication/gs_module.F90
         communication/gs_util.F90
         communication/gs_module.F90
         communication/nn_gather_module.F90
         communication/parallel_info_module.F90
         communication/parallel_util_module.F90
         communication/scatter_module.F90)
set(COMM_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities 
	-I${TruchasExe_SOURCE_DIR}/communication
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(COMM_SOURCE_FILES
                         FILES ${COMM_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${COMM_FPP_FLAGS}
			 PROCESS_TARGET ProcessTruchasCommFiles)
set_source_files_properties(${COMM_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS "-I${PGSLib_MODULE_DIR} -I${Danu_Fortran_MODULE_DIR}")

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${COMM_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ProcessTruchasCommFiles)

