# Truchas files in directory
#   exodus

# List of files to  process
set(EXOMESH_FILES)

# List of files to add to the Truchas library
set(EXOMESH_SOURCE_FILES)

set(EXOMESH_FILES
           exodus/exodus_c_binding.F90
           exodus/exodus_file_type.F90
           exodus/exodus_mesh_type.F90
           exodus/exodus_mesh_io.F90
           exodus/exodus_truchas_hack.F90)

set(EXOMESH_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(EXOMESH_SOURCE_FILES
                         FILES ${EXOMESH_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${EXOMESH_FPP_FLAGS}
			 PROCESS_TARGET ProcessTruchasExoMeshFiles)

# Update the Truchas library file list and targets		       
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${EXOMESH_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ProcessTruchasExoMeshFiles)



