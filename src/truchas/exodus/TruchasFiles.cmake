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

# Source file properties		       
set(EXOMESH_SOURCE_FILES_PROPS COMPILE_FLAGS -I${EXODUS_INCLUDE_DIR})
if ( TARGET ${EXODUS_BUILD_TARGET} )
  list(APPEND EXOMESH_SOURCE_FILE_PROPS 
              OBJECT_DEPENDS ${EXODUS_BUILD_TARGET})
endif()	    
set_source_files_properties(${EXOMESH_SOURCE_FILES} 
                            PROPERTIES ${EXOMESH_SOURCE_FILES_PROPS})
                            

# Update the Truchas library file list and targets		       
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${EXOMESH_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ProcessTruchasExoMeshFiles)



