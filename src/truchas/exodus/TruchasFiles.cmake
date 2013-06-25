# Truchas files in directory
#   exodus

# List of files to  process
set(EXODUS_FILES)

# List of files to add to the Truchas library
set(EXODUS_SOURCE_FILES)

set(EXODUS_FILES
           exodus/exodus.F90
           exodus/exodus_errors.F90
           exodus/exodus_mesh_reader.F90
           exodus/exodus_mesh_type.F90
           exodus/exodus_mesh_utilities.F90
           exodus/exodus_mesh_writer.F90
           exodus/exodus_truchas_hack.F90)

set(EXODUS_FPP_FLAGS 
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(EXODUS_SOURCE_FILES
                         FILES ${EXODUS_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${EXODUS_FPP_FLAGS}
			 PROCESS_TARGET ProcessTruchasExodusFiles)
set_source_files_properties(${EXODUS_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS -I${NetCDF_INCLUDE_DIR})

# Update the Truchas library file list and targets		       
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${EXODUS_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ProcessTruchasExodusFiles)



