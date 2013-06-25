project(ChacoSRC)

message(STATUS "Entering the src directory")

# --- Add include directories
include_directories(${ChacoSRC_SOURCE_DIR}/main)
include_directories(${ChacoSRC_SOURCE_DIR}/internal)
include_directories(${ChacoSRC_SOURCE_DIR}/symlq)
include_directories(${ChacoSRC_SOURCE_DIR}/tinvit)

# --- Source files in subdirectories 

# Define a list of subdirectories
set(ChacoSRC_SUBDIRS
                   assign
		   bpmatch
		   coarsen
		   connect
		   eigen
		   graph
		   inertial
		   input
		   internal
		   klspiff
		   klvspiff
		   main
		   matlab
		   misc
		   optimize
		   refine_map
		   refine_part
		   submain
		   symmlq
		   tinvit
		   util)

# Loop through each subdirectory and append to CHACO_SOURCE_FILES
set(CHACO_SOURCE_FILES)
foreach( subdir ${ChacoSRC_SUBDIRS})
  #message(STATUS "\tAdding files from: src/${subdir}")
  file(GLOB subdir_files "${subdir}/*.c")
  list(APPEND CHACO_SOURCE_FILES ${subdir_files})
endforeach(subdir)

# --- Define the 'chaco' library target and where to install

# - Create the library target
add_library(chaco ${CHACO_SOURCE_FILES})

install(TARGETS chaco
        EXPORT  chaco
        LIBRARY DESTINATION ${Chaco_LIBRARY_INSTALL_DIR}
        ARCHIVE DESTINATION ${Chaco_LIBRARY_INSTALL_DIR})
