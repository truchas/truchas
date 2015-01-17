project(PGSLibMpi)

message(STATUS "Configuring directory par-src")


# --- Include paths

# Add the include directory at this level
include_directories(${PGSLibMpi_SOURCE_DIR}/include)

# Add the MPI include directories
include_directories(${MPI_C_INCLUDE_PATH})

# --- Define source files
set(PGSLibMpi_SOURCE_FILES) 

# gath-scatt fiiles
set(GS_SOURCE_FILES  
              gath-scatt/gather-double.c
	      gath-scatt/scatter-float.c
              gath-scatt/gather-float.c
	      gath-scatt/gs-setup-c.c
	      gath-scatt/scatter-int.c
              gath-scatt/gather-int.c
	      gath-scatt/gs-util-c.c
	      gath-scatt/scatter-log.c
              gath-scatt/gather-log.c
	      gath-scatt/scatter-double.c)
set_source_files_properties(${GS_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS "-I${PGSLibMpi_SOURCE_DIR}/gath-scatt")
list(APPEND PGSLibMpi_SOURCE_FILES ${GS_SOURCE_FILES})			  

# indexing files
set(IDX_SOURCE_FILES    
	      indexing/hash.c
	      indexing/index-c.c) 
list(APPEND PGSLibMpi_SOURCE_FILES ${IDX_SOURCE_FILES})			  

# io files
list(APPEND PGSLibMpi_SOURCE_FILES io/io-c.c)

# reductions files
list(APPEND PGSLibMpi_SOURCE_FILES reductions/redux-c.c)

# scans files
set(SCANS_SOURCE_FILES
	      scans/scan-c-double.c
	      scans/scan-c-float.c
	      scans/scan-c-int.c)
set_source_files_properties(${SCANS_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS "-I${PGSLibMpi_SOURCE_DIR}/scans")
list(APPEND PGSLibMpi_SOURCE_FILES ${SCANS_SOURCE_FILES})			  

# send-recv files
set(SDRV_SOURCE_FILES
              send-rcv/attempt-receive-double.c
	      send-rcv/constrained-send-rcv-log.c
              send-rcv/attempt-receive-float.c
              send-rcv/attempt-receive-int.c
	      send-rcv/random-send-rcv-double.c
              send-rcv/attempt-receive-log.c
	      send-rcv/random-send-rcv-float.c
              send-rcv/random-send-rcv-int.c
              send-rcv/constrained-send-rcv-double.c
	      send-rcv/random-send-rcv-log.c
              send-rcv/constrained-send-rcv-float.c
              send-rcv/constrained-send-rcv-int.c)
set_source_files_properties(${SDRV_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS "-I${PGSLibMpi_SOURCE_DIR}/send-recv")
list(APPEND PGSLibMpi_SOURCE_FILES ${SDRV_SOURCE_FILES})			  

# shift files
set(SHIFT_SOURCE_FILES
              shift/shift-c-double.c
	      shift/shift-c-int.c
              shift/shift-c-float.c
	      shift/shift-c-log.c)
set_source_files_properties(${SHIFT_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS "-I${PGSLibMpi_SOURCE_DIR}/shift")
list(APPEND PGSLibMpi_SOURCE_FILES ${SHIFT_SOURCE_FILES})			  


# utility files
list(APPEND PGSLibMpi_SOURCE_FILES utility/tags-c.c utility/utility-c.c)


# --- Define the parallel library target

# Target name 
add_library(pgslib_mpi ${PGSLibMpi_SOURCE_FILES})
             
# Set output name to the original name scheme
# and define this as a Fortran library
set_target_properties(pgslib_mpi PROPERTIES
                      OUTPUT_NAME       pgslib-par
		      OUTPUT_NAME_DEBUG pgslib-par_g
		      LINK_LANGUAGE     Fortran)

# Define target link libraries		 
target_link_libraries(pgslib_mpi ${MPI_C_LIBRARIES})

# Define the installation location
install(TARGETS pgslib_mpi
        EXPORT  pgslib_mpi
        LIBRARY DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR})

