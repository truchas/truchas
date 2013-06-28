project(TruchasExe)

message(STATUS "configuring the Truchas binary and library directory")

# --- CMake modules
include(FortranPreprocessFiles)

# --- Include paths

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})


# --- Define the source files
set(Truchas_LIBRARY_SOURCE_FILES)
set(Truchas_PROCESS_TARGETS)

# - bcs
message(STATUS "Processing bcs")
include(bcs/TruchasFiles.cmake)

# - communication
message(STATUS "Processing communication")
include(communication/TruchasFiles.cmake)

# - discrete_operators
message(STATUS "Processing discrete_operators")
include(discrete_operators/TruchasFiles.cmake)

# - distributed_mesh
message(STATUS "Processing distributed_mesh")
include(distributed_mesh/TruchasFiles.cmake)

# - drivers
message(STATUS "Processing drivers")
include(drivers/TruchasFiles.cmake)

# - exodus
message(STATUS "Processing exodus")
include(exodus/TruchasFiles.cmake)

# - nka
message(STATUS "Processing nka")
include(nka/TruchasFiles.cmake)

# - functions
message(STATUS "Processing functions")
include(functions/TruchasFiles.cmake)

# - grid_mapping
message(STATUS "Processing grid_mapping")
include(grid_mapping/TruchasFiles.cmake)

# - input
message(STATUS "Processing input")
include(input/TruchasFiles.cmake)

# - materials
message(STATUS "Processing materials")
include(materials/TruchasFiles.cmake)

# - ode
message(STATUS "Processing ode")
include(ode/TruchasFiles.cmake)

# - output
message(STATUS "Processing output")
include(output/TruchasFiles.cmake)

# - physics
message(STATUS "Processing physics")
include(physics/TruchasFiles.cmake)

# - setup
message(STATUS "Processing setup")
include(setup/TruchasFiles.cmake)

# - solver
message(STATUS "Processing solver")
include(solver/TruchasFiles.cmake)

# - utilities
message(STATUS "Processing utilities")
include(utilities/TruchasFiles.cmake)

# --- Custom target to process Fortran files
add_custom_target(ProcessTruchasFiles
                  DEPENDS ${Truchas_PROCESS_TARGETS}
		  COMMENT "Processing all Truchas Fortran files")

# --- Library

# - Define the source files for the library
add_library(truchas ${Truchas_LIBRARY_SOURCE_FILES})

# - Define the dependencies for this library

# Procees targets
set(Truchas_DEPENDS ${Truchas_PROCESS_TARGETS})

# PGSLib
if (ENABLE_PGSLib)
  list(APPEND Truchas_DEPENDS pgslib_fe)
  if (ENABLE_MPI)
    list(APPEND Truchas_DEPENDS pgslib_mpi)
  else()
    list(APPEND Truchas_DEPENDS pgslib_serial)
  endif()  
endif()

# UbikSolve
if (ENABLE_UbikSolve)
  list(APPEND Truchas_DEPENDS ubiksolve)
endif()  

# Chaco
if ( ENABLE_Chaco )
  list(APPEND Truchas_DEPENDS chaco)
endif()  

# Chaparral
if (ENABLE_Chaparral)
  list(APPEND Truchas_DEPENDS chaparral)
endif()  

# Danu
if (ENABLE_Danu)
  list(APPEND Truchas_DEPENDS danu_f)
endif()

# Various build targets
list(APPEND Truchas_DEPENDS 
            ${NetCDF_BUILD_TARGET}
	    ${HYPRE_BUILD_TARGETS})

# Force the completion of the targets defined 
# in Truchas_DEPENDS before compiling source for
# this library
add_dependencies(truchas ${Truchas_DEPENDS})

# - Link libraries
set(Truchas_LINK_LIBS
            ${NetCDF_Fortran_LIBRARIES}
	    ${HYPRE_LIBRARY})

if(ENABLE_PGSLib)
  list(APPEND Truchas_LINK_LIBS pgslib_fe)
  if (ENABLE_MPI)
    list(APPEND Truchas_LINK_LIBS pgslib_mpi)
  else()
    list(APPEND Truchas_LINK_LIBS pgslib_serial)
  endif()
endif()

if (ENABLE_UbikSolve)
  list(APPEND Truchas_LINK_LIBS ubiksolve)
endif() 

if (ENABLE_Chaco)
  list(APPEND Truchas_LINK_LIBS chaco)
endif()

if (ENABLE_Chaparral)
  list(APPEND Truchas_LINK_LIBS vflib)
endif()

if (ENABLE_Danu)
  list(APPEND Truchas_LINK_LIBS danu_f)
endif()

# To pick up the dlopen and dlclose symbols
list(APPEND Truchas_LINK_LIBS ${CMAKE_DL_LIBS})

target_link_libraries(truchas ${Truchas_LINK_LIBS})

# The linking of external packages is broken. Binaries
# and libraries that link to truchas library should not
# need to explicitly add the Truchas_LINK_LIBS. 
# Push this variable definition up one level as
# a work-around until this is fixed.
set(Truchas_LINK_LIBS ${Truchas_LINK_LIBS} PARENT_SCOPE)
       
# --- Main executable

# Process the files
fortran_preprocess_files(MAIN_SOURCE_FILES
                         FILES drivers/main.F90
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${Truchas_FPP_FLAGS}
			 PROCESS_TARGET ProcessTruchasMainFiles)

set(fc_flags -I${NetCDF_INCLUDE_DIR})
if(ENABLE_PGSLib)
  list(APPEND fc_flags -I${PGSLib_MODULE_DIR})
endif()
if(ENABLE_UbikSolve)
  list(APPEND fc_flags -I${UbikSolve_MODULE_DIR})
endif()
if(ENABLE_Danu)
  list(APPEND fc_flags -I${Danu_Fortran_MODULE_DIR})
endif()
build_whitespace_string(MAIN_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${MAIN_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${MAIN_COMPILE_FLAGS})

# Add the executable with required dependencies and libraries		       
add_executable(truchas_exe ${MAIN_SOURCE_FILES})
add_dependencies(truchas_exe ProcessTruchasMainFiles)
target_link_libraries(truchas_exe truchas ${Truchas_LINK_LIBS})

# Kep binary name consistent with old naming convention 
set(truchas_exe_name t-${Truchas_BUILD_ID}-${Truchas_VERSION})
set_target_properties(truchas_exe PROPERTIES
                      OUTPUT_NAME ${truchas_exe_name}
		      EchoString "Built Truchas executable ${truchas_exe_name}")

# --- Install rules 		    
install(TARGETS truchas truchas_exe
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
