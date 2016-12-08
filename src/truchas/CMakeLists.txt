project(TruchasExe C Fortran)

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

include(bcs/TruchasFiles.cmake)
include(communication/TruchasFiles.cmake)
include(discrete_operators/TruchasFiles.cmake)
include(distributed_mesh/TruchasFiles.cmake)
include(drivers/TruchasFiles.cmake)
include(exodus/TruchasFiles.cmake)
include(nka/TruchasFiles.cmake)
include(functions/TruchasFiles.cmake)
include(grid_mapping/TruchasFiles.cmake)
include(input/TruchasFiles.cmake)
include(materials/TruchasFiles.cmake)
include(ustruc/TruchasFiles.cmake)
include(ode/TruchasFiles.cmake)
include(legacy_mesh/TruchasFiles.cmake)
include(output/TruchasFiles.cmake)
include(partitioning/TruchasFiles.cmake)
include(physics/TruchasFiles.cmake)
include(setup/TruchasFiles.cmake)
include(solver/TruchasFiles.cmake)
include(utilities/TruchasFiles.cmake)
include(toolpath/TruchasFiles.cmake)

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
list(APPEND Truchas_DEPENDS pgslib_fe)
if (ENABLE_MPI)
  list(APPEND Truchas_DEPENDS pgslib_mpi)
else()
  list(APPEND Truchas_DEPENDS pgslib_serial)
endif()  

list(APPEND Truchas_DEPENDS ubiksolve)
list(APPEND Truchas_DEPENDS chaco)
list(APPEND Truchas_DEPENDS danu_f)

# Force the completion of the targets defined 
# in Truchas_DEPENDS before compiling source for
# this library
add_dependencies(truchas ${Truchas_DEPENDS})

# - Link libraries
set(Truchas_LINK_LIBS exodus hypre petaca)

list(APPEND Truchas_LINK_LIBS pgslib_fe)
if (ENABLE_MPI)
  list(APPEND Truchas_LINK_LIBS pgslib_mpi)
else()
  list(APPEND Truchas_LINK_LIBS pgslib_serial)
endif()

list(APPEND Truchas_LINK_LIBS ubiksolve)
list(APPEND Truchas_LINK_LIBS chaco)
list(APPEND Truchas_LINK_LIBS danu_f)

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

set(fc_flags -I${Danu_Fortran_MODULE_DIR})
list(APPEND fc_flags -I${PGSLib_MODULE_DIR})
list(APPEND fc_flags -I${UbikSolve_MODULE_DIR})
if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  list(APPEND fc_flags -C)
endif()
build_whitespace_string(MAIN_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${MAIN_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${MAIN_COMPILE_FLAGS})

# Add the executable with required dependencies and libraries         
add_executable(truchas_exe ${MAIN_SOURCE_FILES})
add_dependencies(truchas_exe ProcessTruchasMainFiles)
target_link_libraries(truchas_exe truchas ${Truchas_LINK_LIBS})

# Keep binary name consistent with old naming convention 
set(truchas_exe_name t-${Truchas_BUILD_ID}-${Truchas_VERSION})
set_target_properties(truchas_exe PROPERTIES
                      OUTPUT_NAME ${truchas_exe_name}
                      EchoString "Built Truchas executable ${truchas_exe_name}")

# --- Install rules    
install(TARGETS truchas truchas_exe
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
        ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})

# --- Unit tests
add_subdirectory(solver/Test)
add_subdirectory(physics/heat_species_transport/Test)
add_subdirectory(discrete_operators/Test)
add_subdirectory(utilities/Test)
add_subdirectory(toolpath/Test)
