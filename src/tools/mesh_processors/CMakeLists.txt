project(TruchasToolsMesh Fortran)

# --- CMake modules
include(FortranPreprocessFiles)

# --- Include paths

# Want program modules in their respective build directories.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})

# Add the Danu module directory
include_directories(${Danu_Fortran_MODULE_DIR})

# --- Preprocessing flags
set(fpp_flags_list 
    ${Truchas_FPP_FLAGS} -I${TruchasToolsMesh_SOURCE_DIR})
list(APPEND fpp_flags_list -I${Danu_Fortran_MODULE_DIR})
build_whitespace_string(fpp_flags ${fpp_flags_list})

# --- Add subdirectories 

# addgaps 
add_subdirectory(addgaps)

# scramble 
add_subdirectory(scramble)

# --- Install rules 		    
