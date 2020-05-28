project(TruchasToolsMesh Fortran)

# --- Include paths

# Want program modules in their respective build directories.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})

# --- Add subdirectories

# addgaps
add_subdirectory(addgaps)

# scramble
add_subdirectory(scramble)

# --- Install rules
