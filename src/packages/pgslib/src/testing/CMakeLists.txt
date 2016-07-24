project(PGSLibTesting Fortran)

message(STATUS "Configuring testing directory")

# --- Enable the default CTest targets (disabled)
#include(CTest)

# --- Scripts to run the tests
set(PGSLib_SERIAL_SCRIPT ${PGSLibTesting_SOURCE_DIR}/run-serial-test)
set(PGSLib_PARALLEL_SCRIPT ${PGSLibTesting_SOURCE_DIR}/run-parallel-test)

# --- Remove the default Fortran module directory
#     Do not want test modules in the same directory as library
#     modules.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# --- Include directories for all tests

# Add the library module directory
include_directories(${PGSLib_MODULE_DIR})

# Add the current source and build directories to the include directories 
set(CMAKE_INCLUDE_CURRENT_DIR True)

# --- Test directories
add_subdirectory(assembly)
add_subdirectory(gath-scatt)
add_subdirectory(performance)
add_subdirectory(tests)


