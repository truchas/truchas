##############################################################################
#                                                                             
# Truchas build configuration: Intel compilers, serial, debug with runtime checking
#                                                                             
##############################################################################

# Compilers
set(CMAKE_C_COMPILER icc CACHE STRING "C Compiler")
set(CMAKE_CXX_COMPILER icpc CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER ifort CACHE STRING "Fortran Compiler")

# Compiler flags.  The C/C++ flags default to acceptable CMake-defined
# built-ins that depend on the compiler ID and build type.  Not so for
# Fortran Flags -- we must specify them (or get none).
set(CMAKE_Fortran_FLAGS "-u -C -traceback" CACHE STRING "Fortran compile flags")
set(CMAKE_C_FLAGS "-traceback" CACHE STRING "C compile flags")

# Optimized build
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")

# Serial executable
set(ENABLE_MPI False CACHE BOOL "MPI build control flag")
