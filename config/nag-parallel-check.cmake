##############################################################################
#                                                                             
# Truchas build configuration: Nag compiler, parallel, debug with runtime checking
#                                                                             
##############################################################################

# Compilers
set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_CXX_COMPILER g++ CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER nagfor CACHE STRING "Fortran Compiler")

# Compiler flags.  The C/C++ flags default to acceptable CMake-defined
# built-ins that depend on the compiler ID and build type.  Not so for
# Fortran Flags -- we must specify them (or get none).
set(CMAKE_Fortran_FLAGS "-u -C -gline -nan -maxcontin=99" CACHE STRING "Fortran compile flags")

# Optimized build
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")

# Serial executable
set(ENABLE_MPI True CACHE BOOL "MPI build control flag")
