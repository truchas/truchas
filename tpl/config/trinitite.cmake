# Cray Linux Environment (CLE) -- Trinitite

set(ENABLE_MPI YES CACHE BOOL "Parallel libraries")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")

# Cray compiler wrappers
set(CMAKE_C_COMPILER cc CACHE STRING "C Compiler")
set(CMAKE_CXX_COMPILER CC CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER ftn CACHE STRING "Fortran Compiler")
