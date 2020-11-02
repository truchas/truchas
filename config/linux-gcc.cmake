# Generic Linux with GNU Fortran and GNU C

set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS_RELEASE "-fimplicit-none -ffree-line-length-none -DNDEBUG"
    CACHE STRING "Fortran compile flags")

set(CMAKE_Fortran_FLAGS_DEBUG "-g -fimplicit-none -ffree-line-length-none"
    CACHE STRING "Fortran compile flags")
