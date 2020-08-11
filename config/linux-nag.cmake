# Generic Linux with NAG Fortran and GNU C

set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER nagfor CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
# NB: -O3 optimization fails with NAG 6.0 (1067) with gcc 4.8.3.
set(CMAKE_Fortran_FLAGS_RELEASE "-u -O2 -DNDEBUG" CACHE STRING "Fortran compile flags")

set(CMAKE_Fortran_FLAGS_DEBUG "-u -C -C=dangling -gline -nan" CACHE STRING "Fortran compile flags")
