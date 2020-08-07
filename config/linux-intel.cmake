# Generic Linux with the Intel Compilers

set(CMAKE_C_COMPILER icc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER ifort CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS_RELEASE "-u -traceback -fpe0 -DNDEBUG" CACHE STRING "Fortran compile flags")

set(CMAKE_C_FLAGS_DEBUG "-traceback" CACHE STRING "C compile flags")
set(CMAKE_Fortran_FLAGS_DEBUG "-u -C -check noarg_temp_created -traceback -fpe0" CACHE STRING "Fortran compile flags")
