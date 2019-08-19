# Generic Linux with the Intel Compilers

set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")

set(CMAKE_C_COMPILER icc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER ifort CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_C_FLAGS "-traceback" CACHE STRING "C compile flags")
set(CMAKE_Fortran_FLAGS "-u -C -init=snan -check noarg_temp_created -traceback -fpe0" CACHE STRING "Fortran compile flags")
