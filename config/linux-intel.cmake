# Generic Linux with the Intel Compilers

set(CMAKE_C_COMPILER icc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER ifort CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(Truchas_Fortran_FLAGS="-u -traceback -fpe0")
set(CMAKE_Fortran_FLAGS_RELEASE "${Truchas_Fortran_FLAGS} -O3 -DNDEBUG"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_DEBUG "${Truchas_Fortran_FLAGS} -O0 -C -check noarg_temp_created"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}"
    CACHE STRING "Fortran compile flags")
