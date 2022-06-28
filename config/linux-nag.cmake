# Generic Linux with NAG Fortran and GNU C

set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER nagfor CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS_RELEASE "-u -O3 -DNDEBUG -f2018 -w=uda"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_DEBUG "-u -O0 -C -C=dangling -gline -nan -f2018 -w=uda"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}"
    CACHE STRING "Fortran compile flags")
