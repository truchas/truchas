# Generic Linux with NAG Fortran and GNU C

set(CMAKE_C_COMPILER $ENV{CC} CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER $ENV{FC} CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(Truchas_Fortran_FLAGS "-u -f2018 -w=uda -target=native")
set(CMAKE_Fortran_FLAGS_RELEASE "${Truchas_Fortran_FLAGS} -O3 -DNDEBUG"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_DEBUG "${Truchas_Fortran_FLAGS} -O0 -C -C=dangling -gline -nan"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}"
    CACHE STRING "Fortran compile flags")
