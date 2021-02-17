# Generic Linux with GNU Fortran and GNU C

set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(Truchas_Fortran_FLAGS "-fimplicit-none -ffree-line-length-none"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_RELEASE "${Truchas_Fortran_FLAGS} -O3 -DNDEBUG"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_DEBUG "${Truchas_Fortran_FLAGS} \
-g -O0 -fcheck=bits,bounds,do,mem,pointer -finit-real=nan -finit-integer=-2147483647"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}"
    CACHE STRING "Fortran compile flags")
