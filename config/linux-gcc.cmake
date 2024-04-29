# Generic Linux with GNU Fortran and GNU C

set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -ffpe-trap=zero,invalid,overflow"
    CACHE STRING "Fortran compile flags")
# -fcheck=bits is excluded here since it is not supported on GCC 9.3.0
set(CMAKE_Fortran_FLAGS_DEBUG
  "-g -O0 -ffpe-trap=zero,invalid,overflow -fcheck=bounds,do,mem,pointer -finit-real=nan -finit-integer=-2147483647"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}"
    CACHE STRING "Fortran compile flags")
