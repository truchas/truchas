# Generic Linux with LLVM flang and clang compilers

set(CMAKE_C_COMPILER $ENV{CC} CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER $ENV{FC} CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}"
    CACHE STRING "Fortran compile flags")
