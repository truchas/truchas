# Generic Linux with GNU C/C++/Fortran

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_CXX_COMPILER g++ CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS "-fimplicit-none" CACHE STRING "Fortran compile flags")
