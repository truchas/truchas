# Generic Linux with NAG Fortran and GNU C/C++

set(ENABLE_MPI NO CACHE BOOL "Parallel libraries")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_CXX_COMPILER g++ CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER nagfor CACHE STRING "Fortran Compiler")

set(CMAKE_Fortran_FLAGS "-u -O3" CACHE STRING "Fortran compile flags")
