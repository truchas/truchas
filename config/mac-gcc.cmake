# Generic MacOS with GNU Fortran and GNU C

set(CMAKE_C_COMPILER $ENV{CC} CACHE STRING "C Compiler")
set(CMAKE_CXX_COMPILER $ENV{CXX} CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER $ENV{FC} CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -DNDEBUG"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_DEBUG
    "-g -O0 -fcheck=bounds,do,mem,pointer -finit-real=nan -finit-integer=-2147483647"
    CACHE STRING "Fortran compile flags")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}"
    CACHE STRING "Fortran compile flags")

# Mac has a different default linking model compared to Linux.
# PGSLib depends heavily on the linux model requireing -undefined dynamic_lookup
# -headerpad_max_install_names allow for the shared libraries to be installed
set(CMAKE_SHARED_LINKER_FLAGS
  "-Wl,-undefined -Wl,dynamic_lookup -Wl,-headerpad_max_install_names"
  CACHE STRING "Mac linker flags")
