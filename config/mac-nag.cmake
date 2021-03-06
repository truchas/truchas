# Generic Linux with NAG Fortran and GNU C

set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER nagfor CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS_DEBUG "-u -C -C=dangling -gline -nan" CACHE STRING "Fortran compile flags")

set(CMAKE_Fortran_FLAGS_RELEASE "-u -O2 -DNDEBUG" CACHE STRING "Fortran compile flags")

# Mac has a different default linking model compared to Linux.
# PGSLib depends heavily on the linux model requireing -undefined dynamic_lookup
# -headerpad_max_install_names allow for the shared libraries to be installed
set(CMAKE_SHARED_LINKER_FLAGS
  "-Wl,-undefined -Wl,dynamic_lookup -Wl,-headerpad_max_install_names"
  CACHE STRING "Mac linker flags")
