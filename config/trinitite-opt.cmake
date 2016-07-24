# Cray Linux Environment (CLE) -- Trinitite

set(ENABLE_MPI ON CACHE BOOL "Parallel Truchas")
set(ENABLE_SHARED OFF CACHE BOOL "Build shared libraries")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

# Cray compiler wrappers
set(CMAKE_C_COMPILER cc CACHE STRING "C Compiler")
set(CMAKE_CXX_COMPILER CC CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER ftn CACHE STRING "Fortran Compiler")

# Additional flags to the default CMAKE_<lang>_FLAGS_<build_type> flags
set(CMAKE_Fortran_FLAGS "-u" CACHE STRING "Fortran compile flags")

# These are some really specific hacks for the Cray, and especially for
# the current configuration of Trinitite.  They will likely need updating.

set(CRAY True CACHE BOOL "Cray Linux Environment")
set(ENV{CRAYPE_LINK_TYPE} "dynamic")

# Find_package(MPI) is not able to set these values properly, so hardwire them.
set(MPIEXEC "/opt/cray/alps/6.1.4-18.1/bin/aprun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

# Hardwire path to the numpy include directory. The system python2.7 has
# numpy 1.8 but is missing the header files. A separate "anaconda" python2.7
# module includes a proper numpy 1.9 installation, but is unusable because
# of final link problems due to its own libmpich.a.  So we just use its
# header file (different version!) but the system libraries.
set(NUMPY_INCLUDE_DIRS "/usr/projects/hpcsoft/cle6.0/common/anaconda/2.1.0-python-2.7/lib/python2.7/site-packages/numpy/core/include/" CACHE STRING "NumPy directory on Cray")
