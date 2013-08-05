#
#  Verify HYPRE installations
#
#  File is intended to be included in external/CMakeLists.txt file
#  This file will set HYPRE_VERIFIED to TRUE if HYPRE is located 
#  and it meets the requirements.

#
# HYPRE Requirements
#
#  o Version != 2.6.0b
#  o HYPRE_IS_PARALLEL == ENABLE_MPI flag

# Boolean evaluator
include(BoolEval)

# Default is HYPRE_VERIFIED False
set(HYPRE_VERIFIED False)

# Install prefix
if ( NOT HYPRE_INSTALL_PREFIX )
  set(HYPRE_INSTALL_PREFIX ${TruchasExternal_INSTALL_PREFIX})
endif()  

# Locate HYPRE
set(HYPRE_REQUIRED_VERSION 2.6.0b)
find_package(HYPRE)

# Verify the package
if (HYPRE_FOUND)

  message(STATUS "Verify HYPRE package")

  # Version check
  bool_eval(HYPRE_VERSION_OK 
            ${HYPRE_VERSION} VERSION_EQUAL ${HYPRE_REQUIRED_VERSION})

  # Parallel (MPI) check
  bool_eval(HYPRE_PARALLEL_OK NOT ENABLE_MPI EQUAL HYPRE_IS_PARALLEL)

  # Both checks must pass
  bool_eval(HYPRE_VERIFIED HYPRE_VERSION_OK AND HYPRE_PARALLEL_OK)

endif(HYPRE_FOUND)

if(HYPRE_VERIFIED)
  message(STATUS "Verify HYPRE package -- ok")
else(HYPRE_VERIFIED)
  message(STATUS "Verify HYPRE package -- failed ")
endif(HYPRE_VERIFIED)

