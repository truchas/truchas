#
#  Verify HDF5 installations
#
#  File is intended to be included in external/CMakeLists.txt file
#  This file will set HDF5_VERIFIED to TRUE if HDF5 is located 
#  and it meets the requirements.

#
# HDF5 Requirements
#
#  o Version > 1.8
#  o Need the hdf5_hl (high-level library)
# 

# FindHDF5 relies on HDF5_INSTALL_PREFIX to locate 
# HDF5
if (NOT HDF5_INSTALL_PREFIX)
  set(HDF5_INSTALL_PREFIX ${TruchasExternal_INSTALL_PREFIX})
endif()

# Boolean Evaluator
include(BoolEval)

# Default is HDF5_VERIFIED False
set(HDF5_VERIFIED False)

# Locate HDF5
find_package(HDF5)

# Verify the package
if (HDF5_FOUND)

  message(STATUS "Verify HDF5 package")

  # Version verification
  bool_eval(HDF5_VERSION_OK NOT ${HDF5_VERSION} VERSION_LESS 1.8.3)

  # Need the high level library
  bool_eval(HDF5_VERIFIED HDF5_HL_LIBRARY AND HDF5_VERSION_OK)

endif(HDF5_FOUND)

if(HDF5_VERIFIED)
  message(STATUS "Verify HDF5 package -- ok")
else(HDF5_VERIFIED)
  message(STATUS "Verify HDF5 package -- failed")
endif(HDF5_VERIFIED)

