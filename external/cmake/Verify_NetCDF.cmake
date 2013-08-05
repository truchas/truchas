#
#  Verify NETCDF installations
#
#  File is intended to be included in external/CMakeLists.txt file
#  This file will set NETCDF_VERIFIED to TRUE if NETCDF is located 
#  and it meets the requirements.

#
# NETCDF Requirements
#
#  o Version > 1.8
#  o Need the hdf5_hl (high-level library)
# 

# FindNETCDF relies on NETCDF_INSTALL_PREFIX to locate 
# NETCDF
if (NOT NETCDF_INSTALL_PREFIX)
  set(NETCDF_INSTALL_PREFIX ${TruchasExternal_INSTALL_PREFIX})
endif()

# Boolan Evaluator
include(BoolEval)

# Default is NETCDF_VERIFIED False
set(NETCDF_VERIFIED False)

# Locate NETCDF
find_package(NetCDF COMPONENTS C Fortran)

# Verify the package
if (NETCDF_FOUND)

  message(STATUS "Verify NetCDF package")

  # Version verification
  bool_eval(NETCDF_VERSION_OK NOT ${NETCDF_VERSION} VERSION_LESS 4.1.3)

  bool_eval(NETCDF_VERIFIED ${NETCDF_VERSION_OK} AND ${NETCDF_LARGE_DIMS})

endif(NETCDF_FOUND)

#
set(NetCDF_VERIFIED ${NETCDF_VERIFIED})

if(NETCDF_VERIFIED)
  message(STATUS "Verify NetCDF package -- ok")
else(NETCDF_VERIFIED)
  message(STATUS "Verify NetCDF package -- failed")
endif(NETCDF_VERIFIED)

