#
#  Verify PETACA installations
#
#  File is intended to be included in external/CMakeLists.txt file
#  This file will set PETACA_VERIFIED to TRUE if PETACA is located 
#  and it meets the requirements.

#
# PETACA Requirements
#

# FindPETACA relies on PETACA_INSTALL_PREFIX to locate 
# PETACA
if (NOT PETACA_INSTALL_PREFIX)
  set(PETACA_INSTALL_PREFIX ${TruchasExternal_INSTALL_PREFIX})
endif()

# FindYAJL relies on YAJL_INSTALL_PREFIX to locate 
# YAJL and Petaca require YAJL
if (NOT YAJL_INSTALL_PREFIX)
  set(YAJL_INSTALL_PREFIX ${TruchasExternal_INSTALL_PREFIX})
endif()


# Boolean evaluator
include(BoolEval)

# Locate PETACA
find_package(PETACA)

# Verify the package
if (PETACA_FOUND)

  message(STATUS "Verify PETACA package")
  set(PETACA_VERIFIED TRUE)

endif(PETACA_FOUND)

if(PETACA_VERIFIED)
  message(STATUS "Verify PETACA package -- ok")
else(PETACA_VERIFIED)
  message(STATUS "Verify PETACA package -- failed")
endif(PETACA_VERIFIED)

