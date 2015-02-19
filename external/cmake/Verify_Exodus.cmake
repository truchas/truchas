#
#  Verify Exodus installations
#
#  File is intended to be included in external/CMakeLists.txt file
#  This file will set EXODUS_VERIFIED to TRUE if Exodus is located 
#  and it meets the requirements.

#
# Exodus Requirements
#
#  o Version >= 514
# FindEXODUS relies on EXODUS_INSTALL_PREFIX to locate 

# User defined EXODUS install prefix
if (NOT EXODUS_INSTALL_PREFIX)
  set(EXODUS_INSTALL_PREFIX ${TruchasExternal_INSTALL_PREFIX})
endif()

# Boolean evaluator
include(BoolEval)

# Locate Exodus
find_package(Exodus)

# Verify the package
if(EXODUS_FOUND)
  message(STATUS "Verify EXODUS package")
  bool_eval(EXODUS_VERIFIED NOT "${EXODUS_VERSION}" VERSION_LESS "514")
endif()

set(Exodus_VERIFIED ${EXODUS_VERIFIED})

if(EXODUS_VERIFIED)
  message(STATUS "Verify EXODUS package -- ok")
else()
  message(STATUS "Verify EXODUS package -- failed")
endif()

