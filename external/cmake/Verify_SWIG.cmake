#
#  Verify SWIG installations
#
#  File is intended to be included in external/CMakeLists.txt file
#  This file will set SWIG_VERIFIED to TRUE if SWIG is located 
#  and it meets the requirements.

#
# SWIG Requirements
#
#  o Version > 2.0

# Boolean evaluator
include(BoolEval)

# Default is SWIG_VERIFIED False
set(SWIG_VERIFIED False)

# Locate SWIG
find_package(SWIG)

# Verify the package
if (SWIG_FOUND)

  message(STATUS "Verify SWIG package")

  # Version verification
  bool_eval(SWIG_VERIFIED "${SWIG_VERSION}" VERSION_GREATER "2.0")

endif(SWIG_FOUND)

if(SWIG_VERIFIED)
  message(STATUS "Verify SWIG package -- ok")
else(SWIG_VERIFIED)
  message(STATUS "Verify SWIG package -- failed")
endif(SWIG_VERIFIED)

