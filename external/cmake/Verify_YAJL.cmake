#
#  Verify YAJL installations
#
#  File is intended to be included in external/CMakeLists.txt file
#  This file will set YAJL_VERIFIED to TRUE if YAJL is located 
#  and it meets the requirements.

#
# YAJL Requirements
#
#  o Version >= 2.0.4

# Boolean evaluator
include(BoolEval)

# Locate YAJL
find_package(YAJL)

# Verify the package
if (YAJL_FOUND)

  message(STATUS "Verify YAJL package")

  # Version verification
  bool_eval(YAJL_VERIFIED 
            NOT "${YAJL_VERSION}" VERSION_LESS "2.0.4")

endif(YAJL_FOUND)

if(YAJL_VERIFIED)
  message(STATUS "Verify YAJL package -- ok")
else(YAJL_VERIFIED)
  message(STATUS "Verify YAJL package -- failed")
endif(YAJL_VERIFIED)

