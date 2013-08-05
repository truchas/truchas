#
#  Verify ZLIB installations
#
#  File is intended to be included in external/CMakeLists.txt file
#  This file will set ZLIB_VERIFIED to TRUE if ZLIB is located 
#  and it meets the requirements.

#
# ZLIB Requirements
#
#  o Version >= 1.2.6 

# Boolean evaluator
include(BoolEval)

# Locate ZLIB
find_package(ZLIB)

# Verify the package
if (ZLIB_FOUND)

  message(STATUS "Verify ZLIB package")

  # Version verification
  bool_eval(ZLIB_VERIFIED 
            NOT "${ZLIB_VERSION_STRING}" VERSION_LESS "1.2")

endif(ZLIB_FOUND)

if(ZLIB_VERIFIED)
  message(STATUS "Verify ZLIB package -- ok")
else(ZLIB_VERIFIED)
  message(STATUS "Verify ZLIB package -- failed")
endif(ZLIB_VERIFIED)

