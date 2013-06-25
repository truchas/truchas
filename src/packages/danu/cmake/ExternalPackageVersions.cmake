# ############################################################################ #
#
#   DANU - External Package Requirements  
#  
# ############################################################################ #


# ---------------------------------------------------------------------------- #
# ZLIB
# ---------------------------------------------------------------------------- #
set(ZLIB_VERSION_MAJOR 1)
set(ZLIB_VERSION_MINOR 2)
set(ZLIB_VERSION_PATCH 8)
set(ZLIB_VERSION ${ZLIB_VERSION_MAJOR}.${ZLIB_VERSION_MINOR}.${ZLIB_VERSION_PATCH})
set(ZLIB_ARCHIVE_FILE zlib-${ZLIB_VERSION}.tar.gz)
set(ZLIB_MD5SUM 44d667c142d7cda120332623eab69f40)  
set(ZLIB_DEPENDENCIES)

# ---------------------------------------------------------------------------- #
# HDF5
# ---------------------------------------------------------------------------- #
set(HDF5_VERSION_MAJOR 1)
set(HDF5_VERSION_MINOR 8)
set(HDF5_VERSION_PATCH 8)
set(HDF5_VERSION ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}.${HDF5_VERSION_PATCH})
set(HDF5_ARCHIVE_FILE hdf5-${HDF5_VERSION}.tar.gz)
set(HDF5_MD5SUM  1196e668f5592bfb50d1de162eb16cff)  
set(HDF5_DEPENDENCIES ZLIB)

# ---------------------------------------------------------------------------- #
# SWIG
# ---------------------------------------------------------------------------- #
set(SWIG_VERSION_MAJOR 2)
set(SWIG_VERSION_MINOR 0)
set(SWIG_VERSION_PATCH 4)
set(SWIG_VERSION ${SWIG_VERSION_MAJOR}.${SWIG_VERSION_MINOR}.${SWIG_VERSION_PATCH})
set(SWIG_ARCHIVE_FILE swig-${SWIG_VERSION}.tar.gz)
set(SWIG_MD5SUM  4319c503ee3a13d2a53be9d828c3adc0)  
set(SWIG_DEPENDENCIES)

# ---------------------------------------------------------------------------- #
# Check - C Unit Test Framework 
# ---------------------------------------------------------------------------- #
set(Check_VERSION_MAJOR 0)
set(Check_VERSION_MINOR 9)
set(Check_VERSION_PATCH 8)
set(Check_VERSION ${Check_VERSION_MAJOR}.${Check_VERSION_MINOR}.${Check_VERSION_PATCH})
set(Check_ARCHIVE_FILE check-${Check_VERSION}.tar.gz)
set(Check_MD5SUM  5d75e9a6027cde79d2c339ef261e7470)  
set(Check_DEPENDENCIES)


