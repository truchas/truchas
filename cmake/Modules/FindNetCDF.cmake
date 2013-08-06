# -*- mode: cmake -*-
# Usage:
#    Control the search through NETCDF_INSTALL_PREFIX or setting environment variable
#    NETCDF_ROOT to the NETCDF installation prefix.
#    By default only searches for C library. To search for the C++ library
#    set 'CXX' in the COMPONENTS option in find_package(NETCDF) 
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    NETCDF_FOUND            (BOOL)       Flag indicating if NETCDF was found
#    NETCDF_INCLUDE_DIR      (PATH)       Path to the NETCDF include file
#    NETCDF_INCLUDE_DIRS     (LIST)       List of all required include files
#    NETCDF_C_LIBRARIES      (LIST)       List of all required libraries to link to 
#                                          the NETCDF C library
#    NETCDF_CXX_LIBRARIES    (LIST)       List of all required libraries to link to 
#                                          the NETCDF C++ library (If CXX is in COMPONENTS)
#
#    Additional variables set
#    NETCDF_C_LIBRARY        (FILE)       NETCDF C library
#    NETCDF_CXX_LIBRARY      (FILE)       NETCDF C++ library (If CXX is in the COMPONENTS)
#    NETCDF_LARGE_DIMS       (BOOL)       Checks the header files for size of 
#                                          NC_MAX_DIMS, NC_MAX_VARS and NC_MAX_VARS_DIMS
#                                          Returns TRUE if
#                                          NC_MAX_DIMS >= 655363
#                                          NC_MAX_VARS >= 524288
#                                          NC_MAX_VAR_DIMS >= 8
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# ######################################
# Begin local macros/functions
#
# Macro to handle print out
macro(_netcdf_status mess)
  if(NOT NETCDF_FIND_QUIETLY)
    message(STATUS "${mess}")
  endif()
endmacro()
#
# End local macros/functions
# ######################################


# Search paths
set(_netcdf_search_paths
     ${NETCDF_INSTALL_PREFIX}
     ENV NETCDF_ROOT)

# Locate the config binary
# Search for the nc-config binary
find_program(NETCDF_CONFIG_BINARY nc-config
             HINTS ${NETCDF_INSTALL_PREFIX} ${NETCDF_BIN_DIR}
	     PATH_SUFFIXES bin Bin
	     DOC "NETCDF configuration script")

# Determine if HDF5 is present           
set(NETCDF_HAS_HDF5 False)
if(NETCDF_CONFIG_BINARY)

  execute_process(COMMAND ${NETCDF_CONFIG_BINARY} --has-hdf5
                  RESULT_VARIABLE result
                  OUTPUT_VARIABLE out
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

 if (NOT ${result})              
   if ( "${out}" MATCHES "yes" )
     set(NETCDF_HAS_HDF5 True)
   else("${out}" MATCHES "yes")
     set(NETCDF_HAS_HDF5 False)
   endif("${out}" MATCHES "yes")
 endif(NOT ${result})  

endif(NETCDF_CONFIG_BINARY)  

# Locate a HDF5 installation
if(NETCDF_HAS_HDF5)
  find_package(HDF5)
endif()

# Include file search
if ( NOT NETCDF_INCLUDE_DIRS )

  if(NETCDF_CONFIG_BINARY)

    execute_process(COMMAND ${NETCDF_CONFIG_BINARY} --includedir
                    RESULT_VARIABLE result
                    OUTPUT_VARIABLE dir
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT ${result})              
      set(NETCDF_INCLUDE_DIR ${dir})
    endif(NOT ${result})

  endif(NETCDF_CONFIG_BINARY) 

  find_path(NETCDF_INCLUDE_DIR
            NAMES netcdf.h
            HINTS ${_netcdf_search_paths}
            PATH_SUFFIXES include)

  list(APPEND NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})       

  # Add HDF5 if present
  if( NETCDF_HAS_HDF5)
    if(HDF5_FOUND)
      list(APPEND NETCDF_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    endif(HDF5_FOUND)
  endif(NETCDF_HAS_HDF5)    

endif(NOT NETCDF_INCLUDE_DIRS)

# NETCDF_LIBRARY (NETCDF_C_LIBRARY)
find_library(NETCDF_C_LIBRARY
             NAMES netcdf
             HINTS ${_netcdf_search_paths} ${NETCDF_LIB_DIR}
             PATH_SUFFIXES lib Lib)
set(NETCDF_LIBRARY ${NETCDF_C_LIBRARY})

# NETCDF_CXX_LIBRARY
find_library(NETCDF_CXX_LIBRARY
             NAMES netcdf_c++ 
	     HINTS ${_netcdf_search_paths} ${NETCDF_LIB_DIR}
	     PATH_SUFFIXES lib Lib)

# NETCDF_Fortran_LIBRARY

find_library(NETCDF_Fortran_LIBRARY
             NAMES netcdff 
	     HINTS ${_netcdf_search_paths} ${NETCDF_LIB_DIR}
	     PATH_SUFFIXES lib Lib)

# NETCDF_C_LIBRARIES
if(NETCDF_C_LIBRARY)
  set(NETCDF_C_LIBRARIES ${NETCDF_C_LIBRARY})
  if ( NETCDF_HAS_HDF5 AND HDF5_FOUND )
    list(APPEND NETCDF_C_LIBRARIES ${HDF5_C_LIBRARY} ${HDF5_LINK_LIBRARIES})
  endif()
else()
  set(NETCDF_C_LIBRARIES NETCDF_C_LIBRARIES-NOTFOUND)
endif()  

# NETCDF_CXX_LIBRARIES
if(NETCDF_CXX_LIBRARY)
  set(NETCDF_CXX_LIBRARIES ${NETCDF_CXX_LIBRARY} ${NETCDF_C_LIBRARIES})
else()  
  set(NETCDF_CXX_LIBRARIES NETCDF_CXX_LIBRARIES-NOTFOUND)
endif()  

# NETCDF_Fortran_LIBRARIES
if(NETCDF_Fortran_LIBRARY)
  set(NETCDF_Fortran_LIBRARIES ${NETCDF_Fortran_LIBRARY} ${NETCDF_C_LIBRARIES})
else()  
  set(NETCDF_Fortran_LIBRARIES NETCDF_Fortran_LIBRARIES-NOTFOUND)
endif()  

# NETCDF Version
if(NETCDF_CONFIG_BINARY)
  execute_process(COMMAND "${NETCDF_CONFIG_BINARY}" "--version"
                  RESULT_VARIABLE _ret_code
                  OUTPUT_VARIABLE _stdout
                  ERROR_VARIABLE  _stderr)
  
  string(REGEX REPLACE "[\n\r ]" "" _version_answer "${_stdout}")
  string(REGEX REPLACE "netCDF" "" NETCDF_VERSION "${_version_answer}")
else()  
  set(NETCDF_VERSION NETCDF-NOTFOUND)
endif()   

# Large dimension check
if ( NETCDF_INCLUDE_DIR ) 
       
  set(netcdf_h "${NETCDF_INCLUDE_DIR}/netcdf.h" )
  if ( EXISTS ${netcdf_h} )
    #_netcdf_status( "NetCDF include file ${netcdf_h} will be searched for define values")

    file(STRINGS "${netcdf_h}" netcdf_max_dims_string REGEX "^#define NC_MAX_DIMS")
    string(REGEX REPLACE "[^0-9]" "" netcdf_max_dims "${netcdf_max_dims_string}")

    file(STRINGS "${netcdf_h}" netcdf_max_vars_string REGEX "^#define NC_MAX_VARS")
    string(REGEX REPLACE "[^0-9]" "" netcdf_max_vars "${netcdf_max_vars_string}")

    file(STRINGS "${netcdf_h}" netcdf_max_var_dims_string REGEX "^#define NC_MAX_VAR_DIMS")
    string(REGEX REPLACE "[^0-9]" "" netcdf_max_var_dims "${netcdf_max_var_dims_string}")


    if ( 
         ( (netcdf_max_dims EQUAL 65536)  OR (netcdf_max_dims GREATER 65536)  ) AND
         ( (netcdf_max_vars EQUAL 524288) OR (netcdf_max_vars GREATER 524288) ) AND
         ( (netcdf_max_var_dims EQUAL 8)  OR  (netcdf_max_var_dims GREATER 8) )

       )
         set(NETCDF_LARGE_DIMS TRUE)
    else()
         message(WARNING "The NetCDF found in ${NetCDF_INSTALL_PREFIX} does not have the correct NC_MAX_DIMS, NC_MAX_VARS and NC_MAX_VAR_DIMS\n"
                         "It may not be compatible with other TPL libraries such MOAB and ExodusII\n" )
       set(NETCDF_LARGE_DIMS FALSE)
    endif()

  else() 

    set(NETCDF_LARGE_DIMS NETCDF_LARGE_DIMS-NOTFOUND)

  endif()  
else()
  set(NETCDF_LARGE_DIMS NETCDF_LARGE_DIMS-NOTFOUND)
endif()    

# Now set the component flags need NetCDF_<lang>_FOUND
# to trigger the HANDLE_COMPONENTS correctly
foreach(lang C CXX Fortran)
  set(lib_var NETCDF_${lang}_LIBRARY)
  if(${lib_var})
    set(NetCDF_${lang}_FOUND True)
    set(NETCDF_${lang}_FOUND True)
  else()  
    set(NetCDF_${lang}_FOUND False)
    set(NETCDF_${lang}_FOUND False)
  endif()
endforeach()  

find_package_handle_standard_args(NetCDF
                                  REQUIRED_VARS 
				      NETCDF_INCLUDE_DIR 
				      NETCDF_LIBRARY
				      NETCDF_INCLUDE_DIRS
				      NETCDF_C_LIBRARIES
				  VERSION_VAR NETCDF_VERSION
				  HANDLE_COMPONENTS)



                  




