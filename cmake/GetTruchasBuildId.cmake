#
# GET_TRUCHAS_BUILD_ID
#
# === Usage ===
#
#  get_truchas_build_id( <var> )
#
# Sets the variable <var> to the Truchas build string
# Follows the same pattern as the olg GNUMake style
# OS NAME.ARCH.COMPILER ID.serial|parallel.opt|debug
#
#
MACRO(GET_TRUCHAS_BUILD_ID id)
  set(_tmp ${CMAKE_SYSTEM_NAME}.${CMAKE_SYSTEM_PROCESSOR}.${CMAKE_Fortran_COMPILER_ID})
  # Want all lower case
  string(TOLOWER "${_tmp}" ${id})
ENDMACRO(GET_TRUCHAS_BUILD_ID id)
