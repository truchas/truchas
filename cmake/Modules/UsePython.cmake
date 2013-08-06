# ############################################################################ #
# 
# Useful Python macros/functions
#
# PYTHON_EXECUTABLE must be defined!
#
# Included MACROS
#
# Byte compile a python module
# PYTHON_COMPILE_MODULE(MODULE OUTFILE_VARIABLE)
#
# 
# 
# ############################################################################ #
MACRO(PYTHON_COMPILE_MODULE mod out)
 
  get_filename_component(_mod_path ${mod} PATH)
  get_filename_component(_mod_fullname ${mod} ABSOLUTE)
  get_filename_component(_mod_name ${mod} NAME)
  get_filename_component(_mod_ext ${mod} EXT)

  execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
                  "import py_compile; py_compile.compile(${_mod_fullname})" 
		  RESULT_VARIBLE _mod_result
		  ERROR_QUIET)
  if ( NOT "${_mod_result}" )
    message(FATAL_ERROR "Failed to byte compile Python module ${mod} ")
  endif()

  set(${out} ${CMAKE_CURRENT_BINARY_DIR}/${_mod_name}c)

ENDMACRO(PYTHON_COMPILE_MODULE mod out)




