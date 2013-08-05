#
# BOOL_EVAL(VAR bool_expression)
#  Evaluate bool_expression and set VAR to True|False 
#
# Example:
#  bool_eval(FOO_OK "foo" STREQUAL "foo")
#  sets FOO_OK to True
MACRO(BOOL_EVAL myvar)
  if(${ARGN})
    set(${myvar} True)
  else(${ARGN})
    set(${myvar} False)
  endif()  
ENDMACRO(BOOL_EVAL myvar)  
