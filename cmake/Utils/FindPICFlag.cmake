# - 
#
# Usage
#
#  find_pic_flag(flag_var)
#
# A macro that loops through possible PIC (position in code) flags
# for C and returns the first flag that was successful. Variable set to
# NOTFOUND if all tests fail.
#
MACRO(FIND_PIC_FLAG var)

  include(CheckCCompilerFlag)

  # Possible PIC flags
  set(_test_pic_flags -fpic -fPIC -pic -PIC)

  # Default is NOTFOUND
  set(${var} ${var}-NOTFOUND)

  # Must change result var name for each call 
  # to check_c_compiler_flag
  foreach(_flag ${_test_pic_flags})
    set(res_var flag${_flag}_test)
    check_c_compiler_flag(${_flag} ${res_var})
    if (${${res_var}} )
      set(${var} ${_flag})
      break()
    endif(${${res_var}})
  endforeach()

  if(NOT ${var})
    message(SEND_ERROR "Failed to determine PIC C compiler flag. "
                       "Tried ${_test_pic_flags} and all failed")
  endif()		     

ENDMACRO(FIND_PIC_FLAG)  
