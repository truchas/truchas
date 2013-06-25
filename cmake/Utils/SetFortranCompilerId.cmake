#
# Usage:
#
#  set_fortran_compiler_id()
#
# Called from any CMakeLists.txt file, this macro will set the 
# CMAKE_Fortran_COMPILER_ID variable if it is not already set.
#
#
MACRO(SET_FORTRAN_COMPILER_ID)

  if (CMAKE_Fortran_COMPILER_ID) 
    # Do nothing
  else() 
    if ( "${CMAKE_Fortran_COMPILER}" MATCHES "gfortran" )
      set(CMAKE_Fortran_COMPILER_ID GNU)
    elseif ( "${CMAKE_Fortran_COMPILER}" MATCHES "ifort")
      set(CMAKE_Fortran_COMPILER_ID Intel)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "nagfor")
      set(CMAKE_Fortran_COMPILER_ID Nag)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "pgf77|pgf90|pgf95")
      set(CMAKE_Fortran_COMPILER_ID PGI)
    else() 
      message(WARNING "${CMAKE_Fortran_COMPILER} does not match known Fortran types.")
    endif()
  endif()  

ENDMACRO(SET_FORTRAN_COMPILER_ID)
