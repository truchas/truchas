# ############################################################################ #
#
#
#   DANU CMake 
#    Define compiler flags
#
# ############################################################################ #
include(PrintVariable)
include(DanuGlobalMacros)

# ---------------------------------------------------------------------------- #
# Enable Fortran in CMake
# ---------------------------------------------------------------------------- #
if(ENABLE_Fortran)
  enable_language(Fortran)
endif()

# ---------------------------------------------------------------------------- #
# Set the Fortran Compiler Id
# ---------------------------------------------------------------------------- #
if (ENABLE_Fortran)
  get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  if ( NOT CMAKE_Fortran_COMPILER_ID )

    if(Fortran_COMPILER_NAME MATCHES "gfortran.*$")
      set(CMAKE_Fortran_COMPILER_ID GNU)
    elseif(Fortran_COMPILER_NAME MATCHES "pgf.*$")
      set(CMAKE_Fortran_COMPILER_ID PGI)
    elseif(Fortran_COMPILER_NAME MATCHES "ifort$")  
      set(CMAKE_Fortran_COMPILER_ID Intel)
    elseif(Fortran_COMPILER_NAME MATCHES "nagfor$")  
      set(CMAKE_Fortran_COMPILER_ID NAG)
    else()  
      set(CMAKE_Fortran_COMPILER_ID)
    endif()

  endif()

endif()  

# ---------------------------------------------------------------------------- #
# Fortran Compiler Check
# ---------------------------------------------------------------------------- #
if(ENABLE_Fortran)

  include(FortranCInterface)

  option(ENABLE_FortranCIface_Check ON "Check the Fortran/C interface")
  if ( ENABLE_FortranCIface_Check)
    FortranCInterface_VERIFY()
  endif()

  add_definitions(-DENABLE_FORTRAN)
  include_directories("${Danu_BINARY_DIR}/src")

  # Test the iso_c_binding support
  message(STATUS "Test ${CMAKE_Fortran_COMPILER} for iso_c_binding support")
  set(test_fort_code "${Danu_SOURCE_DIR}/cmake/comp_tests/iso_check.f90")
  try_run(run_result Fortran_COMPILER_SUPPORTS_ISO_C_BINDING
          ${Danu_BINARY_DIR}  ${test_fort_code}
          OUTPUT_VARIABLE out_result)
  if(Fortran_COMPILER_SUPPORTS_ISO_C_BINDING)
    message(STATUS "${CMAKE_Fortran_COMPILER} supports iso_c_binding")
  else()  
    message(ERROR "Output from compile test=${out_result}")
    message(FATAL_ERROR "Fortran interfaces depend on iso_c_binding."
                        " ${CMAKE_Fortran_COMPILER} does not appear"
			" to support this module")
  endif()		      

endif()

# ---------------------------------------------------------------------------- #
# C Compiler Flags
# ---------------------------------------------------------------------------- #

# Add the PIC flag (Position In Code)
if (ENABLE_PIC)
  include(CheckCCompilerFlag)
  set(test_pic_flags -fPIC -fpic -PIC -pic)
  set(pic_flag pic_flag-NOTFOUND)
  set(i 1)
  foreach(flag ${test_pic_flags})
    set(res_var pic_flag_test${i})
    check_c_compiler_flag(${flag} ${res_var})
    if ( ${${res_var}} )
      set(pic_flag ${flag})
      break()
    endif(${${res_var}})  
    math(EXPR i "${i}+1")
  endforeach()

  if(pic_flag)
    message(STATUS "Adding ${pic_flag} to CMAKE_C_FLAGS") 
    set(CMAKE_C_FLAGS "${pic_flag} ${CMAKE_C_FLAGS}")
  else()
    message(FATAL_ERROR "Could not determine the PIC flag for this compiler. "
                        "Tried: ${test_pic_flags} all failed")
  endif()		      

endif()  

#DEBUGif(CMAKE_C_COMPILER_ID MATCHES GNU)
#DEBUG
#DEBUG  # GNU
#DEBUG
#DEBUG  set(C_COMPILER_PIC_OPT "-fPIC")
#DEBUG
#DEBUG  append_set(CMAKE_C_FLAGS "")
#DEBUG  append_set(CMAKE_C_FLAGS_DEBUG "-g -Wall")
#DEBUG  append_set(CMAKE_C_FLAGS_RELEASE "-O3")
#DEBUG
#DEBUGelseif(CMAKE_C_COMPILER_ID MATCHES PGI)
#DEBUG  
#DEBUG  set(C_COMPILER_PIC_OPT "-fpic")
#DEBUG
#DEBUG  # PGI Compiler flags here
#DEBUG  append_set(CMAKE_C_FLAGS_DEBUG "-g -Mbounds -Minfo=all -Minform=inform")
#DEBUG  append_set(CMAKE_C_FLAGS_RELEASE "-fast")
#DEBUG
#DEBUGelseif( CMAKE_C_COMPILER_ID MATCHES Intel )  
#DEBUG  
#DEBUG  set(C_COMPILER_PIC_OPT "-fPIC")
#DEBUG
#DEBUG  # Intel Compiler definitions
#DEBUG  append_set(CMAKE_C_FLAGS "")
#DEBUG  append_set(CMAKE_C_FLAGS_DEBUG "-g -Wall -Wcheck")
#DEBUG  append_set(CMAKE_C_FLAGS_RELEASE "-fast")
#DEBUG
#DEBUGelse()  
#DEBUG
#DEBUG  message(STATUS "Unknown C compiler type")
#DEBUG
#DEBUG  set(C_COMPILER_PIC_OPT "-fPIC")
#DEBUG
#DEBUG  append_set(CMAKE_C_FLAGS "")
#DEBUG  append_set(CMAKE_C_FLAGS_DEBUG   "-g")
#DEBUG  append_set(CMAKE_C_FLAGS_RELEASE "-O3")
#DEBUG
#DEBUGendif()  
#DEBUG
#DEBUG# Python interfaces require PIC (Position In Code) switches
#DEBUGif (ENABLE_Python)
#DEBUG  if(NOT BUILD_SHARED_LIBS)
#DEBUG    message(STATUS "Adding fPIC to C Flags")
#DEBUG    if( NOT ("${CMAKE_C_FLAGS}" MATCHES "${C_COMPILER_PIC_OPT}") )
#DEBUG      append_set(CMAKE_C_FLAGS "${C_OMPILER_PIC_OPT}")
#DEBUG    endif()
#DEBUG  endif()
#DEBUGendif()  
#DEBUG
#DEBUG# ---------------------------------------------------------------------------- #
#DEBUG# Fortran Compiler Flags
#DEBUG# ---------------------------------------------------------------------------- #
#DEBUGset(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
#DEBUG
#DEBUGif(ENABLE_Fortran)
#DEBUG  if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
#DEBUG  
#DEBUG    append_set(CMAKE_Fortran_FLAGS "")
#DEBUG    append_set(CMAKE_Fortran_FLAGS_DEBUG "-g")
#DEBUG    append_set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
#DEBUG  
#DEBUG  elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
#DEBUG  
#DEBUG    # PGI Compiler flags here
#DEBUG    append_set(CMAKE_Fortran_FLAGS_DEBUG "-g -Mbounds -Minfo=all -Minform=inform")
#DEBUG    append_set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
#DEBUG  
#DEBUG  elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)  
#DEBUG  
#DEBUG    #Intel Compiler flags here
#DEBUG
#DEBUG  else()  
#DEBUG  
#DEBUG    message(STATUS "Unknown Fortran compiler type")
#DEBUG    
#DEBUG    append_set(CMAKE_Fortran_FLAGS "")
#DEBUG    append_set(CMAKE_Fortran_FLAGS_DEBUG "-g")
#DEBUG    append_set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
#DEBUG  
#DEBUG  endif()
#DEBUG
#DEBUGendif()  
#DEBUG
#DEBUG# ---------------------------------------------------------------------------- #
#DEBUG# Throw compiler warnings here
#DEBUG# ---------------------------------------------------------------------------- #
#DEBUGif ( ENABLE_Fortran AND ( CMAKE_Fortran_COMPILER_ID MATCHES Intel ) )
#DEBUG  message(WARNING "Danu is NOT compatible with Intel version < 11.0")
#DEBUGendif()    
#DEBUG    
