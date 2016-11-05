project(TruchasToolsRade C Fortran)

# --- CMake modules
include(FortranPreprocessFiles)

# --- Include paths

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Want modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})

# Chaparral
include_directories(${ChaparralVF_SOURCE_DIR})

# MPI
if ( ENABLE_MPI )
  include_directories(${MPI_C_INCLUDE_PATH})
endif()

# --- Compile defines
if(ENABLE_MPI)
  add_definitions(-DUSE_MPI)
else()
  add_definitions(-DVF_NO_MPI)
endif()

# --- Preprocessing flags
set(fpp_flags_list 
    ${Truchas_FPP_FLAGS} -I${TruchasToolsRade_SOURCE_DIR})
build_whitespace_string(fpp_flags ${fpp_flags_list})

# --- Library libre
set(LIBRE_FORTRAN_FILES
          scl.F90
	  f90_assert.F90
	  chaparral_c_binding.F90
	  re_utilities.F90
	  re_encl_type.F90
	  re_dist_vf_type.F90 
          re_exodus_encl.F90
	  re_chaparral_vf.F90
	  re_graphics_gmv.F90)

fortran_preprocess_files(LIBRARY_RADE_FILES
                         FILES ${LIBRE_FORTRAN_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${fpp_flags_list}
			 PROCESS_TARGET ProcessRadeFiles)
list(APPEND LIBRARY_RADE_FILES
            chaparral_ext.c
	    scl_f77.c)
                            
add_library(rade ${LIBRARY_RADE_FILES})
add_dependencies(rade ProcessRadeFiles vflib)
target_link_libraries(rade truchas ${Truchas_LINK_LIBS} vflib)

# --- cmpre binary
fortran_preprocess_files(CMPRE_SOURCE_FILES
                         FILES cmpre.F90
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${fpp_flags_list}
			 PROCESS_TARGET ProcessCmpreFiles)
add_executable(cmpre ${CMPRE_SOURCE_FILES})
add_dependencies(cmpre ProcessCmpreFiles)
target_link_libraries(cmpre rade truchas ${Truchas_LINK_LIBS})
if(ENABLE_MPI)
  set(output_name pcmpre)
else()
  set(output_name cmpre)
endif()
set_target_properties(cmpre PROPERTIES
                      OUTPUT_NAME ${output_name})

# --- genre binary
fortran_preprocess_files(GENRE_SOURCE_FILES
                         FILES genre_command_line.F90 genre.F90
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${fpp_flags_list}
			 PROCESS_TARGET ProcessGenreFiles)
add_executable(genre ${GENRE_SOURCE_FILES})
add_dependencies(genre ProcessGenreFiles)
target_link_libraries(genre rade truchas ${Truchas_LINK_LIBS})
if(ENABLE_MPI)
  set(output_name pgenre)
else()
  set(output_name genre)
endif()
set_target_properties(genre PROPERTIES
                      OUTPUT_NAME ${output_name})

# --- vizre binary
fortran_preprocess_files(VIZRE_SOURCE_FILES
                         FILES vizre_command_line.F90 vizre.F90
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${fpp_flags_list}
			 PROCESS_TARGET ProcessVizreFiles)
add_executable(vizre ${VIZRE_SOURCE_FILES})
add_dependencies(vizre ProcessVizreFiles)
target_link_libraries(vizre rade truchas ${Truchas_LINK_LIBS})
if(ENABLE_MPI)
  set(output_name pvizre)
else()
  set(output_name vizre)
endif()
set_target_properties(vizre PROPERTIES
                      OUTPUT_NAME ${output_name})
 
# --- Install rules 		    
install(TARGETS rade cmpre genre vizre
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
