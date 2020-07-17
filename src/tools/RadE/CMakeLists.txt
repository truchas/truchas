project(TruchasToolsRade C Fortran)

# --- Include paths

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Want modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})

# MPI
include_directories(${MPI_C_INCLUDE_PATH})

# Ensure necessary Intel compiler options are there.
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume realloc_lhs")
endif()

# --- Compile defines
add_definitions(-DUSE_MPI)

# --- Library libre
set(LIBRARY_RADE_FILES
    scl.F90
    f90_assert.F90
    chaparral_c_binding.F90
    re_utilities.F90
    re_encl_type.F90
    re_dist_vf_type.F90
    re_vf_diff_type.F90
    re_exodus_encl.F90
    re_toolpath.F90
    re_chaparral_vf.F90
    re_graphics_gmv.F90
    patching/re_patch_type.F90
    patching/patching_tools.F90
    patching/edge_hash_type.F90
    patching/edge_neighbor_table_type.F90
    patching/vsa_min_heap.F90
    patching/vsa_patch_type.F90
    patching/vsa_patching_type.F90
    patching/vac_min_heap.F90
    patching/vac_patching_type.F90)

list(APPEND LIBRARY_RADE_FILES
            chaparral_ext.c
	    scl_f77.c)

add_library(rade ${LIBRARY_RADE_FILES})
target_link_libraries(rade truchas ${Truchas_LINK_LIBS} chaparral)

# --- cmpre binary
set(CMPRE_SOURCE_FILES cmpre.F90)
add_executable(cmpre ${CMPRE_SOURCE_FILES})
target_link_libraries(cmpre rade truchas ${Truchas_LINK_LIBS})
set(output_name cmpre)
set_target_properties(cmpre PROPERTIES
                      OUTPUT_NAME ${output_name})

# --- genre binary
set(GENRE_SOURCE_FILES
  genre_command_line.F90
  genre.F90)
add_executable(genre ${GENRE_SOURCE_FILES})
target_link_libraries(genre rade truchas ${Truchas_LINK_LIBS})
set(output_name genre)
set_target_properties(genre PROPERTIES
                      OUTPUT_NAME ${output_name})

# --- vizre binary
set(VIZRE_SOURCE_FILES
  vizre_command_line.F90
  vizre.F90)
add_executable(vizre ${VIZRE_SOURCE_FILES})
target_link_libraries(vizre rade truchas ${Truchas_LINK_LIBS})
set(output_name vizre)
set_target_properties(vizre PROPERTIES
                      OUTPUT_NAME ${output_name})

# --- ctest binaries
add_subdirectory(Test)

# --- Install rules
install(TARGETS rade cmpre genre vizre
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
