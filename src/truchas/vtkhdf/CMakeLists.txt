target_sources(truchas PRIVATE
  hdf5_c_binding.F90
  hdf5_c_binding_ext.c
  hl_hdf5.F90
  vtkhdf_file_type.F90
  vtkhdf_ug_type.F90
)

set(MPI_USERS
  hl_hdf5.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  set_source_files_properties(${MPI_USERS}
    TARGET_DIRECTORY truchas PROPERTIES COMPILE_OPTIONS "-mismatch;-w")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set_source_files_properties(${MPI_USERS}
      TARGET_DIRECTORY truchas PROPERTIES COMPILE_OPTIONS "-fallow-argument-mismatch;-w")
endif()

add_subdirectory(test)
