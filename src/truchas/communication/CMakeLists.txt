target_sources(truchas PRIVATE
  parallel_communication.F90
  parallel_communication-broadcast_impl.F90
  parallel_communication-gather_impl.F90
  parallel_communication-scatter_impl.F90
  parallel_communication-reduce_impl.F90
  index_map_type.F90
  index_map_type-gather_offp_impl.F90
  index_map_type-scatter_offp_impl.F90
  index_map_type-gather_impl.F90
  index_map_type-scatter_impl.F90
  index_map_type-localize_impl.F90
)

# We are using MPI's Fortran "mpi" module and call the MPI functions
# with different data types. Some compilers will detect this and report
# compilation errors. Here we add compilation flags that tell the compiler
# to relax and not worry about the apparent contradictory usages.
# NB: This issue goes away if we were to use the "mpi_f08" module instead.

set(MPI_USERS
  parallel_communication.F90
  parallel_communication-broadcast_impl.F90
  parallel_communication-gather_impl.F90
  parallel_communication-scatter_impl.F90
  parallel_communication-reduce_impl.F90
  index_map_type.F90
  index_map_type-gather_offp_impl.F90
  index_map_type-scatter_offp_impl.F90
  index_map_type-gather_impl.F90
  index_map_type-scatter_impl.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  set_source_files_properties(${MPI_USERS}
    TARGET_DIRECTORY truchas PROPERTIES COMPILE_OPTIONS "-mismatch;-w")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10")
    set_source_files_properties(${MPI_USERS}
        TARGET_DIRECTORY truchas PROPERTIES COMPILE_OPTIONS "-fallow-argument-mismatch;-w")
  endif()
endif()

add_subdirectory(test)
