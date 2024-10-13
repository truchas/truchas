target_sources(truchas PRIVATE
  hypre_c_binding.F90
  fhypre.F90
  hypre_ext.c
  hypre_utilities.F90
  hypre_pcg_type.F90
  hypre_hybrid_type.F90
  hypre_ams_type.F90
  csr_graph_type.F90
  pcsr_matrix_type.F90
  pcsr_precon_class.F90
  pcsr_precon_ssor_type.F90
  pcsr_precon_boomer_type.F90
  pcsr_precon_factory.F90
  pbsr_matrix_type.F90
  cg_solver_class.F90
  #fgmres_solver_type.F90 UNUSED + NEEDS UPDATING TO COMPILE
  gmres_left_solver_class.F90
  ridders_class.F90
  nlsol_type.F90
  nlk_solver_class.F90
  upper_packed_matrix_procs.F90
  csr_matrix_type.F90
  csr_precon_boomer_type.F90
  bsr_matrix_type.F90
  #cgsol_type.F90 UNUSED + NEEDS UPDATING TO COMPILE
  vector_class.F90
  complex_lin_op_class.F90
  cs_minres_solver_type.F90
  complex_pcsr_matrix_type.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  set_source_files_properties(cs_minres_solver_type.F90
    TARGET_DIRECTORY truchas PROPERTIES COMPILE_OPTIONS "-mismatch")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10")
    set_source_files_properties(cs_minres_solver_type.F90
        TARGET_DIRECTORY truchas PROPERTIES COMPILE_OPTIONS "-fallow-argument-mismatch")
  endif()
endif()

add_subdirectory(nka)

add_subdirectory(Test)
