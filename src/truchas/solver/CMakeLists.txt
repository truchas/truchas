target_sources(truchas PRIVATE
  hypre_c_binding.F90
  fhypre.F90
  hypre_ext.c
  hypre_pcg_type.F90
  hypre_hybrid_type.F90
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
  bsr_matrix_type.F90
  cgsol_type.F90
  vector_class.F90
)
add_subdirectory(nka)

add_subdirectory(Test)
