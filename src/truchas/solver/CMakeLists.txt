target_sources(truchas PRIVATE
  hypre_c_binding.F90
  fhypre.F90
  hypre_ext.c
  hypre_pcg_type.F90
  hypre_hybrid_type.F90
  pcsr_matrix_type.F90
  pcsr_precon_class.F90
  pcsr_precon_ssor_type.F90
  pcsr_precon_boomer_type.F90
  pcsr_precon_factory.F90
  ridders_class.F90
  nlsol_type.F90
)

add_subdirectory(Test)
