# Truchas files in directory
#   solver

set(SOLVER_SOURCE_FILES
          solver/hypre_c_binding.F90
          solver/fhypre.F90
          solver/hypre_pcg_type.F90
          solver/hypre_hybrid_type.F90
          solver/pcsr_matrix_type.F90
          solver/pcsr_precon_class.F90
          solver/pcsr_precon_ssor_type.F90
          solver/pcsr_precon_boomer_type.F90
          solver/pcsr_precon_factory.F90
          solver/ridders_class.F90
          solver/nlsol_type.F90
          )

set(fc_flags
  -I${PGSLib_MODULE_DIR}
  -I${Truchas_utilities_dir})

build_whitespace_string(SOLVER_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${SOLVER_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${SOLVER_COMPILE_FLAGS})

list(APPEND SOLVER_SOURCE_FILES solver/hypre_ext.c)

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${SOLVER_SOURCE_FILES})
