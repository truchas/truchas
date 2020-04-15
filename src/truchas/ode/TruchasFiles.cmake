# Truchas files in directory
#   ode

set(ODE_SOURCE_FILES
       ode/bdf2/bdf2_controller.F90
       ode/bdf2/bdf2_dae.F90
       ode/bdf2/bdf2_integrator.F90
       ode/bdf2/bdf2_kinds.F90
       ode/bdf2/bdf2_profiling.F90
       ode/bdf2/solution_history.F90)

set_source_files_properties(${ODE_SOURCE_FILES} PROPERTIES
  COMPILE_FLAGS
  "-I${PGSLib_MODULE_DIR} -I${Truchas_utilities_dir}")

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${ODE_SOURCE_FILES})
