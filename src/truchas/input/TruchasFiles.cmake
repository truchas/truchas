# Truchas files in directory
#   input

set(INPUT_SOURCE_FILES
         input/EM_input.F90
         input/bc_input_module.F90
         input/body_input_module.F90
         input/input_driver.F90
         input/lin_solver_input.F90
         input/nonlin_solver_input.F90
         input/numerics_input_module.F90
         input/outputs_input_module.F90
         input/physics_input_module.F90
         input/region_data.F90
         input/region_input_module.F90
         input/solid_mechanics_namelist.F90
)

set(fc_flags
  -I${UbikSolve_MODULE_DIR} -I${PGSLib_MODULE_DIR} -I${Truchas_utilities_dir})
build_whitespace_string(INPUT_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${INPUT_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${INPUT_COMPILE_FLAGS})


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${INPUT_SOURCE_FILES})
