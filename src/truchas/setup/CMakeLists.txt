target_sources(truchas PRIVATE setup_module.F90)

add_subdirectory(base_types)
add_subdirectory(initialize)
add_subdirectory(overwrite)
add_subdirectory(restart)
add_subdirectory(scalars)
add_subdirectory(vof)
