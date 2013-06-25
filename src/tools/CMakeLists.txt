project(TruchasTools)

message(STATUS "configuring the Truchas tools directory")


# Post processing tools (probes, restart)
add_subdirectory(post_processors)

# Mesh processing
add_subdirectory(mesh_processors)

# Radiation enclosure 
add_subdirectory(RadE)

