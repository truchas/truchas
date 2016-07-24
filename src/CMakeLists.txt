message(STATUS "Configuring src/packages/chaco")
add_subdirectory(packages/chaco)

message(STATUS "Configuring src/packages/chaparral")
add_subdirectory(packages/chaparral)

message(STATUS "Configuring src/packages/pgslib")
add_subdirectory(packages/pgslib)

message(STATUS "Configuring src/packages/ubiksolve")
add_subdirectory(packages/ubiksolve)

message(STATUS "Configuring src/packages/danu")
add_subdirectory(packages/danu)

message(STATUS "Configuring src/truchas")
add_subdirectory(truchas)

message(STATUS "Configuring src/tools")
add_subdirectory(tools)

message(STATUS "Configuring src/python")
add_subdirectory(python)
