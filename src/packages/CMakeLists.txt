set(ENABLE_MPI True) # Some packages have both serial and parallel versions

message(STATUS "Configuring src/packages/chaco")
add_subdirectory(chaco)

message(STATUS "Configuring src/packages/chaparral")
add_subdirectory(chaparral)

message(STATUS "Configuring src/packages/pgslib")
add_subdirectory(pgslib)

message(STATUS "Configuring src/packages/ubiksolve")
add_subdirectory(ubiksolve)

message(STATUS "Configuring src/packages/danu")
add_subdirectory(danu)
