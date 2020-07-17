set(ENABLE_MPI True) # Some packages have both serial and parallel versions

message(STATUS "Configuring src/packages/pgslib")
add_subdirectory(pgslib)

message(STATUS "Configuring src/packages/scorpio")
add_subdirectory(scorpio/src)

message(STATUS "Configuring src/packages/ubiksolve")
add_subdirectory(ubiksolve)
