set(ENABLE_MPI True) # Some packages have both serial and parallel versions

message(STATUS "Configuring src/packages/pgslib")
add_subdirectory(pgslib)
