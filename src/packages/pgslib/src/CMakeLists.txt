project(PGSLibSource)

# Serial library
add_subdirectory(serial-src)

# Parallel library
if (ENABLE_MPI)
  add_subdirectory(par-src)
endif()

# Module interface
add_subdirectory(interface)

# Test suite
if (ENABLE_TESTS)
  add_subdirectory(testing)
endif()

