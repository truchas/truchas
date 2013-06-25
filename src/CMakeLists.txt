project(TruchasSrc)

message(STATUS "configuring truchas source directory")

# Chaco
add_subdirectory(packages/chaco)

# Chaparral
add_subdirectory(packages/chaparral)

# PGSLib
if (ENABLE_PGSLib)
  add_subdirectory(packages/pgslib)
endif()

# UbikSolve
add_subdirectory(packages/ubiksolve)

# Danu
add_subdirectory(packages/danu)

# Truchas
add_subdirectory(truchas)

# Tools
add_subdirectory(tools)


