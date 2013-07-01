project(TBrookParserExtensions)

message(STATUS "configuring the TBrook parser extensions")

# For the C wrapper code
include_directories(${PYTHON_INCLUDE_DIRS})

# For the Fortran interface code
include_directories(${Truchas_MODULE_DIR})
include_directories(${NetCDF_INCLUDE_DIR})

# The asciiwriter extension
add_library(asciiwriter SHARED asciiwriter.c)
set_target_properties(asciiwriter PROPERTIES PREFIX "")

# The ensight writer extension
add_library(ensight SHARED ensight.c)
set_target_properties(ensight PROPERTIES PREFIX "")

# The grid-to-grid mapping extension
add_library(gridmap SHARED gridmap_c.c gridmap_f.f90)
target_link_libraries(gridmap truchas)
set_target_properties(gridmap PROPERTIES PREFIX "")

# The Exodus file I/O extension
add_library(exodus SHARED exodus_c.c exodus_f.f90)
target_link_libraries(exodus truchas)
set_target_properties(exodus PROPERTIES PREFIX "")

install(TARGETS asciiwriter ensight gridmap exodus
        LIBRARY DESTINATION ${Truchas_BIN_INSTALL_DIR})
