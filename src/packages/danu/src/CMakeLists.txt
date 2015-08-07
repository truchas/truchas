# ############################################################################ #
#
#
#   DANU CMake 
#      Source File Directory
#
# ############################################################################ #

project(DanuLibSrc)

include(PrintVariable)
include(DanuAddLibrary)

# ---------------------------------------------------------------------------- #
# Add include directories
# ---------------------------------------------------------------------------- #
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${Danu_SOURCE_DIR}/src)

if ( ENABLE_Fortran )
  set(DANU_Fortran_MODULE_DIR ${Danu_BINARY_DIR}/src/mod_files)
  include_directories(${DANU_Fortran_MODULE_DIR})
  global_set(Danu_Fortran_MODULE_DIR ${DANU_Fortran_MODULE_DIR})
endif()


# ---------------------------------------------------------------------------- #
# Library: danu (C interface to HDF5)
# ---------------------------------------------------------------------------- #
set(DANU_SRC
        attribute.c
        dataset.c
        error.c
        file.c
        group.c
        h5_error.c
        h5_object.c
        danu_xdmf_mesh.c
        link.c
        linked_list.c
        memory.c
        slab.c
        utils.c
        mesh.c
        non-series_data.c
        probes.c
	offset.c
        output.c
        series_data.c
        series_group.c
        simulation.c)

set(DANU_HEADERS
         danu.h
         danu_attribute.h
         danu_chunk.h
         danu_dataset.h
         danu_errno.h
         danu_error.h
         danu_file.h
         danu_group.h
         danu_h5_error.h
         danu_h5_object.h
         danu_link.h
         danu_linked_list.h
         danu_memory.h
         danu_xdmf_mesh.h
         danu_slab.h
         danu_types.h
         danu_utils.h
         danu_mesh.h
         danu_non-series.h
         danu_probes.h
         danu_offset.h
         danu_output.h
         danu_series-data.h
         danu_series-group.h
         danu_sim.h)

danu_add_library(danu
                 SOURCE             ${DANU_SRC}
                 HEADERS            ${DANU_HEADERS} 
                 LINK_LIBS          ${HDF5_C_LIBRARY} ${HDF5_LINK_LIBRARIES}
		 DESTINATION        ${Danu_LIB_INSTALL_DIR}
		 HEADER_DESTINATION ${Danu_INCLUDE_INSTALL_DIR}
                 DEPENDENCIES       ${HDF5_BUILD_TARGET})

# ---------------------------------------------------------------------------- #
# Library: danu_f (Fortran interface)
# ---------------------------------------------------------------------------- #
if (ENABLE_Fortran)
  message(STATUS "Building Fortran interface")
  set(DANU_F_SRC 
            module_danu.f90
            module_danu_iface.f90)

  set(DANU_FC_SRC
            fort_error.c
	    fort_transpose.c
            fort_hid.c
            fort_strings.c
            fort_file.c
            fort_output.c
            fort_attribute.c
            fort_group.c
            fort_mesh.c
            fort_sim.c
            fort_non-series.c
            fort_series-group.c
            fort_series-data.c
            fort_probes.c)

  set(DANU_FC_HEADERS                 
            danu_fort_error.h
            danu_fort_transpose.h
            danu_fort_hid.h
            danu_fort_strings.h
            danu_fort_file.h
            danu_fort_attribute.h
            danu_fort_group.h
            danu_fort_mesh.h
            danu_fort_sim.h
            danu_fort_non-series.h
            danu_fort_series-group.h
            danu_fort_series-data.h
            danu_fort_probes.h)

 
  set(CMAKE_Fortran_MODULE_DIRECTORY ${DANU_Fortran_MODULE_DIR})


  danu_add_library(danu_f
                   SOURCE             ${DANU_F_SRC} ${DANU_FC_SRC}
                   HEADERS            ${DANU_FC_HEADERS}
		   DESTINATION        ${Danu_LIB_INSTALL_DIR}
		   HEADER_DESTINATION ${Danu_INCLUDE_INSTALL_DIR}
                   LINK_LIBS danu)

                 
  install(CODE
    "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_Fortran_MODULE_DIRECTORY}\" \"${Danu_MODULE_INSTALL_DIR}\")"
  )

endif()

# Add the Python interface
if ( ENABLE_Python )

  message(STATUS "Building the Python interface")
  add_subdirectory(python)

  # Push these definitions up one level
  set(PYTHON_FOUND ${PYTHON_FOUND} PARENT_SCOPE)
  set(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE} PARENT_SCOPE)
  set(PYTHON_VERSION ${PYTHON_VERSION} PARENT_SCOPE)
  set(PYTHON_LIBRARIES ${PYTHON_LIBRARIES} PARENT_SCOPE)
  set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} PARENT_SCOPE)
  set(NUMPY_VERSION ${NUMPY_VERSION} PARENT_SCOPE)
  set(NUMPY_INCLUDE_DIRS ${NUMPY_INCLUDE_DIRS} PARENT_SCOPE)


endif()
