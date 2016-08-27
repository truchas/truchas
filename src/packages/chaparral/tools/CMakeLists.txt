project(ChaparralTools C)

# --- Include directories
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)
include_directories(${ChaparralVF_SOURCE_DIR})
if (ENABLE_MPI)
  include_directories(${MPI_C_INCLUDE_PATH})
endif()

# --- Directory defines
if (ENABLE_MPI)
  add_definitions(-DVF_MPI)
else()
  add_definitions(-DVF_NO_MPI)
endif()

# exo2inp 
if (ENABLE_ExodusII)

  find_package(NetCDF)
  if ( NETCDF_FOUND )
    add_executable(exo2inp exo2inp.c)
    target_link_libraries(exo2inp vflib) 
    set(exo2inp_link_libs vflib ${EXODUSII_LIBRARIES} ${NETCDF_LIBRARIES})
  endif()
  
endif()

# SS2Mesh
if ( ENABLE_ExodusII)

  add_executable(SS2Mesh SS2Mesh.c)
  target_link_libraries(SS2Mesh vflib ${EXODUSII_LIBRARIES})

endif() 

# VFdriver
if (ENABLE_ExodusII)

  add_executable(VFdriver VFdriver.c)
  target_link_libraries(VFdriver vflib ${EXODUSII_LIBRARIES} ${MPI_C_LIBRARIES})

endif()

# sample - BROKEN missing vf_types.h
#add_executable(sample sample.c)
#target_link_libraries(sample vflib)

# testFF - BROKEN missing VF_Initialize symbol in libVF
#add_executable(testFF testFF.c)
#target_link_libraries(testFF vflib)
#add_test(testFF testFF)





