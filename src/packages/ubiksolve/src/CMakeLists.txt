project(UbikSolveSource)

message(STATUS "Configuring directory src")

include(FortranPreprocessFiles)

# --- Include directories

# - Pick up the module directory for UbikSolve and PGSLib
include_directories(${UbikSolve_MODULE_DIR})
if (ENABLE_PGSLib)
  include_directories(${PGSLib_MODULE_DIR})
  include_directories(${PGSLib_INCLUDE_DIRS})
endif()

# Preprocess flags 
set(UbikSolve_FPP_FLAGS 
              -I${UbikSolveSource_SOURCE_DIR}
              -I${UbikSolve_SOURCE_DIR}/include)

# --- Define *.F -> *.f90 files
set(PGSLibIface_SOURCE_FILES)
set(UbikSolve_F_FILES
              BCGS.F
              BCGSTAB.F
              CG.F
              CSC.F
              Convergence.F
              ELL.F
              ErrorFunction.F
              FlushUnit.F
              Full.F
              GMRES.F
              Matrix.F
              Meta.F
              ParallelInfo.F
              ParseArg.F
              Stationary.F
              TFQMR.F
              Ubik.F
              UbikBase.F
	      UbikConstants.F
              UbikMatrix.F
              UbikSolve.F
              UbikUtil.F
              UbikVector.F
              Vector.F
              WriteString.F
              control_type.F
              numbers.F
              numbers_complete.F
              types.F
              vector_type.F)

# Preprocess the files	    
fortran_preprocess_files(UbikSolve_SOURCE_FILES
                         FILES ${UbikSolve_F_FILES}
			 GREP_FILTERS ^$$
			 FPP_EXECUTABLE ${UbikSolve_PREPROCESSOR}
			 FPP_FLAGS  ${UbikSolve_FPP_FLAGS}
			 PROCESS_TARGET ProcessUbikSrcFiles)

# --- Library

# Define the source files for the library
add_library(ubiksolve ${UbikSolve_SOURCE_FILES})

# Force the preprocessing to occur before building the library
add_dependencies(ubiksolve ProcessUbikSolveSrcFiles pgslib_fe)

# Add the PGSLib target if needed
if ( ENABLE_PGSLib )
  target_link_libraries(ubiksolve ${PGSLib_LIBRARY})
endif()

# Set the library properties
set_target_properties(ubiksolve PROPERTIES
                      OUTPUT_NAME       UbikSolve
	                  OUTPUT_NAME_DEBUG UbikSolve_g
		              LINK_LANGUAGE     Fortran)

# Define the installation location
install(TARGETS ubiksolve
        EXPORT  ubiksolve
        LIBRARY DESTINATION ${UbikSolve_LIBRARY_INSTALL_DIR}
        ARCHIVE DESTINATION ${UbikSolve_LIBRARY_INSTALL_DIR})
