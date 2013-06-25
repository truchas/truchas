project(ChaparralVF)

# --- Include Directories
if (ENABLE_MPI)
  include_directories(${MPI_C_INCLUDE_PATH})
endif()
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

# --- Directory defines
if ( ENABLE_MPI )
  add_definitions(-DVF_MPI)
else()  
  add_definitions(-DVF_NO_MPI)
endif()

if ( NOT ENABLE_ExodusII )
  add_definitions(-DVF_NO_EXODUS_OUTPUT)
endif()

# Check the math functions
include(CheckFunctionExists)
check_function_exists(log log_exists)
if ( NOT log_exists )
  if (UNIX)
    find_library(MATH_LIBRARY m)
  else()
    set(MATH_LIBRARY)
  endif()  
endif()                 


# --- Source files
set(VF_SOURCE_FILES
        Claussenp.c 
        VF_BackFaceCull.c 
        VF_BoundingBox.c 
        VF_BSPtree.c 
        VF_CalcHemicube.c 
        VF_CalcPairwise.c 
        VF_CalcRead.c 
        VF_CalcVF_Analytic.c 
        VF_CalcVF_ComputePair.c 
        VF_CalcVF_Contour.c 
        VF_CalcVF_DoubleArea.c 
        VF_CalcVF_Gauss.c 
        VF_CalcVF_Hemicube.c 
        VF_CalcVF_HemicubeRow.c 
        VF_CalcVF_Hottel.c 
        VF_CalcVF_MonteCarlo.c 
        VF_CalcVF_Occluded.c 
        VF_CalcVF_Unoccluded.c 
        VF_ComputeFluxes.c 
        VF_DefineEnclosure.c 
        VF_DefineTopology.c 
        VF_EnclosureUtils.c 
        VF_FacetAndPolyUtils.c 
        VF_FindCandidates.c 
        VF_FindMinSeperationDist.c 
        VF_GetMatrix.c 
        VF_HemicubeProjectRow.c 
        VF_HemicubeSub.c 
        VF_HemicubeUtils.c 
        VF_Initialize.c 
        VF_MakeMatrixSymmetric.c 
        VF_MatrixCalcHemicube.c 
        VF_MatrixCalcPairwise.c 
        VF_MatrixIO.c 
        VF_MatrixUtils.c 
        VF_MemoryUtils.c 
        VF_OutputBanner.c 
        VF_ParallelUtils.c 
        VF_PartialEnclosure.c 
        VF_PolyClip.c 
        VF_PolyScan.c 
        VF_ProjectOntoHemicube.c 
        VF_RadSolveAZTEC.c 
        VF_RadSolve.c 
        VF_RadSolveCG.c 
        VF_RadSolveGMRES.c 
        VF_SamplePoly.c 
        VF_SetView.c 
        VF_SetViewPort.c 
        VF_ShaftCull.c 
        VF_SmoothMatrix.c 
        VF_TopologyConcat.c 
        VF_TopologyUtils.c 
        VF_Utilities.c 
        VF_Visibility.c 
        VF_WriteExodus.c 
        VF_WriteGenesis.c 
        VF_WriteINP.c)

# --- Library

# Define the library
add_library(vflib ${VF_SOURCE_FILES})

# Define the required link libraries
set(VF_LINK_LIBRARIES ${MATH_LIBRARY})
if (ENABLE_MPI)
  list(APPEND VF_LINK_LIBRARIES ${MPI_C_LIBRARIES})
endif()
if (ENABLE_ExodusII)
  list(APPEND VF_LINK_LIBRARIES ${EXODUSII_LIBRARIES})
endif()
if (VF_LINK_LIBRARIES)
  target_link_libraries(vflib ${VF_LINK_LIBRARIES})
endif()

# Target properties
set_target_properties(vflib PROPERTIES
                      OUTPUT_NAME VF)

# Define the installation location
install(TARGETS vflib
        EXPORT vflib
	LIBRARY DESTINATION ${Chaparral_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Chaparral_LIBRARY_INSTALL_DIR})





