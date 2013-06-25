project(UbikSolveTest)

message(STATUS "configuring directory test")

# --- Remove the default Fortran module path
#     Do not want test modules in that directory
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# --- Include directories

# UbikSolve module directory
include_directories(${UbikSolve_MODULE_DIR})

# PGSLib include directories and module directory
if(ENABLE_PGSLib)
  include_directories(${PGSLib_MODULE_DIR})
  include_directories(${PGSLib_INCLUDE_DIRS})
endif()

# --- Define the source files
set(UbikTest_F_FILES
             src/UbikTest.F
             src/Driver.F
             src/ELL_data.F
             src/Full_data.F
             src/Iterate.F
             src/MatVec.F
             src/Precond.F
             src/SetUpSaad.F
             src/Usage.F
             src/test_data.F)

# Preprocess the files	    
fortran_preprocess_files(UbikTest_SOURCE_FILES
                         FILES ${UbikTest_F_FILES}
			 GREP_FILTERS ^$$
			 FPP_EXECUTABLE ${UbikSolve_PREPROCESSOR}
			 FPP_FLAGS -I${UbikSolveTest_SOURCE_DIR}/include
			 PROCESS_TARGET ProcessUbikTestFiles)

# --- Test harness executable		       
add_executable(UbikTest ${UbikTest_SOURCE_FILES})
add_dependencies(UbikTest ProcessUbikTestFiles)
set(UbikTest_LINK_LIBRARIES ubiksolve)
if (ENABLE_PGSLib)
  list(APPEND UbikTest_LINK_LIBRARIES ${PGSLib_LIBRARY})
endif() 
target_link_libraries(UbikTest ${UbikTest_LINK_LIBRARIES})

# --- Create an output directory for the test suite
set(UbikTest_OUTPUT_DIR ${UbikSolveTest_BINARY_DIR}/out)
add_custom_command(OUTPUT ${UbikTest_OUTPUT_DIR}
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${UbikTest_OUTPUT_DIR}
		   COMMENT "Creating test output directory"
		   VERBATIM)
add_custom_target(CreateUbikTestOutDir ALL DEPENDS ${UbikTest_OUTPUT_DIR})

# --- Add tests
#     UbikSolve had an extensive test suite what is duplicated here 
#     are the tests run from the runall_full target
#     

# - Problem 1 tests

# LU Tests
add_test(LU-1_1 UbikTest problem=1 assymm=0.2 method=LU scale=1)
add_test(LU-2_1 UbikTest problem=1 assymm=0.2 method=LU scale=2)

# BCGS
add_test(BCGS-none-0_1 UbikTest problem=1 assymm=0.2 method=BCGS scale=0)

# BCGSTAB
add_test(BCGSTAB-none-0_1 UbikTest problem=1 assymm=0.2 scale=0 method=BCGSTAB)
add_test(BCGSTAB-JAC-1.2_1 UbikTest problem=1 assymm=0.2 scale=1 method=BCGSTAB pre=JAC relax=1.2)
add_test(BCGSTAB-SSOR-1.2_1 UbikTest problem=1 assymm=0.2 scale=1 method=BCGSTAB pre=SSOR relax=1.2)
add_test(BCGSTAB-ILU-0.0_1 UbikTest problem=1 assymm=0.2 scale=1 method=BCGSTAB pre=ILU relax=0.0)
add_test(BCGSTAB-ILU-1.0_1 UbikTest problem=1 assymm=0.2 scale=1 method=BCGSTAB pre=ILU relax=1.0)
add_test(BCGSTAB-ILU-0.5_1 UbikTest problem=1 assymm=0.2 scale=1 method=BCGSTAB pre=ILU relax=0.5)

# TFQMR 
add_test(TFQMR-JAC-1.2_1 UbikTest problem=1 assymm=0.2 scale=1 method=TFQMR pre=JAC relax=1.2)
add_test(TFQMR-SSOR-1.2_1 UbikTest problem=1 assymm=0.2 scale=1 method=TFQMR pre=SSOR relax=1.2)
add_test(TFQMR-ILU-0.0_1 UbikTest problem=1 assymm=0.2 scale=1 method=TFQMR pre=ILU relax=0.0)
add_test(TFQMR-ILU-1.0_1 UbikTest problem=1 assymm=0.2 scale=1 method=TFQMR pre=ILU relax=1.0)
add_test(TFQMR-ILU-0.5_1 UbikTest problem=1 assymm=0.2 scale=1 method=TFQMR pre=ILU relax=0.5)

# GMRES 
add_test(GMRES-none-0_1 UbikTest problem=1 assymm=0.2 scale=0 method=GMRES)
add_test(GMRES-JAC-1.2_1 UbikTest problem=1 assymm=0.2 scale=1 method=GMRES pre=JAC relax=1.2)
add_test(GMRES-SSOR-1.2_1 UbikTest problem=1 assymm=0.2 scale=1 method=GMRES pre=SSOR relax=1.2)
add_test(GMRES-ILU-0.0_1 UbikTest problem=1 assymm=0.2 scale=1 method=GMRES pre=ILU relax=0.0)
add_test(GMRES-ILU-1.0_1 UbikTest problem=1 assymm=0.2 scale=1 method=GMRES pre=ILU relax=1.0)
add_test(GMRES-ILU-0.5_1 UbikTest problem=1 assymm=0.2 scale=1 method=GMRES pre=ILU relax=0.5)
add_test(GMRES-none-2_1 UbikTest problem=1 assymm=0.2 scale=2 method=GMRES)

# FGMRES 
add_test(FGMRES-none-0_1 UbikTest problem=1 assymm=0.2 scale=0 method=FGMRES)
add_test(FGMRES-JAC-1.2_1 UbikTest problem=1 assymm=0.2 scale=1 method=FGMRES pre=JAC relax=1.2)
add_test(FGMRES-SSOR-1.2_1 UbikTest problem=1 assymm=0.2 scale=1 method=FGMRES pre=SSOR relax=1.2)
add_test(FGMRES-ILU-0.0_1 UbikTest problem=1 assymm=0.2 scale=1 method=FGMRES pre=ILU relax=0.0)
add_test(FGMRES-ILU-1.0_1 UbikTest problem=1 assymm=0.2 scale=1 method=FGMRES pre=ILU relax=1.0)
add_test(FGMRES-ILU-0.5_1 UbikTest problem=1 assymm=0.2 scale=1 method=FGMRES pre=ILU relax=0.5)
add_test(FGMRES-none-2_1 UbikTest problem=1 assymm=0.2 scale=2 method=FGMRES)

# CG
add_test(CG-JAC-1.2_2 UbikTest problem=1 assymm=0 scale=1 method=CG pre=JAC relax=1.2)
add_test(CG-SSOR-1.2_2 UbikTest problem=1 assymm=0 scale=1 method=CG pre=SSOR relax=1.2)

# - Problem 2 tests

# LU Tests
add_test(LU-1_2 UbikTest problem=2 assymm=0.2 method=LU scale=1)
add_test(LU-2_2 UbikTest problem=2 assymm=0.2 method=LU scale=2)

# BCGS
add_test(BCGS-none-0_2 UbikTest problem=2 assymm=0.2 method=BCGS scale=0)

# BCGSTAB
add_test(BCGSTAB-none-0_2 UbikTest problem=2 assymm=0.2 scale=0 method=BCGSTAB)
add_test(BCGSTAB-JAC-1.2_2 UbikTest problem=2 assymm=0.2 scale=1 method=BCGSTAB pre=JAC relax=1.2)
add_test(BCGSTAB-SSOR-1.2_2 UbikTest problem=2 assymm=0.2 scale=1 method=BCGSTAB pre=SSOR relax=1.2)
add_test(BCGSTAB-ILU-0.0_2 UbikTest problem=2 assymm=0.2 scale=1 method=BCGSTAB pre=ILU relax=0.0)
add_test(BCGSTAB-ILU-1.0_2 UbikTest problem=2 assymm=0.2 scale=1 method=BCGSTAB pre=ILU relax=1.0)
add_test(BCGSTAB-ILU-0.5_2 UbikTest problem=2 assymm=0.2 scale=1 method=BCGSTAB pre=ILU relax=0.5)

# TFQMR 
add_test(TFQMR-JAC-1.2_2 UbikTest problem=2 assymm=0.2 scale=1 method=TFQMR pre=JAC relax=1.2)
add_test(TFQMR-SSOR-1.2_2 UbikTest problem=2 assymm=0.2 scale=1 method=TFQMR pre=SSOR relax=1.2)
add_test(TFQMR-ILU-0.0_2 UbikTest problem=2 assymm=0.2 scale=1 method=TFQMR pre=ILU relax=0.0)
add_test(TFQMR-ILU-1.0_2 UbikTest problem=2 assymm=0.2 scale=1 method=TFQMR pre=ILU relax=1.0)
add_test(TFQMR-ILU-0.5_2 UbikTest problem=2 assymm=0.2 scale=1 method=TFQMR pre=ILU relax=0.5)

# GMRES 
add_test(GMRES-none-0_2 UbikTest problem=2 assymm=0.2 scale=0 method=GMRES)
add_test(GMRES-JAC-1.2_2 UbikTest problem=2 assymm=0.2 scale=1 method=GMRES pre=JAC relax=1.2)
add_test(GMRES-SSOR-1.2_2 UbikTest problem=2 assymm=0.2 scale=1 method=GMRES pre=SSOR relax=1.2)
add_test(GMRES-ILU-0.0_2 UbikTest problem=2 assymm=0.2 scale=1 method=GMRES pre=ILU relax=0.0)
add_test(GMRES-ILU-1.0_2 UbikTest problem=2 assymm=0.2 scale=1 method=GMRES pre=ILU relax=1.0)
add_test(GMRES-ILU-0.5_2 UbikTest problem=2 assymm=0.2 scale=1 method=GMRES pre=ILU relax=0.5)
add_test(GMRES-none-2_2 UbikTest problem=2 assymm=0.2 scale=2 method=GMRES)

# FGMRES 
add_test(FGMRES-none-0_2 UbikTest problem=2 assymm=0.2 scale=0 method=FGMRES)
add_test(FGMRES-JAC-1.2_2 UbikTest problem=2 assymm=0.2 scale=1 method=FGMRES pre=JAC relax=1.2)
add_test(FGMRES-SSOR-1.2_2 UbikTest problem=2 assymm=0.2 scale=1 method=FGMRES pre=SSOR relax=1.2)
add_test(FGMRES-ILU-0.0_2 UbikTest problem=2 assymm=0.2 scale=1 method=FGMRES pre=ILU relax=0.0)
add_test(FGMRES-ILU-1.0_2 UbikTest problem=2 assymm=0.2 scale=1 method=FGMRES pre=ILU relax=1.0)
add_test(FGMRES-ILU-0.5_2 UbikTest problem=2 assymm=0.2 scale=1 method=FGMRES pre=ILU relax=0.5)
add_test(FGMRES-none-2_2 UbikTest problem=2 assymm=0.2 scale=2 method=FGMRES)

# CG
add_test(CG-JAC-1.2_2 UbikTest problem=2 assymm=0 scale=1 method=CG pre=JAC relax=1.2)
add_test(CG-SSOR-1.2_2 UbikTest problem=2 assymm=0 scale=1 method=CG pre=SSOR relax=1.2)
