project(TruchasRegTests)

message(STATUS "configuring truchas regressiontests directory")

# Configure the Regression Test harness
set(TestHarness ${TruchasRegTests_BINARY_DIR}/runRegressionTests.py)
configure_file(${TruchasRegTests_SOURCE_DIR}/runRegressionTests.py
               ${TestHarness}
               @ONLY)

# The list of regression tests
set(Truchas_REGRESSION_TESTS
broken-dam
channel_flow
contact_box_close
contact_box_open
diverging_duct
ds1
ds10
ds11
ds2
ds3
ds4
ds5
ds6
ds7
ds8
ds9
em1
em2
em3
gap-rad
gap-rad-flow
htvoid1
htvoid2
htvoid3
htvoid4
hydrostatic
natural_convection
porous_drag
shear_constr_heat
simple_gap
sloshing_flow
static-wall
tm_pc_both_ds
turbulence
vfrad1
vfrad2
viscoplastic_ring
void_collapse
)
build_whitespace_string(Truchas_REGRESSION_TESTS_STRING 
                        ${Truchas_REGRESSION_TESTS})

# Need the name of the executable
get_target_property(EXECUTABLE_NAME truchas_exe OUTPUT_NAME)

# Need the location of the executable, will use the build directory
set(EXECUTABLE_DIR ${TruchasExe_BINARY_DIR})

# Test harness command that runs the list of tests
set(TestHarness_COMMAND
    ${PYTHON_EXECUTABLE} ${TestHarness}
    -b ${EXECUTABLE_DIR} -e ${EXECUTABLE_NAME}
    -o ${TruchasRegTests_BINARY_DIR}/regressiontests_output
)
if(ENABLE_MPI)
  list(APPEND TestHarness_COMMAND -n 4)
endif()

# Now create the custom target that runs the entire test suite
add_custom_target(sharen
                  DEPENDS truchas_exe
                  COMMAND ${TestHarness_COMMAND} -l ${Truchas_REGRESSION_TESTS_STRING}
                  VERBATIM
                  COMMENT "Running old regression test suite")

# Run the tests through the CTest Package
#foreach( test ${Truchas_REGRESSION_TESTS})
#  add_test(NAME ${test} COMMAND ${TestHarness_COMMAND} -l ${test})
#endforeach()
