project(TruchasRegTests)

message(STATUS "configuring truchas regressiontests directory")

# Configure the Regression Test harness
set(TestHarness ${TruchasRegTests_BINARY_DIR}/runRegressionTests.py)
configure_file(${TruchasRegTests_SOURCE_DIR}/runRegressionTests.py
               ${TestHarness}
               @ONLY)

# The list of regression tests
set(Truchas_REGRESSION_TESTS
broken_dam
broken-dam
channel_flow
chemistry
composite_conduction
conduction_45_so
conduction_90
contact_box_close
contact_box_open
diagonal_advection
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
enclosure_radiation
gap-rad
gap-rad-flow
grain_growth
grain_growth_rnt
ht_bc
htvoid1
htvoid2
htvoid3
htvoid4
hydrostatic
natural_convection
phasechange_CnK
phasechange_eutectic
phasechange_mixed
porous_drag
postprocessor
purematerial_solidification
sensitivity_stefan
shear_constr_heat
simple_gap
sloshing_flow
species_advection
static_drop
static-wall
t_dep_temp_bc
tet_shear_heat
thermocapillary_flow
tm_pc_both
tm_pc_both_ds
turbulence
vfrad1
vfrad2
viscoplastic_ring
void_collapse
)
build_whitespace_string(Truchas_REGRESSION_TESTS_STRING 
                        ${Truchas_REGRESSION_TESTS})
print_variable(Truchas_REGRESSION_TESTS_STRING)

# Need the name of the executable
get_target_property(EXECUTABLE_NAME truchas_exe OUTPUT_NAME)

# Need the location of the executable, will use the build directory
set(EXECUTABLE_DIR ${TruchasExe_BINARY_DIR})

# Now create the custom target
add_custom_target(sharen
                  DEPENDS truchas_exe
                  COMMAND ${PYTHON_EXECUTABLE} ${TestHarness} -b ${EXECUTABLE_DIR} -e ${EXECUTABLE_NAME} -o ${TruchasRegTests_BINARY_DIR}/regressiontests_output -l ${Truchas_REGRESSION_TESTS_STRING}
                  VERBATIM
                  COMMENT "Running old regression test suite")


