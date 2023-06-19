project(TruchasTest NONE)

set(CMAKE_MODULE_PATH ${TruchasTest_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

include(AddPytruchasTest)

add_pytruchas_test(enthalpy1 enthalpy/test1.py PROCESSORS 1 LABELS CONDUCTION)
add_pytruchas_test(enthalpy2 enthalpy/test2.py PROCESSORS 1 LABELS CONDUCTION)
add_pytruchas_test(enthalpy3 enthalpy/test3.py PROCESSORS 1 LABELS CONDUCTION)
add_pytruchas_test(enthalpy4 enthalpy/test4.py PROCESSORS 1 LABELS CONDUCTION)
add_pytruchas_test(enthalpy5 enthalpy/test5.py PROCESSORS 1 LABELS CONDUCTION)

add_pytruchas_test(1d-ss-hc-xlinear-hex
                   1d-ss-hc-xlinear-hex/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION HEX)

add_pytruchas_test(1d-ss-hc-xlinear-tet
                   1d-ss-hc-xlinear-tet/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION TET)

add_pytruchas_test(1d-ss-hc-ylinear-hex
                   1d-ss-hc-ylinear-hex/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION HEX)

add_pytruchas_test(1d-ss-hc-ylinear-tet
                   1d-ss-hc-ylinear-tet/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION TET)

add_pytruchas_test(1d-ss-hc-zlinear-hex
                   1d-ss-hc-zlinear-hex/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION HEX)

add_pytruchas_test(1d-ss-hc-zlinear-tet
                   1d-ss-hc-zlinear-tet/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION TET)

add_pytruchas_test(3d-ss-hc-linear-hex
                   3d-ss-hc-linear-hex/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION HEX)

add_pytruchas_test(3d-ss-hc-linear-tet
                   3d-ss-hc-linear-tet/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION TET)

add_pytruchas_test(3d-ss-hc-linear-mixed
                   3d-ss-hc-linear-mixed/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION MIXED)

add_pytruchas_test(oriented-flux-1 oriented-flux/test1.py PROCESSORS 4 LABELS CONDUCTION)
add_pytruchas_test(oriented-flux-2 oriented-flux/test2.py PROCESSORS 4 LABELS CONDUCTION)
add_pytruchas_test(oriented-flux-3 oriented-flux/test3.py PROCESSORS 4 LABELS CONDUCTION)
add_pytruchas_test(oriented-flux-4 oriented-flux/Temp-depen-absorp/test.py PROCESSORS 4 LABELS CONDUCTION)


add_pytruchas_test(two-toolheads toolhead/two-toolheads/test.py PROCESSORS 4 LABELS CONDUCTION)
add_pytruchas_test(laser-power-function-time toolhead/laser-power-function-time/test.py PROCESSORS 4 LABELS CONDUCTION)

add_pytruchas_test(RestartSuite restart/testRestartSuite.py
                   PROCESSORS 4
                   LABELS PHASE_CHANGE)

add_pytruchas_test(remap restart/testRemap.py PROCESSORS 4)

#add_pytruchas_test(em1 em1/testEM1.py
#                   DEPENDS RestartSuite remap # Just a demo of how this works
#		   PROCESSORS 4
#		   LABELS DEMO ELECTROMAGNETIC)

add_pytruchas_test(em1 em1/testEM1.py
                   PROCESSORS 4
                   LABELS ELECTROMAGNETICS)

add_pytruchas_test(em2 em2/testEM2.py
                   PROCESSORS 4
                   LABELS ELECTROMAGNETICS)

add_pytruchas_test(em3 em3/testEM3.py
                   PROCESSORS 4
                   LABELS ELECTROMAGNETICS)

add_pytruchas_test(ds1 ds1/testDS1.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

add_pytruchas_test(ds2 ds2/testDS2.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

add_pytruchas_test(ds3 ds3/testDS3.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

add_pytruchas_test(ds4 ds4/testDS4.py
                   PROCESSORS 4
                   LABELS DIFFUSION PHASE_CHANGE)

add_pytruchas_test(ds5 ds5/testDS5.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

add_pytruchas_test(ds6 ds6/testDS6.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

add_pytruchas_test(ds7 ds7/testDS7.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

#add_pytruchas_test(ds8 ds8/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)

add_pytruchas_test(ds9 ds9/testDS9.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

add_pytruchas_test(ds10 ds10/testDS10.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

add_pytruchas_test(ds11 ds11/testDS11.py
                   PROCESSORS 4
                   LABELS DIFFUSION PHASE_CHANGE)

add_pytruchas_test(pbf-material pbf-material/test.py PROCESSORS 4 LABELS DIFFUSION PHASE_CHANGE)
add_pytruchas_test(diffusion-mtc diffusion/test-mtc.py PROCESSORS 4 LABELS DIFFUSION)

add_pytruchas_test(htvoid1 htvoid1/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)
add_pytruchas_test(htvoid2 htvoid2/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)
#add_pytruchas_test(htvoid3 htvoid3/test.py PROCESSORS 4 LABELS DIFFUSION PHASE_CHANGE FLOW)
add_pytruchas_test(htvoid4 htvoid4/test.py PROCESSORS 4 LABELS DIFFUSION RADIATION FLOW)

add_pytruchas_test(broken-dam broken-dam/test.py PROCESSORS 4 LABELS FLOW)

add_pytruchas_test(evap
                   evap/test.py
                   PROCESSORS 4
                   LABELS CONDUCTION)

add_pytruchas_test(nat-conv nat-conv/test.py PROCESSORS 4 LABELS FLOW CONDUCTION HEX)

add_pytruchas_test(sloshing-flow sloshing-flow/test.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(void void/test.py PROCESSORS 4 LABELS FLOW)

add_pytruchas_test(gap-rad
                   gap-rad/testGapRad.py
                   PROCESSORS 4
                   LABELS DIFFUSION)

add_pytruchas_test(gap-rad-flow gap-rad-flow/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)

add_pytruchas_test(vfrad1
                   vfrad1/testVFRad1.py
                   PROCESSORS 4
                   LABELS RADIATION CONDUCTION HEX)

add_pytruchas_test(vfrad1-tet
                   vfrad1-tet/testVFRad1Tet.py
                   PROCESSORS 4
                   LABELS RADIATION CONDUCTION TET)

add_pytruchas_test(vfrad1-patch
                   vfrad1-patch/testVFRad1Patch.py
                   PROCESSORS 4
           LABELS RADIATION CONDUCTION HEX)

add_pytruchas_test(vfrad2
                   vfrad2/testVFRad2.py
                   PROCESSORS 4
                   LABELS RADIATION CONDUCTION)

add_pytruchas_test(vfrad-moving vfrad-moving/test.py PROCESSORS 4 LABELS RADIATION CONDUCTION)
add_pytruchas_test(vfrad-heater vfrad-heater/test.py PROCESSORS 4 LABELS RADIATION CONDUCTION)

add_pytruchas_test(vfrad3
                   vfrad3/testVFRad3.py
                   PROCESSORS 4
                   LABELS RADIATION CONDUCTION)

add_pytruchas_test(shear
                   sm-shear/test-shear.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-x
                   sm-stretch/test-stretch-x.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-y
                   sm-stretch/test-stretch-y.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-z
                   sm-stretch/test-stretch-z.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-rotated
                   sm-stretch/test-stretch-rotated.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-hex
                   sm-stretch/test-stretch-hex.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-tet
                   sm-stretch/test-stretch-tet.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-wed
                   sm-stretch/test-stretch-wed.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-mixed
                   sm-stretch/test-stretch-mixed.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(stretch-contact
                   sm-stretch/test-stretch-contact.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(thermoelastic
                   sm-thermoelastic/test-thermoelastic.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(contact-0
                   sm-contact/test0.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(contact-1
                   sm-contact/test1.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(contact-2
                   sm-contact/test2.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(contact-box-close
                   contact-box-close/test.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(contact-box-open
                   contact-box-open/test.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(viscoplastic-ring
                   viscoplastic-ring/test.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(viscoplastic-1d
                   sm-viscoplastic/test.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(sm-restart
                   sm-viscoplastic/test-restart.py
                   PROCESSORS 4
                   LABELS SOLID_MECHANICS)

add_pytruchas_test(tangential-surface-tension
                   tangential-surface-tension/test_tst.py
                   PROCESSORS 4
                   LABELS DIFFUSION FLOW)

add_pytruchas_test(tangential-surface-tension-couette
                   tangential-surface-tension/test_tst_couette.py
                   PROCESSORS 4
                   LABELS DIFFUSION FLOW)

add_pytruchas_test(ustruc-ldrd
                   ustruc/test-ldrd.py
                   PROCESSORS 4
                   LABELS MICROSTRUCTURE CONDUCTION PHASE_CHANGE MIXED)

add_pytruchas_test(ustruc-ldrd-restart
                   ustruc/test-ldrd-restart.py
                   PROCESSORS 4
                   LABELS MICROSTRUCTURE CONDUCTION PHASE_CHANGE MIXED RESTART)

add_pytruchas_test(ustruc-gl-temp
                   ustruc/test-gl-temp.py
                   PROCESSORS 4
                   LABELS MICROSTRUCTURE CONDUCTION PHASE_CHANGE MIXED)

add_pytruchas_test(ustruc-gl-restart
                   ustruc/test-gl-restart.py
                   PROCESSORS 4
                   LABELS MICROSTRUCTURE CONDUCTION PHASE_CHANGE MIXED RESTART)

add_pytruchas_test(ustruc-gl-frac
                   ustruc/test-gl-frac.py
                   PROCESSORS 4
                   LABELS MICROSTRUCTURE CONDUCTION PHASE_CHANGE MIXED)

# FREEZING FLOW TESTS -- NEW FLOW SOLVER
add_pytruchas_test(freezing-flow-1 freezing-flow/test1.py PROCESSORS 4 LABELS DIFFUSION PHASE_CHANGE FLOW)
add_pytruchas_test(freezing-flow-2 freezing-flow/test2.py PROCESSORS 4 LABELS DIFFUSION PHASE_CHANGE FLOW)
add_pytruchas_test(freezing-flow-3 freezing-flow/test3.py PROCESSORS 4 LABELS DIFFUSION PHASE_CHANGE FLOW)

# FREE-SURFACE FLOW TESTS
add_pytruchas_test(free-surf-flow-1 free-surf-flow/test1.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-2 free-surf-flow/test2.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-3 free-surf-flow/test3.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-4 free-surf-flow/test4.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-5 free-surf-flow/test5.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-6 free-surf-flow/test6.py PROCESSORS 4 LABELS FLOW)
#add_pytruchas_test(free-surf-flow-7 free-surf-flow/test7.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(free-surf-flow-8 free-surf-flow/test8.py PROCESSORS 4 LABELS FLOW)

# PURE ADVECTION TESTS
add_pytruchas_test(advection-0  advection/test0.py  PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-1a advection/test1a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-1b advection/test1b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-1c advection/test1c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-1d advection/test1d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-2a advection/test2a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-2b advection/test2b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-2c advection/test2c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-3a advection/test3a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-3b advection/test3b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-3c advection/test3c.py PROCESSORS 4 LABELS FLOW)

add_pytruchas_test(advection-3a-exact advection/test3aex.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-3b-exact advection/test3bex.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(advection-3c-exact advection/test3cex.py PROCESSORS 4 LABELS FLOW)


# INVISCID PIPE FLOW TESTS
add_pytruchas_test(inviscid-pipe-flow-1a inviscid-pipe-flow/test1a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(inviscid-pipe-flow-1b inviscid-pipe-flow/test1b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(inviscid-pipe-flow-2a inviscid-pipe-flow/test2a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(inviscid-pipe-flow-2b inviscid-pipe-flow/test2b.py PROCESSORS 4 LABELS FLOW)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(inviscid-pipe-flow-1c inviscid-pipe-flow/test1c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(inviscid-pipe-flow-1d inviscid-pipe-flow/test1d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(inviscid-pipe-flow-2c inviscid-pipe-flow/test2c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(inviscid-pipe-flow-2d inviscid-pipe-flow/test2d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)

# PIPE FLOW TESTS
add_pytruchas_test(pipe-flow-1a pipe-flow/test1a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(pipe-flow-1b pipe-flow/test1b.py PROCESSORS 4 LABELS FLOW)
# see https://gitlab.com/truchas/truchas/issues/212
#add_pytruchas_test(pipe-flow-2a pipe-flow/test2a.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(pipe-flow-2b pipe-flow/test2b.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(pipe-flow-1c pipe-flow/test1c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(pipe-flow-1d pipe-flow/test1d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(pipe-flow-2c pipe-flow/test2c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(pipe-flow-2d pipe-flow/test2d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(pipe-flow-3a pipe-flow/test3a.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(pipe-flow-3b pipe-flow/test3b.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(pipe-flow-3c pipe-flow/test3c.py PROCESSORS 4 LABELS)
add_pytruchas_test(pipe-flow-3d pipe-flow/test3d.py PROCESSORS 4 LABELS)

# DIVERGING DUCT TESTS
add_pytruchas_test(diverging-duct-0  diverging-duct/test0.py  PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(diverging-duct-1a diverging-duct/test1a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(diverging-duct-1b diverging-duct/test1b.py PROCESSORS 4 LABELS FLOW)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(diverging-duct-1c diverging-duct/test1c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(diverging-duct-1d diverging-duct/test1d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)


# HYDROSTATIC TESTS
add_pytruchas_test(hydrostatic-1a hydrostatic/test1a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-1b hydrostatic/test1b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-1c hydrostatic/test1c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-1d hydrostatic/test1d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-2a hydrostatic/test2a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-2b hydrostatic/test2b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-2c hydrostatic/test2c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-2d hydrostatic/test2d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-3a hydrostatic/test3a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-3b hydrostatic/test3b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-3c hydrostatic/test3c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-3d hydrostatic/test3d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-4a hydrostatic/test4a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-4b hydrostatic/test4b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-4c hydrostatic/test4c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-4d hydrostatic/test4d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-5a hydrostatic/test5a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-5b hydrostatic/test5b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-5c hydrostatic/test5c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-5d hydrostatic/test5d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-6a hydrostatic/test6a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-6b hydrostatic/test6b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-6c hydrostatic/test6c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-6d hydrostatic/test6d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-7a hydrostatic/test7a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-7b hydrostatic/test7b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-8a hydrostatic/test8a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-8b hydrostatic/test8b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-8c hydrostatic/test8c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-8d hydrostatic/test8c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-9a hydrostatic/test9a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-9b hydrostatic/test9b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-9c hydrostatic/test9c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-9d hydrostatic/test9d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-10 hydrostatic/test10.py PROCESSORS 4 LABELS FLOW)


add_pytruchas_test(steady-flow-1a steady-flow/test1a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-1b steady-flow/test1b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-1c steady-flow/test1c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-1d steady-flow/test1d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-2a steady-flow/test2a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-2b steady-flow/test2b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-2c steady-flow/test2c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-2d steady-flow/test2d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-3a steady-flow/test3a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-3b steady-flow/test3b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-3c steady-flow/test3c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-3d steady-flow/test3d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-4a steady-flow/test4a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-4b steady-flow/test4b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-4c steady-flow/test4c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-4d steady-flow/test4d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-5a steady-flow/test5a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-5b steady-flow/test5b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-5c steady-flow/test5c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-5d steady-flow/test5d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-6a steady-flow/test6a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(steady-flow-6b steady-flow/test6b.py PROCESSORS 4 LABELS FLOW)

add_pytruchas_test(inflow-bc-1 inflow-bc/test1.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(inflow-bc-2 inflow-bc/test2.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(inflow-bc-3 inflow-bc/test3.py PROCESSORS 4 LABELS FLOW)

# POROUS DRAG TESTS
add_pytruchas_test(drag-1d drag/test-1d.py PROCESSORS 1 LABELS FLOW)
