project(TruchasTest NONE)

set(CMAKE_MODULE_PATH ${TruchasTest_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

include(AddPytruchasTest)

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

add_pytruchas_test(RestartSuite restart/testRestartSuite.py 
                   PROCESSORS 4)

add_pytruchas_test(remap remap/testRemap.py PROCESSORS 4)

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
		   LABELS DIFFUSION)

add_pytruchas_test(ds5 ds5/testDS5.py
		   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(ds6 ds6/testDS6.py
                   PROCESSORS 4 
		   LABELS DIFFUSION)

add_pytruchas_test(ds7 ds7/testDS7.py
                   PROCESSORS 4 
		   LABELS DIFFUSION)

add_pytruchas_test(ds8 ds8/testDS8.py
                   PROCESSORS 4 
		   LABELS DIFFUSION)

add_pytruchas_test(ds9 ds9/testDS9.py
                   PROCESSORS 4 
		   LABELS DIFFUSION)

add_pytruchas_test(ds10 ds10/testDS10.py
                   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(ds11 ds11/testDS11.py
                   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(htvoid1
                   htvoid1/testHTVoid1.py
		   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(htvoid2
                   htvoid2/testHTVoid2.py
		   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(htvoid3
                   htvoid3/testHTVoid3.py
		   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(htvoid4
                   htvoid4/testHTVoid4.py
		   PROCESSORS 4
		   LABELS DIFFUSION RADIATION)

add_pytruchas_test(broken-dam
                   broken-dam/testBrokenDam.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(channel-flow
                   channel-flow/testChannelFlow.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(diverging-duct-td
                   diverging-duct-td/testDivergingDuctTD.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(static-wall
                   static-wall/testStaticWall.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(natural-convection
                   natural-convection/testNaturalConvection.py
		   PROCESSORS 4
		   LABELS OLD_FLOW CONDUCTION HEX)

add_pytruchas_test(natural-conv-tet
                   natural-conv-tet/test.py
		   PROCESSORS 4
		   LABELS OLD_FLOW CONDUCTION TET)

add_pytruchas_test(porous-drag
                   porous-drag/testPorousDrag.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(sloshing-flow
                   sloshing-flow/testSloshingFlow.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(void-collapse
                   void-collapse/testVoidCollapse.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(turbulence
                   turbulence/testTurbulence.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(gap-rad
                   gap-rad/testGapRad.py
		   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(gap-rad-flow
                   gap-rad-flow/testGapRadFlow.py
		   PROCESSORS 4
                   LABELS DIFFUSION OLD_FLOW)

add_pytruchas_test(vfrad1
                   vfrad1/testVFRad1.py
		   PROCESSORS 4
		   LABELS RADIATION CONDUCTION HEX)

add_pytruchas_test(vfrad1-tet
                   vfrad1-tet/testVFRad1Tet.py
		   PROCESSORS 4
		   LABELS RADIATION CONDUCTION TET)

add_pytruchas_test(vfrad2
                   vfrad2/testVFRad2.py
		   PROCESSORS 4
		   LABELS RADIATION CONDUCTION)

add_pytruchas_test(vfrad3
                   vfrad3/testVFRad3.py
		   PROCESSORS 4
		   LABELS RADIATION CONDUCTION)

add_pytruchas_test(contact-box-close
                   contact-box-close/testContactBoxClose.py
		   PROCESSORS 4
		   LABELS SOLID_MECHANICS)

add_pytruchas_test(contact-box-open
                   contact-box-open/testContactBoxOpen.py
		   PROCESSORS 4
		   LABELS SOLID_MECHANICS)

add_pytruchas_test(viscoplastic-ring
                   viscoplastic-ring/testViscoplasticRing.py
		   PROCESSORS 4
		   LABELS SOLID_MECHANICS)

add_pytruchas_test(shear-constr-heat
                   shear-constr-heat/testShearConstrHeat.py
		   PROCESSORS 4
		   LABELS SOLID_MECHANICS)

add_pytruchas_test(simple-gap
                   simple-gap/testSimpleGap.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(tm-pc
                   tm-pc/testTMPC.py
		   PROCESSORS 4
		   LABELS SOLID_MECHANICS)

add_pytruchas_test(tangential-surface-tension
                   tangential-surface-tension/test_tst.py
		   PROCESSORS 4
		   LABELS DIFFUSION OLD_FLOW)

add_pytruchas_test(ustruc1
                   ustruc1/test.py
		   PROCESSORS 4
		   LABELS MICROSTRUCTURE CONDUCTION MIXED)

add_pytruchas_test(ustruc2
                   ustruc2/test.py
		   PROCESSORS 4
		   LABELS MICROSTRUCTURE CONDUCTION MIXED RESTART)

add_pytruchas_test(diag-advect
                   diag-advect/test.py
		   PROCESSORS 4
		   LABELS ADVECTION)

# PIPE FLOW TESTS -- OLD FLOW SOLVER
# Commented-out tests will fail, but shouldn't, due to known solver issues
# see https://gitlab.com/truchas/truchas/issues/163
add_pytruchas_test(pipe-flow-old-inviscid-1 pipe-flow-old/test1.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(pipe-flow-old-inviscid-2 pipe-flow-old/test2.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-inviscid-3 pipe-flow-old/test3.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-inviscid-4 pipe-flow-old/test4.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(pipe-flow-old-1 pipe-flow-old/test5.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(pipe-flow-old-2 pipe-flow-old/test6.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-3 pipe-flow-old/test7.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-4 pipe-flow-old/test8.py PROCESSORS 4 LABELS OLD_FLOW)

# PIPE FLOW TESTS -- NEW FLOW SOLVER
add_pytruchas_test(pipe-flow-inviscid-1 pipe-flow/test1.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(pipe-flow-inviscid-2 pipe-flow/test2.py PROCESSORS 4 LABELS FLOW)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(pipe-flow-inviscid-3 pipe-flow/test3.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(pipe-flow-inviscid-4 pipe-flow/test4.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(pipe-flow-1 pipe-flow/test5.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(pipe-flow-2 pipe-flow/test6.py PROCESSORS 4 LABELS FLOW)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(pipe-flow-3 pipe-flow/test7.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(pipe-flow-4 pipe-flow/test8.py PROCESSORS 4 LABELS FLOW NOT_WORKING)

# DIVERGING DUCT TESTS -- OLD FLOW SOLVER
add_pytruchas_test(diverging-duct-old-0 diverging-duct-old-0/test0.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(diverging-duct-old-1 diverging-duct-old-1/test1.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(diverging-duct-old-2 diverging-duct-old-2/test2.py PROCESSORS 4 LABELS OLD_FLOW)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(diverging-duct-old-3 diverging-duct-old-3/test3.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(diverging-duct-old-4 diverging-duct-old-4/test4.py PROCESSORS 4 LABELS OLD_FLOW)

# DIVERGING DUCT TESTS -- NEW FLOW SOLVER
add_pytruchas_test(diverging-duct-0 diverging-duct-0/test0.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(diverging-duct-1 diverging-duct-1/test1.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(diverging-duct-2 diverging-duct-2/test2.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(diverging-duct-3 diverging-duct-3/test3.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(diverging-duct-4 diverging-duct-4/test4.py PROCESSORS 4 LABELS FLOW NOT_WORKING)

# HYDROSTATIC TESTS -- OLD FLOW SOLVER
add_pytruchas_test(hydrostatic hydrostatic-old/testHydrostatic.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-1a hydrostatic-old/test1a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-1b hydrostatic-old/test1b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-1c hydrostatic-old/test1c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-1d hydrostatic-old/test1d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-2a hydrostatic-old/test2a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-2b hydrostatic-old/test2b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-2c hydrostatic-old/test2c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-2d hydrostatic-old/test2d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-3a hydrostatic-old/test3a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-3b hydrostatic-old/test3b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-3c hydrostatic-old/test3c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-3d hydrostatic-old/test3d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-4a hydrostatic-old/test4a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-4b hydrostatic-old/test4b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-4c hydrostatic-old/test4c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-4d hydrostatic-old/test4d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-5a hydrostatic-old/test5a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-5b hydrostatic-old/test5b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-5c hydrostatic-old/test5c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-5d hydrostatic-old/test5d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-6a hydrostatic-old/test6a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-6b hydrostatic-old/test6b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-6c hydrostatic-old/test6c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-6d hydrostatic-old/test6d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-7a hydrostatic-old/test7a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(hydrostatic-old-7b hydrostatic-old/test7b.py PROCESSORS 4 LABELS OLD_FLOW)

# HYDROSTATIC TESTS -- NEW FLOW SOLVER
add_pytruchas_test(hydrostatic-1a hydrostatic/test1a.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-1b hydrostatic/test1b.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-1c hydrostatic/test1c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-1d hydrostatic/test1d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-2a hydrostatic/test2a.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-2b hydrostatic/test2b.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-2c hydrostatic/test2c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-2d hydrostatic/test2d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-3a hydrostatic/test3a.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-3b hydrostatic/test3b.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-3c hydrostatic/test3c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-3d hydrostatic/test3d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-4a hydrostatic/test4a.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-4b hydrostatic/test4b.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-4c hydrostatic/test4c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-4d hydrostatic/test4d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-5a hydrostatic/test5a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-5b hydrostatic/test5b.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-5c hydrostatic/test5c.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-5d hydrostatic/test5d.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-6a hydrostatic/test6a.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-6b hydrostatic/test6b.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-6c hydrostatic/test6c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-6d hydrostatic/test6d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(hydrostatic-7a hydrostatic/test7a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(hydrostatic-7b hydrostatic/test7b.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
