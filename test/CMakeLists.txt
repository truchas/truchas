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

add_pytruchas_test(ds8-old ds8/test-old.py PROCESSORS 4 LABELS DIFFUSION OLD_FLOW)
#add_pytruchas_test(ds8 ds8/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)

add_pytruchas_test(ds9 ds9/testDS9.py
                   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(ds10 ds10/testDS10.py
                   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(ds11 ds11/testDS11.py
                   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(htvoid1-old htvoid1/test-old.py PROCESSORS 4 LABELS DIFFUSION OLD_FLOW)
add_pytruchas_test(htvoid2-old htvoid2/test-old.py PROCESSORS 4 LABELS DIFFUSION OLD_FLOW)
add_pytruchas_test(htvoid3-old htvoid3/test-old.py PROCESSORS 4 LABELS DIFFUSION OLD_FLOW)
add_pytruchas_test(htvoid4-old htvoid4/test-old.py PROCESSORS 4 LABELS DIFFUSION OLD_FLOW)

add_pytruchas_test(htvoid1 htvoid1/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)
add_pytruchas_test(htvoid2 htvoid2/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)
#add_pytruchas_test(htvoid3 htvoid3/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)
add_pytruchas_test(htvoid4 htvoid4/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)

add_pytruchas_test(broken-dam-old broken-dam-old/test.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(broken-dam broken-dam/test.py PROCESSORS 4 LABELS FLOW)

add_pytruchas_test(evap
                   evap/test.py
		   PROCESSORS 4
		   LABELS CONDUCTION)

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

add_pytruchas_test(nat-conv-old nat-conv-old/test1.py PROCESSORS 4 LABELS OLD_FLOW CONDUCTION HEX)
add_pytruchas_test(nat-conv-tet-old nat-conv-old/test2.py PROCESSORS 4 LABELS OLD_FLOW CONDUCTION TET)

add_pytruchas_test(nat-conv nat-conv/test.py PROCESSORS 4 LABELS FLOW CONDUCTION HEX)

add_pytruchas_test(porous-drag
                   porous-drag/testPorousDrag.py
		   PROCESSORS 4
		   LABELS OLD_FLOW)

add_pytruchas_test(sloshing-flow-old sloshing-flow-old/test.py PROCESSORS 4 LABELS OLD_FLOW)

#add_pytruchas_test(sloshing-flow sloshing-flow/test.py PROCESSORS 4 LABELS FLOW NOT_WORKING)

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

add_pytruchas_test(gap-rad-flow gap-rad-flow/test.py PROCESSORS 4 LABELS DIFFUSION FLOW)
add_pytruchas_test(gap-rad-flow-old gap-rad-flow/test-old.py PROCESSORS 4 LABELS DIFFUSION OLD_FLOW)

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
		   LABELS OLD_FLOW SOLID_MECHANICS)

add_pytruchas_test(tm-pc
                   tm-pc/testTMPC.py
		   PROCESSORS 4
		   LABELS SOLID_MECHANICS)

add_pytruchas_test(tangential-surface-tension-old
                   tangential-surface-tension/test_tst_old.py
		   PROCESSORS 4
		   LABELS DIFFUSION OLD_FLOW)

add_pytruchas_test(tangential-surface-tension
                   tangential-surface-tension/test_tst.py
		   PROCESSORS 4
		   LABELS DIFFUSION FLOW)

add_pytruchas_test(tangential-surface-tension-couette
                   tangential-surface-tension/test_tst_couette.py
		   PROCESSORS 4
		   LABELS DIFFUSION FLOW)

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

# FREEZING FLOW TESTS -- NEW FLOW SOLVER
add_pytruchas_test(freezing-flow-1 freezing-flow/test1.py PROCESSORS 4 LABELS DIFFUSION FLOW)
add_pytruchas_test(freezing-flow-2 freezing-flow/test2.py PROCESSORS 4 LABELS DIFFUSION FLOW)
add_pytruchas_test(freezing-flow-3 freezing-flow/test3.py PROCESSORS 4 LABELS DIFFUSION FLOW)

# FREE-SURFACE FLOW TESTS -- OLD FLOW SOLVER
add_pytruchas_test(free-surf-flow-old-1 free-surf-flow-old/test1.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(free-surf-flow-old-2 free-surf-flow-old/test2.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(free-surf-flow-old-3 free-surf-flow-old/test3.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(free-surf-flow-old-4 free-surf-flow-old/test4.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(free-surf-flow-old-5 free-surf-flow-old/test5.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(free-surf-flow-old-6 free-surf-flow-old/test6.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(free-surf-flow-old-7 free-surf-flow-old/test7.py PROCESSORS 4 LABELS OLD_FLOW NOT_WORKING)

# FREE-SURFACE FLOW TESTS
add_pytruchas_test(free-surf-flow-1 free-surf-flow/test1.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-2 free-surf-flow/test2.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-3 free-surf-flow/test3.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-4 free-surf-flow/test4.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-5 free-surf-flow/test5.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(free-surf-flow-6 free-surf-flow/test6.py PROCESSORS 4 LABELS FLOW)
#add_pytruchas_test(free-surf-flow-7 free-surf-flow/test7.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
add_pytruchas_test(free-surf-flow-8 free-surf-flow/test8.py PROCESSORS 4 LABELS FLOW)

# PURE ADVECTION TESTS -- OLD FLOW SOLVER
add_pytruchas_test(advection-old-0  advection-old/test0.py  PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(advection-old-1a advection-old/test1a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(advection-old-1b advection-old/test1b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(advection-old-1c advection-old/test1c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(advection-old-1d advection-old/test1d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(advection-old-2a advection-old/test2a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(advection-old-2b advection-old/test2b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(advection-old-2c advection-old/test2c.py PROCESSORS 4 LABELS OLD_FLOW)

# PURE ADVECTION TESTS -- NEW FLOW SOLVER
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


# INVISCID PIPE FLOW TESTS -- OLD FLOW SOLVER
# Commented-out tests will fail, but shouldn't, due to known solver issues
# see https://gitlab.com/truchas/truchas/issues/163
add_pytruchas_test(inviscid-pipe-flow-old-1a inviscid-pipe-flow-old/test1a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(inviscid-pipe-flow-old-1b inviscid-pipe-flow-old/test1b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(inviscid-pipe-flow-old-2a inviscid-pipe-flow-old/test2a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(inviscid-pipe-flow-old-2b inviscid-pipe-flow-old/test2b.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(inviscid-pipe-flow-old-1c inviscid-pipe-flow-old/test1c.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(inviscid-pipe-flow-old-1d inviscid-pipe-flow-old/test1d.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(inviscid-pipe-flow-old-2c inviscid-pipe-flow-old/test2c.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(inviscid-pipe-flow-old-2d inviscid-pipe-flow-old/test2d.py PROCESSORS 4 LABELS OLD_FLOW)

# PIPE FLOW TESTS -- OLD FLOW SOLVER
# Commented-out tests will fail, but shouldn't, due to known solver issues
# see https://gitlab.com/truchas/truchas/issues/163
add_pytruchas_test(pipe-flow-old-1a pipe-flow-old/test1a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(pipe-flow-old-1b pipe-flow-old/test1b.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-1c pipe-flow-old/test1c.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-1d pipe-flow-old/test1d.py PROCESSORS 4 LABELS OLD_FLOW)
# see https://gitlab.com/truchas/truchas/issues/212
#add_pytruchas_test(pipe-flow-old-2a pipe-flow-old/test2a.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-2b pipe-flow-old/test2b.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-2c pipe-flow-old/test2c.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(pipe-flow-old-2d pipe-flow-old/test2d.py PROCESSORS 4 LABELS OLD_FLOW)


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

# DIVERGING DUCT TESTS -- OLD FLOW SOLVER
add_pytruchas_test(diverging-duct-old-0 diverging-duct-old/test0.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(diverging-duct-old-1a diverging-duct-old/test1a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(diverging-duct-old-1b diverging-duct-old/test1b.py PROCESSORS 4 LABELS OLD_FLOW)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(diverging-duct-old-1c diverging-duct-old/test1c.py PROCESSORS 4 LABELS OLD_FLOW)
#add_pytruchas_test(diverging-duct-old-1d diverging-duct-old/test1d.py PROCESSORS 4 LABELS OLD_FLOW)

# DIVERGING DUCT TESTS -- NEW FLOW SOLVER
add_pytruchas_test(diverging-duct-0  diverging-duct/test0.py  PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(diverging-duct-1a diverging-duct/test1a.py PROCESSORS 4 LABELS FLOW)
add_pytruchas_test(diverging-duct-1b diverging-duct/test1b.py PROCESSORS 4 LABELS FLOW)
# These should work after the discretization is fixed
# see https://gitlab.com/truchas/truchas/issues/163
#add_pytruchas_test(diverging-duct-1c diverging-duct/test1c.py PROCESSORS 4 LABELS FLOW NOT_WORKING)
#add_pytruchas_test(diverging-duct-1d diverging-duct/test1d.py PROCESSORS 4 LABELS FLOW NOT_WORKING)

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

add_pytruchas_test(hydrostatic-old-10 hydrostatic-old/test10.py PROCESSORS 4 LABELS OLD_FLOW)

# HYDROSTATIC TESTS -- NEW FLOW SOLVER
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

add_pytruchas_test(steady-flow-old-1a steady-flow-old/test1a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-1b steady-flow-old/test1b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-1c steady-flow-old/test1c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-1d steady-flow-old/test1d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-2a steady-flow-old/test2a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-2b steady-flow-old/test2b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-2c steady-flow-old/test2c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-2d steady-flow-old/test2d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-3a steady-flow-old/test3a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-3b steady-flow-old/test3b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-3c steady-flow-old/test3c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-3d steady-flow-old/test3d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-4a steady-flow-old/test4a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-4b steady-flow-old/test4b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-4c steady-flow-old/test4c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-4d steady-flow-old/test4d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-5a steady-flow-old/test5a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-5b steady-flow-old/test5b.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-5c steady-flow-old/test5c.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-5d steady-flow-old/test5d.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-6a steady-flow-old/test6a.py PROCESSORS 4 LABELS OLD_FLOW)
add_pytruchas_test(steady-flow-old-6b steady-flow-old/test6b.py PROCESSORS 4 LABELS OLD_FLOW)


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
