project(TruchasTest)

message(STATUS "configuring the Truchas test directory")

set(CMAKE_MODULE_PATH ${TruchasTest_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

include(AddPytruchasTest)

add_pytruchas_test(RestartSuite restart/testRestartSuite.py 
                   PROCESSORS 4)

add_pytruchas_test(remap remap/testRemap.py PROCESSORS 4)

add_pytruchas_test(em1 em1/testEM1.py
                   DEPENDS RestartSuite remap # Just a demo of how this works
		   PROCESSORS 4
		   LABELS DEMO ELECTROMAGNETIC)


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

# Not working -- problems with coordinates and connectivity functions
add_pytruchas_test(ds9 ds9/testDS9.py
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
		   LABELS DIFFUSION)

add_pytruchas_test(broken-dam
                   broken-dam/testBrokenDam.py
		   PROCESSORS 4
		   LABELS FLOW)

add_pytruchas_test(diverging-duct
                   diverging-duct/testDivergingDuct.py
		   PROCESSORS 4
		   LABELS FLOW)

add_pytruchas_test(hydrostatic
                   hydrostatic/testHydrostatic.py
		   PROCESSORS 4
		   LABELS FLOW)

add_pytruchas_test(static_wall
                   static-wall/testStaticWall.py
		   PROCESSORS 4
		   LABELS FLOW)

add_pytruchas_test(gap-rad
                   gap-rad/testGapRad.py
		   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(gap-rad-flow
                   gap-rad-flow/testGapRadFlow.py
		   PROCESSORS 4
                   LABELS DIFFUSION FLOW)

add_pytruchas_test(vfrad1
                   vfrad1/testVFRad1.py
		   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(vfrad2
                   vfrad2/testVFRad2.py
		   PROCESSORS 4
		   LABELS DIFFUSION)

add_pytruchas_test(sloshing-flow
                   sloshing-flow/testSloshingFlow.py
		   PROCESSORS 4
		   LABELS FLOW)

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
		   LABELS FLOW)

add_pytruchas_test(tm-pc
                   tm-pc/testTMPC.py
		   PROCESSORS 4
		   LABELS SOLID_MECHANICS)

