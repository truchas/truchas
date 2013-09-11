project(TruchasTest)

message(STATUS "configuring the Truchas test directory")

include(AddPythonTest)

add_python_test(RestartSuite
                ${TruchasTest_SOURCE_DIR}/restart/testRestartSuite.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(remap
                ${TruchasTest_SOURCE_DIR}/remap/testRemap.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(em1
                ${TruchasTest_SOURCE_DIR}/em1/testEM1.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds1
                ${TruchasTest_SOURCE_DIR}/ds1/testDS1.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds2
                ${TruchasTest_SOURCE_DIR}/ds2/testDS2.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds3
                ${TruchasTest_SOURCE_DIR}/ds3/testDS3.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds4
                ${TruchasTest_SOURCE_DIR}/ds4/testDS4.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds5
                ${TruchasTest_SOURCE_DIR}/ds5/testDS5.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds6
                ${TruchasTest_SOURCE_DIR}/ds6/testDS6.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds7
                ${TruchasTest_SOURCE_DIR}/ds7/testDS7.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds8
                ${TruchasTest_SOURCE_DIR}/ds8/testDS8.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

# Not working -- problems with coordinates and connectivity functions
#add_python_test(ds9
#                ${TruchasTest_SOURCE_DIR}/ds9/testDS9.py
#		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(htvoid1
                ${TruchasTest_SOURCE_DIR}/htvoid1/testHTVoid1.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(htvoid2
                ${TruchasTest_SOURCE_DIR}/htvoid2/testHTVoid2.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(htvoid3
                ${TruchasTest_SOURCE_DIR}/htvoid3/testHTVoid3.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(htvoid4
                ${TruchasTest_SOURCE_DIR}/htvoid4/testHTVoid4.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(broken-dam
                ${TruchasTest_SOURCE_DIR}/broken-dam/testBrokenDam.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(ds11
                ${TruchasTest_SOURCE_DIR}/ds11/testDS11.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(static_wall
                ${TruchasTest_SOURCE_DIR}/static-wall/testStaticWall.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(gap-rad
                ${TruchasTest_SOURCE_DIR}/gap-rad/testGapRad.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(vfrad1
                ${TruchasTest_SOURCE_DIR}/vfrad1/testVFRad1.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(vfrad2
                ${TruchasTest_SOURCE_DIR}/vfrad2/testVFRad2.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(sloshing-flow
                ${TruchasTest_SOURCE_DIR}/sloshing-flow/testSloshingFlow.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(contact-box-close
                ${TruchasTest_SOURCE_DIR}/contact-box-close/testContactBoxClose.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(contact-box-open
                ${TruchasTest_SOURCE_DIR}/contact-box-open/testContactBoxOpen.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

