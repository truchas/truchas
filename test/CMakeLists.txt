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

add_python_test(ds11
                ${TruchasTest_SOURCE_DIR}/ds11/testDS11.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(static_wall
                ${TruchasTest_SOURCE_DIR}/static-wall/testStaticWall.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(gap_rad
                ${TruchasTest_SOURCE_DIR}/gap-rad/testGapRad.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 
