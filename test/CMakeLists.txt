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

add_python_test(ds11
                ${TruchasTest_SOURCE_DIR}/ds11/testDS11.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(static_wall
                ${TruchasTest_SOURCE_DIR}/static-wall/testStaticWall.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 
