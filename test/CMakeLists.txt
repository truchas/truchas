project(TruchasTest)

message(STATUS "configuring the Truchas test directory")

include(AddPythonTest)

add_python_test(restart 
                ${TruchasTest_SOURCE_DIR}/restart/restart.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 

add_python_test(remap
                ${TruchasTest_SOURCE_DIR}/remap/testRemap.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 
