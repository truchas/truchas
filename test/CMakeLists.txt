project(TruchasTest)

message(STATUS "configuring the Truchas test directory")

include(AddPythonTest)

add_python_test(restart 
                ${TruchasTest_SOURCE_DIR}/restart/restart.py
		PYTHON_PATH ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR}) 
