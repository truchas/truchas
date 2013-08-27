project(PGSLibSerial)

message(STATUS "Configuring serial-src directory")

include(FortranPreprocessFiles)

include_directories(${PGSLib_MODULE_DIR})

# --- Define source files for the library
set(PGSLibSerial_F_FILES
                     gath-scatt/gs-f.F
		     indexing/index-f.F
		     io/input-f.F
		     reductions/redux-f.F
		     scans/scan-f.F
		     utility/utility-f.F)

fortran_preprocess_files(PGSLibSerial_SOURCE_FILES
                         FILES ${PGSLibSerial_F_FILES}
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibSerial_SOURCE_DIR}/scans
			 GREP_FILTERS !!CPP!!
			 PROCESS_TARGET ProcessPGSLibSerialFiles)

# --- Define the library target
add_library(pgslib_serial ${PGSLibSerial_SOURCE_FILES})
add_dependencies(pgslib_serial ProcessSerialFiles)

# Change name if Debug build
set_target_properties(pgslib_serial PROPERTIES
                      OUTPUT_NAME       pgslib-ser
		      OUTPUT_NAME_DEBUG pgslib-ser_g)

# Define the installation location
install(TARGETS pgslib_serial
        EXPORT  pgslib_serial
        LIBRARY DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR})


