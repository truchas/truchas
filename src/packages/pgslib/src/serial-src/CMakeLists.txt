project(PGSLibSerial Fortran)

include_directories(${PGSLib_MODULE_DIR})

# --- Define source files for the library
set(PGSLibSerial_SOURCE_FILES
                     gath-scatt/gs-f.f90
		     indexing/index-f.f90
		     io/input-f.f90
		     reductions/redux-f.f90
		     scans/scan-f.f90
		     utility/utility-f.f90)

# --- Define the library target
add_library(pgslib_serial ${PGSLibSerial_SOURCE_FILES})

# Change name if Debug build
set_target_properties(pgslib_serial PROPERTIES
                      OUTPUT_NAME       pgslib-ser
		      OUTPUT_NAME_DEBUG pgslib-ser_g)

# Define the installation location
install(TARGETS pgslib_serial
        EXPORT  pgslib_serial
        LIBRARY DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR})
