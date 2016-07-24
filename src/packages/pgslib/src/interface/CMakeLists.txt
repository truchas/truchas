project(PGSLibIface C Fortran)

include(FortranPreprocessFiles)

# --- Include directories

# - Pick up the module directory
include_directories(${PGSLib_MODULE_DIR})

# --- Define *.F -> *.f90 files
set(PGSLibIface_SOURCE_FILES pgslib_c_binding.f90)

# gath-scatt files
set(GS_F_FILES
           gath-scatt/pgslib_gs_comm_module.F
	   gath-scatt/pgslib_gs_setup_module.F
           gath-scatt/pgslib_gs_module.F
	   gath-scatt/pgslib_gs_util_module.F)
fortran_preprocess_files(GS_SOURCE_FILES
                         FILES ${GS_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS  -I${PGSLibIface_SOURCE_DIR}/gath-scatt})
list(APPEND PGSLibIface_SOURCE_FILES ${GS_SOURCE_FILES})

# global files
set(GLOB_F_FILES globals/pgslib_globals_module.F globals/pgslib_type_module.F)
fortran_preprocess_files(GLOB_SOURCE_FILES
                         FILES ${GLOB_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR})
list(APPEND PGSLibIface_SOURCE_FILES ${GLOB_SOURCE_FILES})

# indexing files
set(IDX_F_FILES
            indexing/pgslib_index_gid_module.F
	    indexing/pgslib_index_table_module.F
            indexing/pgslib_index_module.F)
fortran_preprocess_files(IDX_SOURCE_FILES
                         FILES ${IDX_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR})
list(APPEND PGSLibIface_SOURCE_FILES ${IDX_SOURCE_FILES})

# io files
set(IO_F_FILES
           io/pgslib_io_bcast_module.F
	   io/pgslib_io_dist_module.F
           io/pgslib_io_collate_module.F
	   io/pgslib_io_module.F)
fortran_preprocess_files(IO_SOURCE_FILES
                         FILES ${IO_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS  -I${PGSLibIface_SOURCE_DIR}/io)
list(APPEND PGSLibIface_SOURCE_FILES ${IO_SOURCE_FILES})

# misc files
set(MISC_F_FILES
             misc/pgslib_decomp_module.F
	     misc/pgslib_permute_module.F
             misc/pgslib_misc_module.F 
	     misc/pgslib_shift_module.F)
fortran_preprocess_files(MISC_SOURCE_FILES
                         FILES ${MISC_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS  -I${PGSLibIface_SOURCE_DIR}/misc)
list(APPEND PGSLibIface_SOURCE_FILES ${MISC_SOURCE_FILES})

# partition files
set(PART_F_FILES
             partitions/graph_partition_data_types.F
             partitions/partition_constants.F
	     partitions/partition_data_types.F
	     partitions/partition_module.F)
fortran_preprocess_files(PART_SOURCE_FILES
                         FILES ${PART_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR})
list(APPEND PGSLibIface_SOURCE_FILES ${PART_SOURCE_FILES})

# reductions files
set(REDX_F_FILES
             reductions/pgslib_red_minmaxloc_module.F
             reductions/pgslib_red_numeric_module.F
             reductions/pgslib_reductions_module.F)
fortran_preprocess_files(REDX_SOURCE_FILES
                         FILES ${REDX_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS  -I${PGSLibIface_SOURCE_DIR}/reductions)
list(APPEND PGSLibIface_SOURCE_FILES ${REDX_SOURCE_FILES})

# scans files
set(SCANS_F_FILES
              scans/pgslib_scan_seg_bit_module.F
              scans/pgslib_scan_seg_module.F 
              scans/pgslib_scan_module.F 
              scans/pgslib_scan_no_seg_module.F)
fortran_preprocess_files(SCANS_SOURCE_FILES
                         FILES ${SCANS_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS  -I${PGSLibIface_SOURCE_DIR}/scans)
list(APPEND PGSLibIface_SOURCE_FILES ${SCANS_SOURCE_FILES})

# sort files
set(SORT_F_FILES sort/pgslib_grade_module.F  sort/pgslib_sort_module.F)
fortran_preprocess_files(SORT_SOURCE_FILES
                         FILES ${SORT_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR})
list(APPEND PGSLibIface_SOURCE_FILES ${SORT_SOURCE_FILES})

# user files
set(USER_F_FILES
             user/pgslib_scatter_log.F
	     user/pgslib_scatter_sum_module.F
             user/pgslib_scatter_minmax.F
	     user/pgslib_user_gs_module.F
             user/pgslib_gather_module.F
	     user/pgslib_scatter_module.F)
fortran_preprocess_files(USER_SOURCE_FILES
                         FILES ${USER_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS  -I${PGSLibIface_SOURCE_DIR}/user)
list(APPEND PGSLibIface_SOURCE_FILES ${USER_SOURCE_FILES})

# utility files
set(UTIL_F_FILES
             utility/pgslib_instrument.F
	     utility/pgslib_stats.F
             utility/pgslib_error_module.F
	     utility/pgslib_misc_utility.F
	     utility/pgslib_timing_module.F
             utility/pgslib_init.F
	     utility/pgslib_process_info.F
	     utility/pgslib_utility_module.F)
fortran_preprocess_files(UTIL_SOURCE_FILES
                         FILES ${UTIL_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS  -I${PGSLibIface_SOURCE_DIR}/utility)
list(APPEND PGSLibIface_SOURCE_FILES ${UTIL_SOURCE_FILES})

# main (module) file
set(MAIN_F_FILES main/pgslib_module.F)
fortran_preprocess_files(MAIN_SOURCE_FILES
                         FILES ${MAIN_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR})
list(APPEND PGSLibIface_SOURCE_FILES ${MAIN_SOURCE_FILES})


# Additional C source files
list(APPEND PGSLibIface_SOURCE_FILES utility/get_process_info.c)

# --- Library

# Force the preprocessing to occur before building the library
add_library(pgslib_fe ${PGSLibIface_SOURCE_FILES})
add_custom_target(ProcessPGSLibIfaceFiles ALL DEPENDS ${PGSLibIface_SOURCE_FILES})
add_dependencies(pgslib_fe ProcessPGSLibIfaceFiles)

# Set the library properties
set_target_properties(pgslib_fe PROPERTIES
                      OUTPUT_NAME       pgslib_fe
		      OUTPUT_NAME_DEBUG pgslib_fe_g
		      LINK_LANGUAGE     Fortran)

# Define the installation location
install(TARGETS pgslib_fe
        EXPORT  pgslib_fe
        LIBRARY DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR})


