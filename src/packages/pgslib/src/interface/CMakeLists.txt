project(PGSLibIface C Fortran)

# --- Include directories
#include_directories(${PGSLib_MODULE_DIR})

# --- Define *.F -> *.f90 files
set(PGSLibIface_SOURCE_FILES pgslib_c_binding.F90)

# gath-scatt files
set(GS_SOURCE_FILES
           gath-scatt/pgslib_gs_comm_module.F90
	   gath-scatt/pgslib_gs_setup_module.F90
           gath-scatt/pgslib_gs_module.F90
	   gath-scatt/pgslib_gs_util_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${GS_SOURCE_FILES})

# global files
set(GLOB_SOURCE_FILES
  globals/pgslib_globals_module.F90
  globals/pgslib_type_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${GLOB_SOURCE_FILES})

# indexing files
set(IDX_SOURCE_FILES
            indexing/pgslib_index_gid_module.F90
	    indexing/pgslib_index_table_module.F90
            indexing/pgslib_index_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${IDX_SOURCE_FILES})

# io files
set(IO_SOURCE_FILES
           io/pgslib_io_bcast_module.F90
	   io/pgslib_io_dist_module.F90
           io/pgslib_io_collate_module.F90
	   io/pgslib_io_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${IO_SOURCE_FILES})

# misc files
set(MISC_SOURCE_FILES
             misc/pgslib_decomp_module.F90
	     misc/pgslib_permute_module.F90
             misc/pgslib_misc_module.F90
	     misc/pgslib_shift_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${MISC_SOURCE_FILES})

# partition files
set(PART_SOURCE_FILES
             partitions/graph_partition_data_types.F90
             partitions/partition_constants.F90
	     partitions/partition_data_types.F90
	     partitions/partition_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${PART_SOURCE_FILES})

# reductions files
set(REDX_SOURCE_FILES
             reductions/pgslib_red_minmaxloc_module.F90
             reductions/pgslib_red_numeric_module.F90
             reductions/pgslib_reductions_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${REDX_SOURCE_FILES})

# scans files
set(SCANS_SOURCE_FILES
              scans/pgslib_scan_seg_bit_module.F90
              scans/pgslib_scan_seg_module.F90
              scans/pgslib_scan_module.F90
              scans/pgslib_scan_no_seg_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${SCANS_SOURCE_FILES})

# sort files
set(SORT_SOURCE_FILES
  sort/pgslib_grade_module.F90
  sort/pgslib_sort_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${SORT_SOURCE_FILES})

# user files
set(USER_SOURCE_FILES
             user/pgslib_scatter_log.F90
	     user/pgslib_scatter_sum_module.F90
             user/pgslib_scatter_minmax.F90
	     user/pgslib_user_gs_module.F90
             user/pgslib_gather_module.F90
	     user/pgslib_scatter_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${USER_SOURCE_FILES})

# utility files
set(UTIL_SOURCE_FILES
             utility/pgslib_instrument.F90
	     utility/pgslib_stats.F90
             utility/pgslib_error_module.F90
	     utility/pgslib_misc_utility.F90
	     utility/pgslib_timing_module.F90
             utility/pgslib_init.F90
	     utility/pgslib_process_info.F90
	     utility/pgslib_utility_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${UTIL_SOURCE_FILES})

# main (module) file
set(MAIN_SOURCE_FILES main/pgslib_module.F90)
list(APPEND PGSLibIface_SOURCE_FILES ${MAIN_SOURCE_FILES})


# Additional C source files
list(APPEND PGSLibIface_SOURCE_FILES utility/get_process_info.c)

# --- Library

# Force the preprocessing to occur before building the library
add_library(pgslib_fe ${PGSLibIface_SOURCE_FILES})
# - Pick up the module directory
target_include_directories(pgslib_fe PUBLIC $<TARGET_PROPERTY:pgslib_fe,Fortran_MODULE_DIRECTORY>)


add_library(truchas::pgslib_fe ALIAS pgslib_fe)

# Set the library properties
set_target_properties(pgslib_fe PROPERTIES
                      OUTPUT_NAME       pgslib_fe
		      OUTPUT_NAME_DEBUG pgslib_fe_g
		      LINKER_LANGUAGE     Fortran)

target_link_libraries(pgslib_fe truchas::pgslib_mpi)

# Define the installation location
install(TARGETS pgslib_fe
        EXPORT  pgslib_fe
        LIBRARY DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${PGSLib_LIBRARY_INSTALL_DIR})
