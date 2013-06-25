/* This include translates names from base form into system dependent form.
   The most common translation is to append an underscore (_) to all
   routine names.  This is required on some systems, since on some sytems
   the fortran compiler appends an _ to all external names. */

/* $Id: pgslib-names.h,v 1.2 2002/09/12 20:52:34 lally Exp $ */

#ifndef PGSLIB_NAMES_H__
#define PGSLIB_NAMES_H__

/* This file has two parts.  The first is the translator macro.  The second
   is the definition of all the routines names.*/

/**********************************************************************/
/* Translator Macro                                                   */
/**********************************************************************/

/* FortranCInterface_names.h is created by CMake. It contians
   macros that manage the Fortran to C name mangling. The
   macro TR_ROUTINE_GLOBAL_ handles global routines with underscores
   in the name */
#include <FortranCInterface_names.h>


/**********************************************************************/
/* List of all C routines in PGSLib, by source directory              */
/**********************************************************************/

/* Routines in gath-scatt directory */
#define pgslib_setup_n_duplicate_c	TR_ROUTINE_GLOBAL_(pgslib_setup_n_duplicate_c,PGSLIB_SETUP_N_DUPLICATE_C)
#define pgslib_setup_duplicate_buffer_c	TR_ROUTINE_GLOBAL_(pgslib_setup_duplicate_buffer_c,PGSLIB_SETUP_DUPLICATE_BUFFER_C)
#define pgslib_prep_supplement_c	TR_ROUTINE_GLOBAL_(pgslib_prep_supplement_c,PGSLIB_PREP_SUPPLEMENT_C)
#define pgslib_gs_init_trace_c		TR_ROUTINE_GLOBAL_(pgslib_gs_init_trace_c,PGSLIB_GS_INIT_TRACE_C)
#define pgslib_gs_release_trace_c	TR_ROUTINE_GLOBAL_(pgslib_gs_release_trace_c,PGSLIB_GS_RELEASE_TRACE_C)
#define pgslib_gather_buf_int_c		TR_ROUTINE_GLOBAL_(pgslib_gather_buf_int_c,PGSLIB_GATHER_BUF_INT_C)
#define pgslib_gather_buf_float_c	TR_ROUTINE_GLOBAL_(pgslib_gather_buf_float_c,PGSLIB_GATHER_BUF_FLOAT_C)
#define pgslib_gather_buf_double_c	TR_ROUTINE_GLOBAL_(pgslib_gather_buf_double_c,PGSLIB_GATHER_BUF_DOUBLE_C)
#define pgslib_gather_buf_log_c		TR_ROUTINE_GLOBAL_(pgslib_gather_buf_log_c,PGSLIB_GATHER_BUF_LOG_C)
#define pgslib_scatter_buf_int_c	TR_ROUTINE_GLOBAL_(pgslib_scatter_buf_int_c,PGSLIB_SCATTER_BUF_INT_C)
#define pgslib_scatter_buf_float_c	TR_ROUTINE_GLOBAL_(pgslib_scatter_buf_float_c,PGSLIB_SCATTER_BUF_FLOAT_C)
#define pgslib_scatter_buf_double_c	TR_ROUTINE_GLOBAL_(pgslib_scatter_buf_double_c,PGSLIB_SCATTER_BUF_DOUBLE_C)
#define pgslib_scatter_buf_log_c	TR_ROUTINE_GLOBAL_(pgslib_scatter_buf_log_c,PGSLIB_SCATTER_BUF_LOG_C)

#define pgslib_trace_degree_c	        TR_ROUTINE_GLOBAL_(pgslib_trace_degree_c,PGSLIB_TRACE_DEGREE_C)

/* Routines in indexing directory */
#define pgslib_init_access_table_c	TR_ROUTINE_GLOBAL_(pgslib_init_access_table_c,PGSLIB_INIT_ACCESS_TABLE_C)
#define pgslib_free_access_table_c	TR_ROUTINE_GLOBAL_(pgslib_free_access_table_c,PGSLIB_FREE_ACCESS_TABLE_C)
#define pgslib_add_item_to_table_c	TR_ROUTINE_GLOBAL_(pgslib_add_item_to_table_c,PGSLIB_ADD_ITEM_TO_TABLE_C)
#define pgslib_count_items_in_table_c	TR_ROUTINE_GLOBAL_(pgslib_count_items_in_table_c,PGSLIB_COUNT_ITEMS_IN_TABLE_C)
#define pgslib_items_from_table_c	TR_ROUTINE_GLOBAL_(pgslib_items_from_table_c,PGSLIB_ITEMS_FROM_TABLE_C)
#define pgslib_item_index_from_table_c	TR_ROUTINE_GLOBAL_(pgslib_item_index_from_table_c,PGSLIB_ITEM_INDEX_FROM_TABLE_C)

/* Routines in io directory */
#define pgslib_bcast_int_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_int_scalar_c,PGSLIB_BCAST_INT_SCALAR_C)
#define pgslib_bcast_float_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_float_scalar_c,PGSLIB_BCAST_FLOAT_SCALAR_C)
#define pgslib_bcast_double_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_double_scalar_c,PGSLIB_BCAST_DOUBLE_SCALAR_C)
#define pgslib_bcast_log_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_log_scalar_c,PGSLIB_BCAST_LOG_SCALAR_C)
#define pgslib_bcast_char_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_char_scalar_c,PGSLIB_BCAST_CHAR_SCALAR_C)
#define pgslib_bcast_int_vector_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_int_vector_c,PGSLIB_BCAST_INT_VECTOR_C)
#define pgslib_bcast_float_vector_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_float_vector_c,PGSLIB_BCAST_FLOAT_VECTOR_C)
#define pgslib_bcast_double_vector_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_double_vector_c,PGSLIB_BCAST_DOUBLE_VECTOR_C)
#define pgslib_bcast_log_vector_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_log_vector_c,PGSLIB_BCAST_LOG_VECTOR_C)
#define pgslib_bcast_char_vector_c	TR_ROUTINE_GLOBAL_(pgslib_bcast_char_vector_c,PGSLIB_BCAST_CHAR_VECTOR_C)

#define pgslib_dist_int_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_dist_int_scalar_c,PGSLIB_DIST_INT_SCALAR_C)
#define pgslib_dist_float_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_dist_float_scalar_c,PGSLIB_DIST_FLOAT_SCALAR_C)
#define pgslib_dist_double_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_dist_double_scalar_c,PGSLIB_DIST_DOUBLE_SCALAR_C)
#define pgslib_dist_log_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_dist_log_scalar_c,PGSLIB_DIST_LOG_SCALAR_C)
#define pgslib_dist_int_vector_c	TR_ROUTINE_GLOBAL_(pgslib_dist_int_vector_c,PGSLIB_DIST_INT_VECTOR_C)
#define pgslib_dist_float_vector_c	TR_ROUTINE_GLOBAL_(pgslib_dist_float_vector_c,PGSLIB_DIST_FLOAT_VECTOR_C)
#define pgslib_dist_double_vector_c	TR_ROUTINE_GLOBAL_(pgslib_dist_double_vector_c,PGSLIB_DIST_DOUBLE_VECTOR_C)
#define pgslib_dist_log_vector_c	TR_ROUTINE_GLOBAL_(pgslib_dist_log_vector_c,PGSLIB_DIST_LOG_VECTOR_C)

#define pgslib_collate_int_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_collate_int_scalar_c,PGSLIB_COLLATE_INT_SCALAR_C)
#define pgslib_collate_float_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_collate_float_scalar_c,PGSLIB_COLLATE_FLOAT_SCALAR_C)
#define pgslib_collate_double_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_collate_double_scalar_c,PGSLIB_COLLATE_DOUBLE_SCALAR_C)
#define pgslib_collate_log_scalar_c	TR_ROUTINE_GLOBAL_(pgslib_collate_log_scalar_c,PGSLIB_COLLATE_LOG_SCALAR_C)
#define pgslib_collate_int_vector_c	TR_ROUTINE_GLOBAL_(pgslib_collate_int_vector_c,PGSLIB_COLLATE_INT_VECTOR_C)
#define pgslib_collate_float_vector_c	TR_ROUTINE_GLOBAL_(pgslib_collate_float_vector_c,PGSLIB_COLLATE_FLOAT_VECTOR_C)
#define pgslib_collate_double_vector_c	TR_ROUTINE_GLOBAL_(pgslib_collate_double_vector_c,PGSLIB_COLLATE_DOUBLE_VECTOR_C)
#define pgslib_collate_log_vector_c	TR_ROUTINE_GLOBAL_(pgslib_collate_log_vector_c,PGSLIB_COLLATE_LOG_VECTOR_C)
#define pgslib_collate_char_vector_c    TR_ROUTINE_GLOBAL_(pgslib_collate_char_vector_c,PGSLIB_COLLATE_CHAR_VECTOR_C)

/* Routines in reductions directory */
#define pgslib_global_min_int_c		TR_ROUTINE_GLOBAL_(pgslib_global_min_int_c,PGSLIB_GLOBAL_MIN_INT_C)
#define pgslib_global_min_float_c	TR_ROUTINE_GLOBAL_(pgslib_global_min_float_c,PGSLIB_GLOBAL_MIN_FLOAT_C)
#define pgslib_global_min_double_c	TR_ROUTINE_GLOBAL_(pgslib_global_min_double_c,PGSLIB_GLOBAL_MIN_DOUBLE_C)
#define pgslib_global_max_int_c		TR_ROUTINE_GLOBAL_(pgslib_global_max_int_c,PGSLIB_GLOBAL_MAX_INT_C)
#define pgslib_global_max_float_c	TR_ROUTINE_GLOBAL_(pgslib_global_max_float_c,PGSLIB_GLOBAL_MAX_FLOAT_C)
#define pgslib_global_max_double_c	TR_ROUTINE_GLOBAL_(pgslib_global_max_double_c,PGSLIB_GLOBAL_MAX_DOUBLE_C)
#define pgslib_global_sum_int_c		TR_ROUTINE_GLOBAL_(pgslib_global_sum_int_c,PGSLIB_GLOBAL_SUM_INT_C)
#define pgslib_global_sum_float_c	TR_ROUTINE_GLOBAL_(pgslib_global_sum_float_c,PGSLIB_GLOBAL_SUM_FLOAT_C)
#define pgslib_global_sum_double_c	TR_ROUTINE_GLOBAL_(pgslib_global_sum_double_c,PGSLIB_GLOBAL_SUM_DOUBLE_C)
#define pgslib_global_all_log_c		TR_ROUTINE_GLOBAL_(pgslib_global_all_log_c,PGSLIB_GLOBAL_ALL_LOG_C)
#define pgslib_global_any_log_c		TR_ROUTINE_GLOBAL_(pgslib_global_any_log_c,PGSLIB_GLOBAL_ANY_LOG_C)
#define pgslib_global_minloc_int_c	TR_ROUTINE_GLOBAL_(pgslib_global_minloc_int_c,PGSLIB_GLOBAL_MINLOC_INT_C)
#define pgslib_global_minloc_float_c	TR_ROUTINE_GLOBAL_(pgslib_global_minloc_float_c,PGSLIB_GLOBAL_MINLOC_FLOAT_C)
#define pgslib_global_minloc_double_c	TR_ROUTINE_GLOBAL_(pgslib_global_minloc_double_c,PGSLIB_GLOBAL_MINLOC_DOUBLE_C)
#define pgslib_global_maxloc_int_c	TR_ROUTINE_GLOBAL_(pgslib_global_maxloc_int_c,PGSLIB_GLOBAL_MAXLOC_INT_C)
#define pgslib_global_maxloc_float_c	TR_ROUTINE_GLOBAL_(pgslib_global_maxloc_float_c,PGSLIB_GLOBAL_MAXLOC_FLOAT_C)
#define pgslib_global_maxloc_double_c	TR_ROUTINE_GLOBAL_(pgslib_global_maxloc_double_c,PGSLIB_GLOBAL_MAXLOC_DOUBLE_C)

/* Routines in the scans directory */
#define off_node_sum_prefix_int_c	TR_ROUTINE_GLOBAL_(off_node_sum_prefix_int_c,OFF_NODE_SUM_PREFIX_INT_C)
#define off_node_sum_prefix_float_c	TR_ROUTINE_GLOBAL_(off_node_sum_prefix_float_c,OFF_NODE_SUM_PREFIX_FLOAT_C)
#define off_node_sum_prefix_double_c	TR_ROUTINE_GLOBAL_(off_node_sum_prefix_double_c,OFF_NODE_SUM_PREFIX_DOUBLE_C)
#define off_node_sum_suffix_int_c	TR_ROUTINE_GLOBAL_(off_node_sum_suffix_int_c,OFF_NODE_SUM_SUFFIX_INT_C)
#define off_node_sum_suffix_float_c	TR_ROUTINE_GLOBAL_(off_node_sum_suffix_float_c,OFF_NODE_SUM_SUFFIX_FLOAT_C)
#define off_node_sum_suffix_double_c	TR_ROUTINE_GLOBAL_(off_node_sum_suffix_double_c,OFF_NODE_SUM_SUFFIX_DOUBLE_C)

/* Routine in send-rcv directory */
#define pgslib_attempt_receive_int_c	TR_ROUTINE_GLOBAL_(pgslib_attempt_receive_int_c,PGSLIB_ATTEMPT_RECEIVE_INT_C)
#define pgslib_attempt_receive_float_c	TR_ROUTINE_GLOBAL_(pgslib_attempt_receive_float_c,PGSLIB_ATTEMPT_RECEIVE_FLOAT_C)
#define pgslib_attempt_receive_double_c	TR_ROUTINE_GLOBAL_(pgslib_attempt_receive_double_c,PGSLIB_ATTEMPT_RECEIVE_DOUBLE_C)
#define pgslib_attempt_receive_log_c	TR_ROUTINE_GLOBAL_(pgslib_attempt_receive_log_c,PGSLIB_ATTEMPT_RECEIVE_LOG_C)
#define pgslib_cnstd_send_rcv_int_c	TR_ROUTINE_GLOBAL_(pgslib_cnstd_send_rcv_int_c,PGSLIB_CNSTD_SEND_RCV_INT_C)
#define pgslib_cnstd_send_rcv_float_c	TR_ROUTINE_GLOBAL_(pgslib_cnstd_send_rcv_float_c,PGSLIB_CNSTD_SEND_RCV_FLOAT_C)
#define pgslib_cnstd_send_rcv_double_c	TR_ROUTINE_GLOBAL_(pgslib_cnstd_send_rcv_double_c,PGSLIB_CNSTD_SEND_RCV_DOUBLE_C)
#define pgslib_cnstd_send_rcv_log_c	TR_ROUTINE_GLOBAL_(pgslib_cnstd_send_rcv_log_c,PGSLIB_CNSTD_SEND_RCV_LOG_C)
#define pgslib_random_send_rcv_int_c	TR_ROUTINE_GLOBAL_(pgslib_random_send_rcv_int_c,PGSLIB_RANDOM_SEND_RCV_INT_C)
#define pgslib_random_send_rcv_float_c	TR_ROUTINE_GLOBAL_(pgslib_random_send_rcv_float_c,PGSLIB_RANDOM_SEND_RCV_FLOAT_C)
#define pgslib_random_send_rcv_double_c	TR_ROUTINE_GLOBAL_(pgslib_random_send_rcv_double_c,PGSLIB_RANDOM_SEND_RCV_DOUBLE_C)
#define pgslib_random_send_rcv_log_c	TR_ROUTINE_GLOBAL_(pgslib_random_send_rcv_log_c,PGSLIB_RANDOM_SEND_RCV_LOG_C)

/* Routines in the scans directory */
#define physical_c_shift_up_int_c	TR_ROUTINE_GLOBAL_(physical_c_shift_up_int_c,PHYSICAL_C_SHIFT_UP_INT_C)
#define physical_c_shift_up_float_c	TR_ROUTINE_GLOBAL_(physical_c_shift_up_float_c,PHYSICAL_C_SHIFT_UP_FLOAT_C)
#define physical_c_shift_up_double_c	TR_ROUTINE_GLOBAL_(physical_c_shift_up_double_c,PHYSICAL_C_SHIFT_UP_DOUBLE_C)
#define physical_c_shift_down_int_c	TR_ROUTINE_GLOBAL_(physical_c_shift_down_int_c,PHYSICAL_C_SHIFT_DOWN_INT_C)
#define physical_c_shift_down_float_c	TR_ROUTINE_GLOBAL_(physical_c_shift_down_float_c,PHYSICAL_C_SHIFT_DOWN_FLOAT_C)
#define physical_c_shift_down_double_c	TR_ROUTINE_GLOBAL_(physical_c_shift_down_double_c,PHYSICAL_C_SHIFT_DOWN_DOUBLE_C)

/* Routines in utility directory */
#define pgslib_initialize_c		TR_ROUTINE_GLOBAL_(pgslib_initialize_c,PGSLIB_INITIALIZE_C)
#define pgslib_mpi_init			TR_ROUTINE_GLOBAL_(pgslib_mpi_init,PGSLIB_MPI_INIT)
#define pgslib_finalize_c		TR_ROUTINE_GLOBAL_(pgslib_finalize_c,PGSLIB_FINALIZE_C)
#define pgslib_error_c			TR_ROUTINE_GLOBAL_(pgslib_error_c,PGSLIB_ERROR_C)
#define pgslib_fatal_error_c		TR_ROUTINE_GLOBAL_(pgslib_fatal_error_c,PGSLIB_FATAL_ERROR_C)
#define pgslib_abort_c			TR_ROUTINE_GLOBAL_(pgslib_abort_c,PGSLIB_ABORT_C)
#define pgslib_output_c			TR_ROUTINE_GLOBAL_(pgslib_output_c,PGSLIB_OUTPUT_C)
#define pgslib_flush_output_c		TR_ROUTINE_GLOBAL_(pgslib_flush_output_c,PGSLIB_FLUSH_OUTPUT_C)
#define pgslib_close_output_c		TR_ROUTINE_GLOBAL_(pgslib_close_output_c,PGSLIB_CLOSE_OUTPUT_C)
#define pgslib_check_malloc_c		TR_ROUTINE_GLOBAL_(pgslib_check_malloc_c,PGSLIB_CHECK_MALLOC_C)
#define pgslib_barrier_c		TR_ROUTINE_GLOBAL_(pgslib_barrier_c,PGSLIB_BARRIER_C)

#define pgslib_get_argc			TR_ROUTINE_GLOBAL_(pgslib_get_argc,PGSLIB_GET_ARGC)
#define pgslib_get_argv			TR_ROUTINE_GLOBAL_(pgslib_get_argv,PGSLIB_GET_ARGV)
#define pgslib_cleanup_argv		TR_ROUTINE_GLOBAL_(pgslib_cleanup_argv,PGSLIB_CLEANUP_ARGV)

#define pgslib_barrier_time_c           TR_ROUTINE_GLOBAL_(pgslib_barrier_time_c,PGSLIB_BARRIER_TIME_C)
#define pgslib_sr_time_c           	TR_ROUTINE_GLOBAL_(pgslib_sr_time_c,PGSLIB_SR_TIME_C)

#endif
