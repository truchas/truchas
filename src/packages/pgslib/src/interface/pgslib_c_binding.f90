!!
!! PGSLIB_C_BINDING
!!
!! Bindings to the low-level C code used by the top-level Fortran interface.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! 1. The lower level C component uses ints to represent logical values.
!! Fortran logical variables are not interoperable with C int variables,
!! however, and so we must resort to old-school hackery to interface Fortran
!! logical arguments to C int arguments.  This involves handling two key
!! issues: the data size and the representation of logical values (see Note 2).
!! Here we make the assumption (potentially not portable) that integer and
!! logical types with the same type parameter are the same size; that is,
!! logical(c_int) and integer(c_int) are the same size.  Then in the interface
!! where a "logical" argument is being passed we lie and declare it as
!! logical(c_int) to match the actual Fortran usage (the compiler will see it
!! as an error if we did otherwise) even though it ought to be declared as
!! integer(c_int) to match the true C function argument.
!!
!! 2. While Fortran uses integer-like data elements to store logical values,
!! it does not specify how logical values are represented, and thus logical
!! values are not generally interoperable with C.  [Intel Fortran, for example,
!! uses the state of a single bit by default, so that an even value is false
!! and an odd value is true.  This can be changed to use the C representation
!! with the compiler flag "-fpscomp logicals".]  Here we assume that the C
!! functions are merely moving/copying integer data around and never attempt
!! to interpret the data as logical values, and thus are free to entirely
!! ignore the representation issue.  The only exception is for the functions
!! pgslib_global_{any,all}_log_c.  Here the "logical" arguments are declared
!! as integers, as they truly are, and the client code explicitly converts
!! its logical values to corresponding C integer values that are passed to
!! the functions.
!!
!! 3. These interface blocks are doing two things: defining the interface
!! to external C functions and defining a generic procedure that has them
!! as its specifics.  I really do not care for this style, and would prefer
!! to define the generic in a separate interface block (which would refer
!! only to the specific procedure names).  However, the Intel 14.0 compiler
!! has a bug that prevents one from doing this (fixed in 15.0).
!!
!! 4. In the Fortran clients it is permissible, and often the case, for data
!! arrays to be 0-sized.  Such arrays are not strictly interoperable with C
!! (being addressed in F2015).  We have entirely ignored that issue and
!! assumed that the generated code will gracefully handle deciding what the
!! address of a 0-sized array is (this is what actually gets passed to C).
!! It could pass anything, as the C code will know via other arguments that
!! that address is a dummy and will not dereference it.
!!

module pgslib_c_binding

  use,intrinsic :: iso_c_binding, only: c_ptr, c_int, c_float, c_double, c_char
  implicit none
  public

  !! Functions from par-src/utility/utility-c.c
  interface
    subroutine PGSLib_Initialize_C (nPE, thisPE, IO_ROOT_PE, File_Per_PE, File_Prefix) &
        bind(c, name="pgslib_initialize_c")
      import c_int, c_char
      integer(c_int), intent(out) :: nPE, thisPE
      integer(c_int), intent(inout) :: IO_ROOT_PE
      integer(c_int), intent(in) :: File_Per_PE
      character(kind=c_char), intent(in) :: File_Prefix(*)
    end subroutine
    subroutine PGSLib_Finalize_C () bind(c, name="pgslib_finalize_c")
    end subroutine
    subroutine pgslib_get_argc (argc) bind(c)
      import c_int
      integer(c_int) :: argc
    end subroutine
    subroutine pgslib_get_argv (i, string_l, string) bind(c)
      import c_int, c_char
      integer(c_int) :: i, string_l
      character(kind=c_char) :: string(*)
    end subroutine
    subroutine pgslib_cleanup_argv () bind(c)
    end subroutine
    subroutine PGSLib_Output_C (Message) bind(c, name="pgslib_output_c")
      import c_char
      character(kind=c_char), intent(in) :: Message(*)
    end subroutine
    subroutine PGSLib_Error_C (Estring) bind(c, name="pgslib_error_c")
      import c_char
      character(kind=c_char) :: Estring(*)
    end subroutine
    subroutine PGSLib_Abort_C() bind(c, name="pgslib_abort_c")
    end subroutine
    subroutine PGSLib_Flush_Output_C() bind(c, name="pgslib_flush_output_c")
    end subroutine
    subroutine PGSLib_Barrier_C() bind(c, name="pgslib_barrier_c")
    end subroutine
    subroutine PGSLib_Barrier_Time_C (BT) bind(c, name="pgslib_barrier_time_c")
      import c_float
      real(c_float) :: bt
    end subroutine
    subroutine PGSLib_SR_Time_C (t) bind(c, name="pgslib_sr_time_c")
      import c_float
      real(c_float) :: t
    end subroutine
  end interface

  !! Functions from par-src/gath-scatt/gs-util-c.c
  interface
    subroutine PGSLib_GS_Init_Trace_C (GS_Trace) &
        bind(c, name="pgslib_gs_init_trace_c")
      import c_ptr
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_GS_Release_Trace_C (GS_Trace) &
        bind(c, name="pgslib_gs_release_trace_c")
      import c_ptr
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Trace_Degree_C (Scatter_Degree, Gather_Degree, GS_Trace) &
        bind(c, name="pgslib_trace_degree_c")
      import c_int, c_ptr
      integer(c_int) :: Scatter_Degree, Gather_Degree
      type(c_ptr) :: GS_Trace
    end subroutine
  end interface

  !! Functions from par-src/gath-scatt/gs-setup-c.c
  interface
    subroutine PGSLib_Setup_N_Duplicate_C (N_Duplicate, GS_Trace) &
        bind(c, name="pgslib_setup_n_duplicate_c")
      import c_int, c_ptr
      integer(c_int) :: N_Duplicate
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Prep_Supplement_C (N_Supplement, PEs, GS_Trace) &
        bind(c, name="pgslib_prep_supplement_c")
      import c_int, c_ptr
      integer(c_int) :: N_Supplement, PEs(*)
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Setup_Duplicate_Buffer_C (GS_Trace) &
        bind(c, name="pgslib_setup_duplicate_buffer_c")
      import c_ptr
      type(c_ptr) :: GS_trace
    end subroutine
  end interface

  !! Functions from par-src/gath-scatt/gather-*.c
  interface
    subroutine PGSLib_Gather_Buf_INT_C (Supplement, Duplicate, BlockSize, GS_Trace) &
        bind(c, name="pgslib_gather_buf_int_c")
      import c_ptr, c_int
      integer(c_int) :: Supplement(*), Duplicate(*)
      integer(c_int) :: BlockSize
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Gather_Buf_FLOAT_C (Supplement, Duplicate, BlockSize, GS_Trace) &
        bind(c, name="pgslib_gather_buf_float_c")
      import c_ptr, c_int, c_float
      real(c_float) :: Supplement(*), Duplicate(*)
      integer(c_int) :: BlockSize
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Gather_Buf_DOUBLE_C (Supplement, Duplicate, BlockSize, GS_Trace) &
        bind(c, name="pgslib_gather_buf_double_c")
      import c_ptr, c_int, c_double
      real(c_double) :: Supplement(*), Duplicate(*)
      integer(c_int) :: BlockSize
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Gather_Buf_LOG_C (Supplement, Duplicate, BlockSize, GS_Trace) &
        bind(c, name="pgslib_gather_buf_log_c")
      import c_ptr, c_int
      logical(c_int) :: Supplement(*), Duplicate(*) ! A LIE: SEE NOTE 1
      integer(c_int) :: BlockSize
      type(c_ptr) :: GS_Trace
    end subroutine
  end interface

  !! Functions from par-src/gath-scatt/scatter-*.c
  interface
    subroutine PGSLib_Scatter_Buf_INT_C (Duplicate, Supplement, BlockSize, GS_Trace) &
        bind(c, name="pgslib_scatter_buf_int_c")
      import c_ptr, c_int
      integer(c_int) :: Duplicate(*), Supplement(*)
      integer(c_int) :: BlockSize
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Scatter_Buf_FLOAT_C (Duplicate, Supplement, BlockSize, GS_Trace) &
        bind(c, name="pgslib_scatter_buf_float_c")
      import c_ptr, c_int, c_float
      real(c_float) :: Duplicate(*), Supplement(*)
      integer(c_int) :: BlockSize
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Scatter_Buf_DOUBLE_C (Duplicate, Supplement, BlockSize, GS_Trace) &
        bind(c, name="pgslib_scatter_buf_double_c")
      import c_ptr, c_int, c_double
      real(c_double) :: Duplicate(*), Supplement(*)
      integer(c_int) :: BlockSize
      type(c_ptr) :: GS_Trace
    end subroutine
    subroutine PGSLib_Scatter_Buf_LOG_C (Duplicate, Supplement, BlockSize, GS_Trace) &
        bind(c, name="pgslib_scatter_buf_log_c")
      import c_ptr, c_int
      logical(c_int) :: Duplicate(*), Supplement(*) ! A LIE: SEE NOTE 1
      integer(c_int) :: BlockSize
      type(c_ptr) :: GS_Trace
    end subroutine
  end interface

  !! Functions from par-src/io/io-c.c
  interface
    subroutine PGSLib_bcast_int_scalar_c (scalar) &
        bind(c, name="pgslib_bcast_int_scalar_c")
      import c_int
      integer(c_int) :: scalar
    end subroutine
    subroutine PGSLib_bcast_float_scalar_c (scalar) &
        bind(c, name="pgslib_bcast_float_scalar_c")
      import c_float
      real(c_float) :: scalar
    end subroutine
    subroutine PGSLib_bcast_double_scalar_c (scalar) &
        bind(c, name="pgslib_bcast_double_scalar_c")
      import c_double
      real(c_double) :: scalar
    end subroutine
    subroutine PGSLib_bcast_log_scalar_c (scalar) &
        bind(c, name="pgslib_bcast_log_scalar_c")
      import c_int
      logical(c_int) :: scalar ! A LIE: SEE NOTE 1
    end subroutine
    subroutine PGSLib_bcast_int_vector_c (vector, vec_len) &
        bind(c, name="pgslib_bcast_int_vector_c")
      import c_int
      integer(c_int) :: vector(*)
      integer(c_int), intent(in) :: vec_len
    end subroutine
    subroutine PGSLib_bcast_float_vector_c (vector, vec_len) &
        bind(c, name="pgslib_bcast_float_vector_c")
      import c_int, c_float
      real(c_float) :: vector(*)
      integer(c_int), intent(in) :: vec_len
    end subroutine
    subroutine PGSLib_bcast_double_vector_c (vector, vec_len) &
        bind(c, name="pgslib_bcast_double_vector_c")
      import c_int, c_double
      real(c_double) :: vector(*)
      integer(c_int), intent(in) :: vec_len
    end subroutine
    subroutine PGSLib_bcast_log_vector_c (vector, vec_len) &
        bind(c, name="pgslib_bcast_log_vector_c")
      import c_int
      logical(c_int) :: vector(*) ! A LIE: SEE NOTE 1
      integer(c_int), intent(in) :: vec_len
    end subroutine
    subroutine PGSLib_bcast_char_vector_c (vector, vec_len) &
        bind(c, name="pgslib_bcast_char_vector_c")
      import c_int, c_char
      character(kind=c_char) :: vector(*)
      integer(c_int), intent(in) :: vec_len
    end subroutine
  end interface

  !! Functions from par-src/io/io-c.c
  interface PGSLib_Collate_Scalar_C ! SEE NOTE 3
    subroutine PGSLib_Collate_int_scalar_c (scalarv_out, scalar_in) &
        bind(c, name="pgslib_collate_int_scalar_c")
      import c_int
      integer(c_int), intent(out) :: scalarv_out(*)
      integer(c_int), intent(in)  :: scalar_in
    end subroutine
    subroutine PGSLib_Collate_float_scalar_c (scalarv_out, scalar_in) &
        bind(c, name="pgslib_collate_float_scalar_c")
      import c_float
      real(c_float), intent(out) :: scalarv_out(*)
      real(c_float), intent(in)  :: scalar_in
    end subroutine
    subroutine PGSLib_Collate_double_scalar_c (scalarv_out, scalar_in) &
        bind(c, name="pgslib_collate_double_scalar_c")
      import c_double
      real(c_double), intent(out) :: scalarv_out(*)
      real(c_double), intent(in)  :: scalar_in
    end subroutine
    subroutine PGSLib_Collate_log_scalar_c (scalarv_out, scalar_in) &
        bind(c, name="pgslib_collate_log_scalar_c")
      import c_int
      logical(c_int), intent(out) :: scalarv_out(*) ! A LIE: SEE NOTE 1
      logical(c_int), intent(in)  :: scalar_in      ! A LIE: SEE NOTE 1
    end subroutine
  end interface
  private :: PGSLib_Collate_int_scalar_c, PGSLib_Collate_float_scalar_c, &
             PGSLib_Collate_double_scalar_c, PGSLib_Collate_log_scalar_c

  !! Functions from par-src/io/io-c.c
  interface PGSLib_Collate_Vector_c ! SEE NOTE 3
    subroutine PGSLib_Collate_int_vector_c (vector_out, Lengths, vector_in, In_Length) &
        bind(c, name="pgslib_collate_int_vector_c")
      import c_int
      integer(c_int), intent(out) :: vector_out(*)
      integer(c_int), intent(in)  :: vector_in(*)
      integer(c_int), intent(in)  :: Lengths(*), In_Length
    end subroutine
    subroutine PGSLib_Collate_float_vector_c (vector_out, Lengths, vector_in, In_Length) &
        bind(c, name="pgslib_collate_float_vector_c")
      import c_int, c_float
      real(c_float),  intent(out) :: vector_out(*)
      real(c_float),  intent(in)  :: vector_in(*)
      integer(c_int), intent(in)  :: Lengths(*), In_Length
    end subroutine
    subroutine PGSLib_Collate_double_vector_c (vector_out, Lengths, vector_in, In_Length) &
        bind(c, name="pgslib_collate_double_vector_c")
      import c_int, c_double
      real(c_double), intent(out) :: vector_out(*)
      real(c_double), intent(in)  :: vector_in(*)
      integer(c_int), intent(in)  :: Lengths(*), In_Length
    end subroutine
    subroutine PGSLib_Collate_log_vector_c (vector_out, Lengths, vector_in, In_Length) &
        bind(c, name="pgslib_collate_log_vector_c")
      import c_int
      logical(c_int), intent(out) :: vector_out(*) ! A LIE: SEE NOTE 1
      logical(c_int), intent(in)  :: vector_in(*)  ! A LIE: SEE NOTE 1
      integer(c_int), intent(in)  :: Lengths(*), In_Length
    end subroutine
    subroutine PGSLib_Collate_char_vector_c (vector_out, Lengths, vector_in, In_Length) &
        bind(c, name="pgslib_collate_char_vector_c")
      import c_int, c_char
      character(kind=c_char), intent(out) :: vector_out(*)
      character(kind=c_char), intent(in)  :: vector_in(*)
      integer(c_int), intent(in) :: Lengths(*), In_Length
    end subroutine
  end interface
  private :: PGSLib_Collate_int_vector_c, PGSLib_Collate_float_vector_c, &
             PGSLib_Collate_double_vector_c, PGSLib_Collate_log_vector_c, &
             PGSLib_Collate_char_vector_c

  !! Functions from par-src/io/io-c.c
  interface PGSLib_Dist_Scalar_c  ! SEE NOTE 3
    subroutine PGSLib_Dist_int_scalar_c (scalar_out, scalarv_in) &
        bind(c, name="pgslib_dist_int_scalar_c")
      import c_int
      integer(c_int), intent(out) :: scalar_out
      integer(c_int), intent(in)  :: scalarv_in(*)
    end subroutine
    subroutine PGSLib_Dist_float_scalar_c (scalar_out, scalarv_in) &
        bind(c, name="pgslib_dist_float_scalar_c")
      import c_float
      real(c_float), intent(out) :: scalar_out
      real(c_float), intent(in)  :: scalarv_in(*)
    end subroutine
    subroutine PGSLib_Dist_double_scalar_c (scalar_out, scalarv_in) &
        bind(c, name="pgslib_dist_double_scalar_c")
      import c_double
      real(c_double), intent(out) :: scalar_out
      real(c_double), intent(in)  :: scalarv_in(*)
    end subroutine
    subroutine PGSLib_Dist_log_scalar_c (scalar_out, scalarv_in) &
        bind(c, name="pgslib_dist_log_scalar_c")
      import c_int
      logical(c_int), intent(out) :: scalar_out    ! A LIE: SEE NOTE 1
      logical(c_int), intent(in)  :: scalarv_in(*) ! A LIE: SEE NOTE 1
    end subroutine
  end interface
  private :: PGSLib_Dist_int_scalar_c, PGSLib_Dist_float_scalar_c, &
             PGSLib_Dist_double_scalar_c, PGSLib_Dist_log_scalar_c

  !! Functions from par-src/io/io-c.c
  interface PGSLib_Dist_Vector_c  ! SEE NOTE 3
    subroutine PGSLib_dist_int_vector_c (vector_out, out_len, vector_in, lengths) &
        bind(c, name="pgslib_dist_int_vector_c")
      import c_int
      integer(c_int), intent(out) :: vector_out(*)
      integer(c_int), intent(in)  :: vector_in(*)
      integer(c_int), intent(in)  :: out_len, lengths(*)
    end subroutine
    subroutine PGSLib_dist_float_vector_c (vector_out, out_len, vector_in, lengths) &
        bind(c, name="pgslib_dist_float_vector_c")
      import c_int, c_float
      real(c_float), intent(out) :: vector_out(*)
      real(c_float), intent(in)  :: vector_in(*)
      integer(c_int), intent(in) :: out_len, lengths(*)
    end subroutine
    subroutine PGSLib_dist_double_vector_c (vector_out, out_len, vector_in, lengths) &
        bind(c, name="pgslib_dist_double_vector_c")
      import c_int, c_double
      real(c_double), intent(out) :: vector_out(*)
      real(c_double), intent(in)  :: vector_in(*)
      integer(c_int), intent(in) :: out_len, lengths(*)
    end subroutine
    subroutine PGSLib_dist_log_vector_c (vector_out, out_len, vector_in, lengths) &
        bind(c, name="pgslib_dist_log_vector_c")
      import c_int
      logical(c_int), intent(out) :: vector_out(*) ! A LIE: SEE NOTE 1
      logical(c_int), intent(in)  :: vector_in(*)  ! A LIE: SEE NOTE 1
      integer(c_int), intent(in)  :: out_len, lengths(*)
    end subroutine
  end interface
  private :: PGSLib_dist_int_vector_c, PGSLib_dist_float_vector_c, &
             PGSLib_dist_double_vector_c, PGSLib_dist_log_vector_c

  !! Functions from par-src/scans/scan-c-*.c
  interface Off_Node_SUM_PREFIX ! SEE NOTE 3
    subroutine Off_Node_SUM_PREFIX_INT_C (Dest_Data, Dest_Seg, Src_Data, Src_Seg) &
        bind(c, name="off_node_sum_prefix_int_c")
      import c_int
      integer(c_int) :: Dest_Data, Dest_Seg, Src_Data, Src_Seg
    end subroutine
    subroutine Off_Node_SUM_PREFIX_FLOAT_C (Dest_Data, Dest_Seg, Src_Data, Src_Seg) &
        bind(c, name="off_node_sum_prefix_float_c")
      import c_int, c_float
      real(c_float) :: Dest_Data, Src_Data
      integer(c_int) :: Dest_Seg, Src_Seg
    end subroutine
    subroutine Off_Node_SUM_PREFIX_DOUBLE_C (Dest_Data, Dest_Seg, Src_Data, Src_Seg) &
        bind(c, name="off_node_sum_prefix_double_c")
      import c_int, c_double
      real(c_double) :: Dest_Data, Src_Data
      integer(c_int) :: Dest_Seg, Src_Seg
    end subroutine
  end interface
  private :: Off_Node_SUM_PREFIX_INT_C, Off_Node_SUM_PREFIX_FLOAT_C, Off_Node_SUM_PREFIX_DOUBLE_C

  !! Functions from par-src/scans/scan-c-*.c
  interface Off_Node_SUM_SUFFIX ! SEE NOTE 3
    subroutine Off_Node_SUM_SUFFIX_INT_C (Dest_Data, Dest_Seg, Src_Data, Src_Seg) &
        bind(c, name="off_node_sum_suffix_int_c")
      import c_int
      integer(c_int) :: Dest_Data, Dest_Seg, Src_Data, Src_Seg
    end subroutine
    subroutine Off_Node_SUM_SUFFIX_FLOAT_C (Dest_Data, Dest_Seg, Src_Data, Src_Seg) &
        bind(c, name="off_node_sum_suffix_float_c")
      import c_int, c_float
      real(c_float) :: Dest_Data, Src_Data
      integer(c_int) :: Dest_Seg, Src_Seg
    end subroutine
    subroutine Off_Node_SUM_SUFFIX_DOUBLE_C (Dest_Data, Dest_Seg, Src_Data, Src_Seg) &
        bind(c, name="off_node_sum_suffix_double_c")
      import c_int, c_double
      real(c_double) :: Dest_Data, Src_Data
      integer(c_int) :: Dest_Seg, Src_Seg
    end subroutine
  end interface
  private :: Off_Node_SUM_SUFFIX_INT_C, Off_Node_SUM_SUFFIX_FLOAT_C, Off_Node_SUM_SUFFIX_DOUBLE_C

  !! Functions from par-src/indexing/index-c.c
  interface
    subroutine PGSLib_Init_Access_Table_C (Access_Table) &
        bind(c, name="pgslib_init_access_table_c")
      import c_ptr
      type(c_ptr) :: Access_Table
    end subroutine
    subroutine PGSLib_Free_Access_Table_C (Access_Table) &
        bind(c, name="pgslib_free_access_table_c")
      import c_ptr
      type(c_ptr) :: Access_Table
    end subroutine
    subroutine PGSLib_Add_Item_To_Table_C (Item, PE, Access_Table, ierror) &
        bind(c, name="pgslib_add_item_to_table_c")
      import c_ptr, c_int
      integer(c_int), intent(in) :: Item, PE
      type(c_ptr) :: Access_Table
      integer(c_int), intent(out) :: ierror
    end subroutine
    subroutine PGSLib_Count_Items_In_Table_C (Count, Access_Table) &
        bind(c, name="pgslib_count_items_in_table_c")
      import c_ptr, c_int
      integer(c_int), intent(out) :: Count
      type(c_ptr) :: Access_Table
    end subroutine
    subroutine PGSLib_Items_From_Table_C (Items, PEs, Count, Access_Table, ierror) &
        bind(c, name="pgslib_items_from_table_c")
      import c_ptr, c_int
      integer(c_int), intent(out) :: Items(*), PEs(*), ierror
      integer(c_int), intent(in)  :: Count
      type(c_ptr) :: Access_Table
    end subroutine
    subroutine PGSLib_Item_Index_From_Table_C (Index, Item, PE, Access_Table) &
        bind(c, name="PGSLib_Item_Index_From_Table_C")
      import c_ptr, c_int
      integer(c_int), intent(out) :: Index
      integer(c_int), intent(in)  :: Item, PE
      type(c_ptr) :: Access_Table
    end subroutine
  end interface

  !! Functions from par-src/reductions/redux-c.c
  interface PGSLib_Global_MIN_c ! SEE NOTE 3
    subroutine PGSLib_global_min_int_c (MinA) &
        bind(c, name="pgslib_global_min_int_c")
      import c_int
      integer(c_int), intent(inout):: MinA
    end subroutine
    subroutine PGSLib_global_min_float_c (MinA) &
        bind(c, name="pgslib_global_min_float_c")
      import c_float
      real(c_float), intent(inout):: MinA
    end subroutine
    subroutine PGSLib_global_min_double_c (MinA) &
        bind(c, name="pgslib_global_min_double_c")
      import c_double
      real(c_double), intent(inout):: MinA
    end subroutine
  end interface
  private :: PGSLib_global_min_int_c, PGSLib_global_min_float_c, PGSLib_global_min_double_c

  !! Functions from par-src/reductions/redux-c.c
  interface PGSLib_Global_MAX_C ! SEE NOTE 3
    subroutine PGSLib_global_max_int_c (MaxA) &
        bind(c, name="pgslib_global_max_int_c")
      import c_int
      integer(c_int), intent(inout):: MaxA
    end subroutine
    subroutine PGSLib_global_max_float_c (MaxA) &
        bind(c, name="pgslib_global_max_float_c")
      import c_float
      real(c_float), intent(inout):: MaxA
    end subroutine
    subroutine PGSLib_global_max_double_c (MaxA) &
        bind(c, name="pgslib_global_max_double_c")
      import c_double
      real(c_double), intent(inout):: MaxA
    end subroutine
  end interface
  private :: PGSLib_global_max_int_c, PGSLib_global_max_float_c, PGSLib_global_max_double_c

  !! Functions from par-src/reductions/redux-c.c
  interface PGSLib_Global_SUM_C ! SEE NOTE 3
    subroutine PGSLib_global_sum_int_c (SumA) &
        bind(c, name="pgslib_global_sum_int_c")
      import c_int
      integer(c_int), intent(inout):: SumA
    end subroutine
    subroutine PGSLib_global_sum_float_c (SumA) &
        bind(c, name="pgslib_global_sum_float_c")
      import c_float
      real(c_float), intent(inout):: SumA
    end subroutine
    subroutine PGSLib_global_sum_double_c (SumA) &
        bind(c, name="pgslib_global_sum_double_c")
      import c_double
      real(c_double), intent(inout):: SumA
    end subroutine
  end interface
  private :: PGSLib_global_sum_int_c, PGSLib_global_sum_float_c, PGSLib_global_sum_double_c

  !! Functions from par-src/reductions/redux-c.c
  interface
    subroutine PGSLib_Global_ALL_C (AllA) &
        bind(c, name="pgslib_global_all_log_c")
      import c_int
      integer(c_int), intent(inout) :: AllA ! SEE NOTE 2
    end subroutine
    subroutine PGSLib_Global_ANY_C (AnyA) &
        bind(c, name="pgslib_global_any_log_c")
      import c_int
      integer(c_int), intent(inout) :: AnyA ! SEE NOTE 2
    end subroutine
  end interface

  !! Functions from par-src/reductions/redux-c.c
  interface PGSLib_Global_MINLOC_C  ! SEE NOTE 3
    subroutine PGSLib_Global_MINLOC_int_C (MinV, GlobalIndex) &
        bind(c, name="pgslib_global_minloc_int_c")
      import c_int
      integer(c_int), intent(inout) :: GlobalIndex
      integer(c_int), intent(inout) :: MinV
    end subroutine
    subroutine PGSLib_Global_MINLOC_float_C (MinV, GlobalIndex) &
        bind(c, name="pgslib_global_minloc_float_c")
      import c_int, c_float
      integer(c_int), intent(inout) :: GlobalIndex
      real(c_float),  intent(inout) :: MinV
    end subroutine
    subroutine PGSLib_Global_MINLOC_double_C (MinV, GlobalIndex) &
        bind(c, name="pgslib_global_minloc_double_c")
      import c_int, c_double
      integer(c_int), intent(inout) :: GlobalIndex
      real(c_double), intent(inout) :: MinV
    end subroutine
  end interface
  private :: PGSLib_Global_MINLOC_int_C, PGSLib_Global_MINLOC_float_C, PGSLib_Global_MINLOC_double_C

  !! Functions from par-src/reductions/redux-c.c
  interface PGSLib_Global_MAXLOC_C  ! SEE NOTE 3
    subroutine PGSLib_Global_MAXLOC_int_C (MaxV, GlobalIndex) &
        bind(c, name="pgslib_global_maxloc_int_c")
      import c_int
      integer(c_int), intent(inout) :: GlobalIndex
      integer(c_int), intent(inout) :: MaxV
    end subroutine
    subroutine PGSLib_Global_MAXLOC_float_C (MaxV, GlobalIndex) &
        bind(c, name="pgslib_global_maxloc_float_c")
      import c_int, c_float
      integer(c_int), intent(inout) :: GlobalIndex
      real(c_float),  intent(inout) :: MaxV
    end subroutine
    subroutine PGSLib_Global_MAXLOC_double_C (MaxV, GlobalIndex) &
        bind(c, name="pgslib_global_maxloc_double_c")
      import c_int, c_double
      integer(c_int), intent(inout) :: GlobalIndex
      real(c_double), intent(inout) :: MaxV
    end subroutine
  end interface
  private :: PGSLib_Global_MAXLOC_int_C, PGSLib_Global_MAXLOC_float_C, PGSLib_Global_MAXLOC_double_C

end module pgslib_c_binding
