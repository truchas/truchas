MODULE Gather_Test_Module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE -
  !    Test the gather routines provided by PGSLib.
  !    The routines in this module are designed to be called from a driver.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE PGSLib_MODULE
  USE Test_gather_module
  implicit none

  PRIVATE
  PUBLIC :: Test_Gather
  
CONTAINS
  subroutine Test_Gather(node_src_db, elements_nodes, elements_nodes_global)
    real (KIND(1.0D0)),                   &
         &                intent(IN   ),  &
         &                dimension(  :) :: node_src_db
    integer,                              &
         &                intent(INOUT ),  &
         &                dimension(:,:) :: elements_nodes, elements_nodes_global

    !Local variables
    integer :: n_nodes, n_nodes_tot, n_elems, n_elems_tot, n_elem_1
    integer :: mem_error, n, e, error_count
    logical :: fatal, local_error
    integer, dimension(2) :: error_index
    character (LEN=1024) :: output_string
    real (KIND(node_src_db)) :: Tol


    real (KIND(node_src_db)),   &
         &                dimension(size(elements_nodes,1), size(elements_nodes,2)) :: elem_dest_db
    real (KIND(node_src_db)), pointer, dimension(  :) :: node_src_db_tot
    integer,                  pointer, dimension(:,:) :: elements_nodes_tot
    real (KIND(node_src_db)), pointer, dimension(:,:) :: elem_dest_db_tot, elem_dest_db_tot_expected

    ! Find mesh sizes
    n_nodes   = size(node_src_db,1)
    n_elems   = size(elements_nodes,2)
    n_elem_1 = size(elements_nodes,1)
    n_nodes_tot = PGSLib_Global_SUM(n_nodes)
    n_elems_tot = PGSLib_Global_SUM(n_elems)
    
    ! Set error tolerance
    Tol = n_elems*EPSILON(Tol)
    ! Allocate arrays
    if (PGSLib_Inquire_IO_P()) then
       fatal = .false.
       allocate(node_src_db_tot(n_nodes_tot), STAT=mem_error)
       fatal = fatal .or. (mem_error /= 0)
       allocate(elements_nodes_tot(n_elem_1, n_elems_tot), STAT=mem_error)
       fatal = fatal .or. (mem_error /= 0)
       allocate(elem_dest_db_tot  (n_elem_1, n_elems_tot), STAT=mem_error)
       fatal = fatal .or. (mem_error /= 0)
       allocate(elem_dest_db_tot_expected(n_elem_1, n_elems_tot), STAT=mem_error)
       fatal = fatal .or. (mem_error /= 0)
     else
        ALLOCATE(node_src_db_tot(1))
        ALLOCATE(elements_nodes_tot(1,1))
        ALLOCATE(elem_dest_db_tot(1,1))
        ALLOCATE(elem_dest_db_tot_expected(1,1))
        fatal = .false.
     end if

     call pgslib_check_error(fatal, "Could not allocate memory for large arrays in test-gather")
     
     elem_dest_db = 0
     call Test_en_gather(elem_dest_db, node_src_db, elements_nodes)

#ifdef DEBUG_TEST_GS
    do e = 1, size(elements_nodes,2)
       do n = 1, size(elements_nodes,1)
          write(output_string,*) 'n,e, elements_nodes, e_global, e_dest: ', &
               &    n,e, elements_nodes(n,e), &
               &    elements_nodes_global(n,e), elem_dest_db(n,e)
          call pgslib_output(output_string)
       end do
    end do
    call pgslib_flush_output()
#endif

     do n = 1, n_elem_1
        call pgslib_collate(elem_dest_db_tot(n,:)  , elem_dest_db(n,:)         )
        call pgslib_collate(elements_nodes_tot(n,:), elements_nodes_global(n,:))
     end do
     call pgslib_collate(node_src_db_tot, node_src_db)

     fatal = .false.
     if (PGSLib_Inquire_IO_P()) then
        elem_dest_db_tot_expected = 0.0
        do e=1,n_elems_tot
           do n = 1,n_elem_1
              elem_dest_db_tot_expected(n,e) = node_src_db_tot(elements_nodes_tot(n,e))
           end do
        end do
        fatal = ANY(elem_dest_db_tot_expected /= elem_dest_db_tot)
     end if
     
     fatal = PGSLib_Global_ANY(fatal)
     local_error = PGSLib_Global_ANY(fatal)

    IF (local_error) then
       call pgslib_output("FAILED gather test in test-gather")
       call pgslib_error("FAILED gather test in test-gather")
       if (pgslib_inquire_IO_P()) then
          error_index = MAXLOC(ABS(elem_dest_db_tot - elem_dest_db_tot_expected))
          error_count = COUNT(ABS(elem_dest_db_tot - elem_dest_db_tot_expected) > Tol)
          write(output_string,'("Error_count = ",i8,"  Max error = ",e10.4," at ",i8,i8)') &
               &                 error_count, &
               &  ABS(elem_dest_db_tot(error_index(1),error_index(2)) -                 &
               &      elem_dest_db_tot_expected(error_index(1), error_index(2))), &
               &  error_index(1), error_index(2)
          call pgslib_error(output_string)
          call pgslib_output(output_string)
       end if
    ELSE
       call pgslib_output("REJOICE: PASSED simple gather test.")
    end IF

     deallocate(elem_dest_db_tot_expected, elem_dest_db_tot,       &
          &     elements_nodes_tot, node_src_db_tot)

     return
   end subroutine Test_Gather

 end MODULE Gather_Test_Module
 
