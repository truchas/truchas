MODULE Scatter_Test_Module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE -
  !    Test the scatter routines provided by PGSLib.
  !    The routines in this module are designed to be called from a driver.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE PGSLib_MODULE
  USE Test_scatter_module
  implicit none

  PRIVATE
  PUBLIC :: Test_Scatter_SUM, Test_Scatter_MAX

CONTAINS
  subroutine Test_Scatter_SUM(elem_src_db, n_nodes, elements_nodes, elements_nodes_global, MASK)
    real (KIND(1.0D0)),                   &
         &                intent(IN   ),  &
         &                dimension(:,:) :: elem_src_db
    integer,                              &
         &                intent(INOUT ),  &
         &                dimension(:,:) :: elements_nodes, elements_nodes_global
    logical,                              &
         &                intent(IN   ),  &
         &                optional     ,  &
         &                dimension(:,:) :: Mask

    integer :: n_nodes
    !Local variables
    integer :: n_nodes_tot, n_elems, n_elems_tot, n_elem_1
    integer :: mem_error, n, e, error_count
    integer, dimension(1) :: error_index
    logical :: fatal, local_error
    character (LEN=1024) :: output_string
    real (KIND(elem_src_db)),   &
         &                dimension(n_nodes)           :: node_dest_db
    real (KIND(elem_src_db)), pointer, dimension(:,:) :: elem_src_db_tot
    integer,                  pointer, dimension(:,:) :: elements_nodes_tot
    real (KIND(elem_src_db)), pointer, dimension(  :) :: node_dest_db_tot, node_dest_db_tot_expected
    logical                 , pointer, dimension(:,:) :: Mask_Tot

    ! Use this to determine tolerances
    integer (KIND(elements_nodes)),   &
         &                dimension(n_nodes)           :: node_degree
    integer (KIND(elements_nodes)),   &
         &                dimension(SIZE(elements_nodes,1),  &
         &                          SIZE(elements_nodes,2)) :: elem_src_int
    real (KIND(elem_src_db)),   &
         &                dimension(n_nodes)           :: Tol_Local
    
    real (KIND(elem_src_db)),   POINTER, &
         &                dimension(:)           :: Tol_Tot
    

    ! Find mesh sizes
    n_elems   = size(elements_nodes,2)
    n_elem_1 = size(elements_nodes,1)
    n_nodes_tot = PGSLib_Global_SUM(n_nodes)
    n_elems_tot = PGSLib_Global_SUM(n_elems)
    
    ! Set error tolerance
    ! Set tolerance as a function of node degree
    elem_src_int = 1
    node_degree = 0
    if (PRESENT(MASK)) then
       call TEST_EN_Scatter_SUM(node_degree, elem_src_int, &
            &                   elements_nodes, MASK)
    else
       call TEST_EN_Scatter_SUM(node_degree, elem_src_int, &
            &                   elements_nodes)
    end if
    
    Tol_Local = 5.0*node_degree*EPSILON(node_dest_db)
    write(output_string,'("Maximum node degree = ",i8, "  Maximum Tolerance = ", e12.4)') &
         &               pgslib_global_maxval(node_degree),    &
         &               pgslib_global_maxval(Tol_local)
    call pgslib_output(output_string)

    if (PGSLib_Inquire_IO_P()) then
       allocate(elem_src_db_tot(n_elem_1, n_elems_tot), &
            &   elements_nodes_tot(n_elem_1, n_elems_tot))
       allocate(node_dest_db_tot(n_nodes_tot), &
            &   node_dest_db_tot_expected(n_nodes_tot) )
       allocate(Mask_Tot(n_elem_1, n_elems_tot))
       allocate(Tol_Tot(n_nodes_tot))
    else
       allocate(elem_src_db_tot(1, 1), &
            &   elements_nodes_tot(1, 1))
       allocate(node_dest_db_tot(1), &
            &   node_dest_db_tot_expected(1) )
       allocate(Mask_Tot(1,1))
       allocate(Tol_Tot(1))
    end if
    
       
    node_dest_db = 0.0

#ifdef DEBUG_TEST_GS
    if (present(MASK))then
       do e = 1, size(elements_nodes,2)
          do n = 1, size(elements_nodes,1)
             write(output_string,*) 'n,e, elements_nodes, mask, e_global: ', &
                  &    n,e, elements_nodes(n,e), mask(n,e), &
                  &    elements_nodes_global(n,e), elem_src_db(n,e)
             call pgslib_output(output_string)
          end do
       end do
    else
       do e = 1, size(elements_nodes,2)
          do n = 1, size(elements_nodes,1)
             write(output_string,*) 'n,e, elements_nodes, e_global: ', &
                  &    n,e, elements_nodes(n,e), &
                  &    elements_nodes_global(n,e), elem_src_db(n,e)
             call pgslib_output(output_string)
          end do
       end do
    end if
    call pgslib_flush_output()
#endif
    
    if (PRESENT(Mask)) then
       call TEST_EN_Scatter_SUM(node_dest_db, elem_src_db, elements_nodes, MASK)
    else
       call TEST_EN_Scatter_SUM(node_dest_db, elem_src_db, elements_nodes)
    end if

    do n = 1, n_elem_1
       call pgslib_collate(elem_src_db_tot(n,:), elem_src_db(n,:))
       call pgslib_collate(elements_nodes_tot(n,:), elements_nodes_global(n,:))
       if (present(MASK)) then
          call pgslib_collate(Mask_Tot(n,:), MASK(n,:))
       else
          Mask_Tot = .TRUE.
       endif
    end do
       
    node_dest_db_tot = 0.0
    call pgslib_collate(node_dest_db_tot, node_dest_db)
    Tol_Tot = 0.0
    call pgslib_collate(Tol_Tot, Tol_local)
       
    if (pgslib_inquire_IO_P()) then
       node_dest_db_tot_expected = 0.0
       do e = 1, size(elem_src_db_tot,2)
          do n = 1, size(elem_src_db_tot,1)
             if (.NOT. Mask_Tot(n,e)) cycle
             node_dest_db_tot_expected(elements_nodes_tot(n,e)) =                &
                  &         node_dest_db_tot_expected(elements_nodes_tot(n,e)) + &
                  &         elem_src_db_tot(n,e)
          end do
       end do
       
       local_error = ANY(ABS(node_dest_db_tot - node_dest_db_tot_expected) > Tol_Tot)
    else
       local_error = .false.
    end if

    local_error = pgslib_global_any(local_error)
  
    IF (local_error) then
       call pgslib_output("FAILED: Scatter test in test_scatter")
       call pgslib_error("FAILED: Scatter test in test_scatter")
       if (pgslib_inquire_IO_P()) then
          error_index = MAXLOC(ABS(node_dest_db_tot - node_dest_db_tot_expected))
          error_count = COUNT(ABS(node_dest_db_tot - node_dest_db_tot_expected) > Tol_Tot)
          write(output_string,'("Error_count = ",i8,"  Max error = ",e10.4," at ",i8)') &
               &                 error_count, &
               &  ABS(node_dest_db_tot(error_index(1)) - node_dest_db_tot_expected(error_index(1))), &
               &  error_index
          call pgslib_error(output_string)
          call pgslib_output(output_string)
       end if
    ELSE
       call pgslib_output(' REJOICE: Scatter passed a simple test.'   )
    end IF



     deallocate(node_dest_db_tot_expected, node_dest_db_tot,       &
          &     elements_nodes_tot, elem_src_db_tot)

     return
   end subroutine Test_Scatter_SUM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Test Scatter_Max
  subroutine Test_Scatter_Max(elem_src_db, n_nodes, elements_nodes, elements_nodes_global, MASK)
    real (KIND(1.0D0)),                   &
         &                intent(IN   ),  &
         &                dimension(:,:) :: elem_src_db
    integer,                              &
         &                intent(INOUT ),  &
         &                dimension(:,:) :: elements_nodes, elements_nodes_global

    logical,                              &
         &                intent(IN   ),  &
         &                optional     ,  &
         &                dimension(:,:) :: Mask

    integer :: n_nodes
    !Local variables
    integer :: n_nodes_tot, n_elems, n_elems_tot, n_elem_1
    integer :: mem_error, n, e, error_count
    integer, dimension(1) :: error_index
    logical :: fatal, local_error
    character (LEN=1024) :: output_string
    real (KIND(elem_src_db)) :: Tol

    real (KIND(elem_src_db)),   &
         &                dimension(n_nodes)           :: node_dest_db
    real (KIND(elem_src_db)), pointer, dimension(:,:) :: elem_src_db_tot
    integer,                  pointer, dimension(:,:) :: elements_nodes_tot
    real (KIND(elem_src_db)), pointer, dimension(  :) :: node_dest_db_tot, node_dest_db_tot_expected
    logical                 , pointer, dimension(:,:) :: Mask_Tot


    ! Find mesh sizes
    n_elems   = size(elements_nodes,2)
    n_elem_1 = size(elements_nodes,1)
    n_nodes_tot = PGSLib_Global_SUM(n_nodes)
    n_elems_tot = PGSLib_Global_SUM(n_elems)
    
    ! Set error tolerance
    Tol = n_elems*EPSILON(Tol)
    if (PGSLib_Inquire_IO_P()) then
       allocate(elem_src_db_tot(n_elem_1, n_elems_tot), &
            &   elements_nodes_tot(n_elem_1, n_elems_tot))
       allocate(node_dest_db_tot(n_nodes_tot), &
            &   node_dest_db_tot_expected(n_nodes_tot) )
       allocate(Mask_Tot(n_elem_1, n_elems_tot))
    else
       allocate(elem_src_db_tot(1, 1), &
            &   elements_nodes_tot(1, 1))
       allocate(node_dest_db_tot(1), &
            &   node_dest_db_tot_expected(1) )
       allocate(Mask_Tot(1,1))
    end if
    
       
    node_dest_db = -HUGE(node_dest_db)/2.

    if (PRESENT(Mask)) then
       call TEST_EN_Scatter_Max(node_dest_db, elem_src_db, elements_nodes, MASK)
    else
       call TEST_EN_Scatter_Max(node_dest_db, elem_src_db, elements_nodes)
    end if

    call pgslib_output("finished TEST_EN_Scatter_Max")
    call pgslib_flush_output()
    do n = 1, n_elem_1
       call pgslib_collate(elem_src_db_tot(n,:), elem_src_db(n,:))
       call pgslib_collate(elements_nodes_tot(n,:), elements_nodes_global(n,:))
       if (present(MASK)) then
          call pgslib_collate(Mask_Tot(n,:), MASK(n,:))
       else
          Mask_Tot = .TRUE.
       endif
    end do
       
    node_dest_db_tot = -HUGE(node_dest_db_tot)/2.

    call pgslib_flush_output()

    call pgslib_collate(node_dest_db_tot, node_dest_db)
       
    call pgslib_flush_output()

    if (pgslib_inquire_IO_P()) then
       node_dest_db_tot_expected = -HUGE(node_dest_db_tot_expected)/2.
       do e = 1, size(elem_src_db_tot,2)
          do n = 1, size(elem_src_db_tot,1)
             if (.NOT. Mask_Tot(n,e)) cycle
             node_dest_db_tot_expected(elements_nodes_tot(n,e)) =                    &
                  &         MAX(node_dest_db_tot_expected(elements_nodes_tot(n,e)),  &
                  &             elem_src_db_tot(n,e))
          end do
       end do
       
       local_error = ANY((ABS(node_dest_db_tot - node_dest_db_tot_expected) > Tol) .AND. &
            &                (node_dest_db_tot > -HUGE(node_dest_db_tot)/2.))
    else
       local_error = .false.
    end if

    local_error = pgslib_global_any(local_error)

    IF (local_error) then
       call pgslib_output("FAILED: Scatter MAX test in test_scatter")
       call pgslib_error("FAILED: Scatter MAX test in test_scatter")
       if (pgslib_inquire_IO_P()) then
          error_index = MAXLOC((ABS(node_dest_db_tot - node_dest_db_tot_expected)), &
            &                MASK=(node_dest_db_tot > -HUGE(node_dest_db_tot)/2.))
          error_count = COUNT((ABS(node_dest_db_tot - node_dest_db_tot_expected) > Tol).AND. &
            &                (node_dest_db_tot > -HUGE(node_dest_db_tot)/2.))
          write(output_string,'("Error_count = ",i8,"  Max error = ",e10.4," at ",i8)') &
               &                 error_count, &
               &  ABS(node_dest_db_tot(error_index(1)) - node_dest_db_tot_expected(error_index(1))), &
               &  error_index
          call pgslib_error(output_string)
          call pgslib_output(output_string)
       end if
    ELSE
       call pgslib_output(' REJOICE: Scatter MAX passed a test.'   )
    end IF



     deallocate(node_dest_db_tot_expected, node_dest_db_tot,       &
          &     elements_nodes_tot, elem_src_db_tot)

     return
   end subroutine Test_Scatter_MAX

 end MODULE Scatter_Test_Module
 
