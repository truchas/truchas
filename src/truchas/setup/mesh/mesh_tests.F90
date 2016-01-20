!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Mesh_Tests


  !=======================================================================
  ! Purpose(s):
  !
  !   Routines to test consistency of mesh and mesh-related derived types.
  !
  !   Public Interface:
  !
  !     * call Test_All_Neighbors()
  !
  !        Check the Mesh%Ngbr_Cells_All component for consistency
  !	
  ! Contains: Test_All_Neighbors
  !
  ! Author(s): Robert Ferrell (ferrell.cpca.com)
  !
  !======================================================================

  use truchas_logging_services
  implicit none
  PRIVATE

  PUBLIC :: Test_All_Neighbors

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  function  Test_All_Neighbors() RESULT(Passed)
    !====================================================================
    ! Purpose(s):
    !
    !   Test the correctness of Mesh%Ngbr_Cells_All.  Cell neighbors
    !   share at least one node.
    !
    !   Test is done in two ways.
    !   1. Check that all identified neighbors share at least one
    !      node.  (All identified neighbors are actual neighbors.)
    !   2. Check that all cells which share a node are identified
    !      as neighbors. (Found all actual neighbors.)
    !
    !   If either condition fails, then the test fails.
    !   This routine returns TRUE if the test passes, otherwise
    !   returns FALSE.
    !
    !   The same value is returned on all processors.
    !
    !   The tests are performend iff debug is turned on a run-time
    !   (run with -d).  Otherwise the tests are not done,
    !   and the return value is TRUE.
    !====================================================================
    use debug_control_data
    use gs_module,            only: EE_Gather
    use mesh_module,          only: Mesh
    use mesh_parameter_module, only: ncells, nvc
    use pgslib_module,        only: PGSLib_Global_All
    use var_vector_module

    ! Arguments and return values
    logical :: Passed

    ! Local variables
    integer, dimension(ncells)  :: donor_vrtx
    type (int_var_vector), dimension(ncells) :: gathered_vrtx
    type (log_var_vector), dimension(ncells) :: confirmed_ngbrs
    integer, POINTER, dimension(:) :: ngbr_vrtx
    integer, POINTER, dimension(:) :: ngbr_cells_orig
    logical, POINTER, dimension(:) :: ngbr_check

    integer :: donor_v, owner_c, owner_v, c, last_false, last_false_cell, n
    integer :: given_size, found_size, bad_neighbor
    integer :: out_stat

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    out_stat = 0
    ! Check if debugging is on, otherwise return.
    if (debug == DEBUG_NONE) then
       Passed = .TRUE.
       RETURN
    else
       if (debug >= DEBUG_NOISY) then
          call TLS_info ('')
          call TLS_info (' Testing all-neighbor connectivity ... ', advance=.false.)
       end if
    end if
    

    ! ***** Test Part 1, all neighbors are actual neighbors? *****

    ! Need a var_vector array to hold the result of the gather
    Call CREATE(ARRAY = gathered_vrtx,             &
                SIZES = SIZES(Mesh%Ngbr_Cells_All) )

    ! Need a var_vector array to hold the result of the tests
    Call CREATE(ARRAY = confirmed_ngbrs,           &
                SIZES = SIZES(Mesh%Ngbr_Cells_All) )

    ! Default is that we didn't confirm this as a neighbor.
    ngbr_check => FLATTEN(confirmed_ngbrs)
    ngbr_check = .FALSE.

    ! Loop over all donor vertices
    do donor_v = 1,nvc
      ! Pluck out the donor_v vertex, note that need original vertex numbers
      donor_vrtx = Mesh%Ngbr_Vrtx_Orig(donor_v)

      ! Gather into the neighbor list
      call EE_Gather(DEST   = gathered_vrtx, &
                     SOURCE = donor_vrtx     )

      ! Loop over all cells to check that for each listed neighbor,
      ! the gathered vertex is the same as one of the owner cell's vertices.
      do owner_c = 1, ncells
         ngbr_check => FLATTEN(confirmed_ngbrs(owner_c))
         ngbr_vrtx  => FLATTEN(gathered_vrtx  (owner_c))
         do owner_v = 1, nvc
            ngbr_check = ngbr_check .OR.               &
                         ( Mesh(owner_c)%Ngbr_Vrtx_Orig(owner_v) == ngbr_vrtx )
         end do
      end do
      
    end do
      
      
    ! Now test to see that all identified neighbors got flagged as 
    ! an actual neighbor.


    Passed = .TRUE.
    CELL_LOOP: do c = 1, ncells
       ngbr_check => FLATTEN(confirmed_ngbrs(c))
       ngbr_cells_orig => FLATTEN(Mesh(c)%Ngbr_Cells_ALL_Orig)
       NGBR: do n = 1, SIZE(ngbr_check)
          Passed = Passed .AND. ngbr_check(n)
          if (.NOT. ngbr_check(n) ) then
             ! Put these two lines in because totalview was giving alarming (false) information
             given_size = SIZES(confirmed_ngbrs(c))
             found_size = SIZE(ngbr_check)

             last_false = n
             last_false_cell  = c
             bad_neighbor = ngbr_cells_orig(n)
             

          end if
       end do NGBR
     end do CELL_LOOP


    ! ***** Finished Test Part 1, all neighbors are actual neighbors? *****



     ! Return same value on all processors
     Passed = PGSLib_Global_ALL(Passed)

     if (debug >= DEBUG_NOISY) then
        call TLS_info ('done.')
     end if

  end function Test_All_Neighbors

end MODULE Mesh_tests
