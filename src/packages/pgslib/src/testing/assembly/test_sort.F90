MODULE TEST_SORT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE -
  !   Test the sorting/ranking routines
  !
  ! INPUT
  !   None
  ! OUTPUT
  !   Results of tests, warning_error and fatal_error
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  ! $Id: test_sort.F,v 1.3 2001/04/20 19:36:15 ferrell Exp $

  USE PGSLib_MODULE
  implicit none
  
  PRIVATE
  PUBLIC  :: Test_Sorting

CONTAINS
  subroutine test_sorting(local_error)
    implicit none
    logical :: local_error

    logical :: fatal_error, warning_error

    ! Misc scalars
    integer :: K1, K2, MaxSize, N_Loc
    integer :: thisPE, nPE

    ! Initialize error reporting variables
    local_error = .false.

    ! These are used for psuedo random number
    k1 = 1000003
    k2 = 262147
    MaxSize = 10111
    thisPE = PGSLib_Inquire_thisPE()
    nPE    = PGSLib_Inquire_nPE()

    ! First test, Array size is random
    N_Loc = MOD(thispe*K1 + K2, MaxSize)

    fatal_error = .false.
    warning_error = .false.
    call test_local_rank_int_1D(fatal_error, warning_error, N_LOCAL=N_Loc)
    local_error = local_error .or. (fatal_error .or. warning_error)
    
    ! Now test with some arrays having 0 size
    if (MOD(thisPE,2) /= 0) N_Loc = 0

    fatal_error = .false.
    warning_error = .false.
    call test_local_rank_int_1D(fatal_error, warning_error, N_LOCAL=N_Loc)
    local_error = local_error .or. (fatal_error .or. warning_error)

    ! Finally, test with all arrays having 0 size
    N_Loc = 0
    fatal_error = .false.
    warning_error = .false.
    call test_local_rank_int_1D(fatal_error, warning_error, N_LOCAL=N_Loc)
    local_error = local_error .or. (fatal_error .or. warning_error)

    return
  end subroutine test_sorting

  subroutine test_local_rank_int_1D(fatal_error, warning_error, N_LOCAL)
#define _DATA_TYPE_ integer
#define _OP_ID_ 0
#define _DATA_TYPE_STRING_ 'integer'
#define _RAND_OP_   
#define _START_ 1
#define _STOP_ SIZE(Dest_Expected_Tot, 1)
#define _STEP_  1
#define _PGSLib_SCAN_ROUTINE_ PGSLib_SUM_PREFIX
#define _RANK_TEST_STRING_    'PGSLib_GRADE_UP'

    implicit none
    ! Subroutine arguments
    logical :: fatal_error, warning_error
    integer, intent(IN) :: N_LOCAL

    ! Arrays sizes, local size for distributed arrays, total size for total arrays
    integer :: N_Loc, N_Tot, Start, Stop, Step
    ! All arrays come in two sizes, distributed and total
    _DATA_TYPE_, pointer, dimension(:) :: Src, Dest
    _DATA_TYPE_, pointer, dimension(:) :: Src_Tot, Dest_Tot, Dest_Expected_Tot
    _DATA_TYPE_, pointer, dimension(:) :: Sort_Index_Tot, Rank_Tot, Rank_Expected_Tot
    _DATA_TYPE_, pointer, dimension(:) :: Dest_Error
    integer, pointer, dimension(:)     :: Rank
    integer, pointer, dimension(:)     :: Temp_Array
    _DATA_TYPE_ :: TOLERANCE, Temp
    logical, pointer, dimension(:) :: Seg, TMask
    logical, pointer, dimension(:) :: Seg_Tot, TMask_Tot

    ! Misc scalars
    integer, dimension(1) :: error_index
    integer :: error_count, Temp_Index
    integer :: thisPE, nPE, I, J, K1, K2, MaxSize
    logical :: error_flag
    character (LEN=256):: out_string
    character (LEN=512):: output_string


    ! These are used for psuedo random number
    k1 = 1000003
    k2 = 262147
    MaxSize = 10111

    thisPE = PGSLib_Inquire_thisPE()
    nPE    = PGSLib_Inquire_nPE()
    TOLERANCE = 100*EPSILON(1.0)

    ! Array size is input
    N_Loc = N_LOCAL
    N_Tot = PGSLib_Global_SUM(N_Loc)

    write(output_string,*) 'N_Loc, N_Tot =', N_Loc, N_Tot
    call pgslib_output(output_string)
    ! Allocate arrays
    ALLOCATE(  Src(N_Loc),                    &
         &    Dest(N_Loc),                    &
         &    Rank(N_Loc),                    &
         &     Seg(N_Loc),                    &
         &   TMask(N_Loc),                    &
         &   Temp_Array(N_Loc))
    if (PGSLib_Inquire_IO_P()) then
       ALLOCATE(           Src_Tot(N_Tot),                    &
            &             Dest_Tot(N_Tot),                    &
            &    Dest_Expected_Tot(N_Tot),                    &
            &           Dest_Error(N_Tot),                    &
            &             Rank_Tot(N_Tot),                    &
            &    Rank_Expected_Tot(N_Tot),                    &
            &       Sort_Index_Tot(N_Tot),                    &
            &              Seg_Tot(N_Tot),                    &
            &            TMask_Tot(N_Tot))
    else
       ALLOCATE(           Src_Tot(1),                    &
            &             Dest_Tot(1),                    &
            &    Dest_Expected_Tot(1),                    &
            &           Dest_Error(1),                    &
            &             Rank_Tot(1),                    &
            &    Rank_Expected_Tot(1),                    &
            &       Sort_Index_Tot(1),                    &
            &              Seg_Tot(1),                    &
            &            TMask_Tot(1))
    end if

    ! Initialize arrays
    Src = (/ (I, I=1,N_Loc) /)
    Src = (MOD((thisPE*INT(Src))*K1 + K2, MaxSize))
    Src = _RAND_OP_(Src)

    call PGSLib_Collate(Src_Tot, Src)


    ! Get the ranking
    Rank = PGSlib_Grade_Up(Src)
    WRITE(out_string,*) 'MINVAL(Rank), MAXVAL(Rank) =', MINVAL(Rank), MAXVAL(Rank)
    call pgslib_output(out_string)
    ! Sort the data
    Call PGSLib_Permute(Dest, Src, Rank)
    ! Rank again just for grins
!    Rank = PGSlib_Grade_Up(Dest)

    ! We will need a simple serial sort for comparison
    Dest_Expected_Tot = Src_Tot
    Sort_Index_Tot = (/ (I, I=1,SIZE(Sort_Index_Tot,1) ) /)
    DO J = 2, SIZE(Dest_Expected_Tot,1)
       Temp = Dest_Expected_Tot(J)
       Temp_Index = Sort_Index_Tot(J)
       DO I = J-1, 1, -1
          IF (Dest_Expected_Tot(I) <= Temp) Exit
          Dest_Expected_Tot(I+1) = Dest_Expected_Tot(I)
          Sort_Index_Tot(I+1)    = Sort_Index_Tot(I)
       END DO
       Dest_Expected_Tot(I+1) = Temp
       Sort_Index_Tot(I+1)    = Temp_Index
    END DO

    DO I = 1, SIZE(Sort_Index_Tot,1)
       Rank_Expected_Tot(Sort_Index_Tot(I)) = I
    ENDDO

    ! Compare
    call pgslib_collate(Dest_Tot, Dest)

    call pgslib_collate(Rank_Tot, Rank)

#ifdef DEBUG_ASSEMBLY
    if (PGSLib_Inquire_IO_P()) then
       open(UNIT=20, FILE='Dest-out')
       open(UNIT=21, FILE='Rank-out')
       DO I = 1, SIZE(Dest_Tot,1)
          if (Dest_Tot(I) == Dest_Expected_Tot(I)) then
             WRITE(20,*) I, Dest_Tot(I), Dest_Expected_Tot(I), '  SAME'
          else
             WRITE(20,*) I, Dest_Tot(I), Dest_Expected_Tot(I), '  Diff = ',Dest_Tot(I)-Dest_Expected_Tot(I)
20           FORMAT(1x,i8, i8, i8, a, i8)
          end if
          if (Rank_Tot(I) == Rank_Expected_Tot(I)) then
             WRITE(21,21) I, Src_Tot(I), Rank_Tot(I), Rank_Expected_Tot(I), '  SAME'
          else
             WRITE(21,21) I, Src_Tot(I), Rank_Tot(I), Rank_Expected_Tot(I), '  Diff = ',Rank_Tot(I)-Rank_Expected_Tot(I)
          end if
21        FORMAT(1x,i8, i8, i8, i8, a, i8)
       END DO
       close(20)
       close(21)
    ENDIF
#endif

    IF (PGSLib_Inquire_IO_P()) then
       Dest_Error = Dest_Expected_Tot - Dest_Tot
       error_flag = ANY(ABS(Dest_Expected_Tot - Dest_Tot) .GT. TOLERANCE)
    else
       error_flag = .FALSE.
    end IF

    error_flag = PGSLib_Global_ANY(Error_flag)
    if (error_flag) then
       call pgslib_output('Failed '// _DATA_TYPE_STRING_ // ' no segment, no mask ' // _RANK_TEST_STRING_ // ' test')
       call pgslib_error ('Failed '// _DATA_TYPE_STRING_ // ' no segment, no mask ' // _RANK_TEST_STRING_ // ' test')
       if (pgslib_inquire_IO_P()) then
          WRITE(output_string,*) 'Comparing Destination arrays'
          call pgslib_error(output_string)
          call pgslib_output(output_string)
          error_index = MAXLOC((ABS(Dest_Tot - Dest_Expected_Tot)), &
            &                MASK=(Dest_Tot > -HUGE(Dest_Tot)/2.))
          error_count = COUNT((ABS(Dest_Tot - Dest_Expected_Tot) > Tolerance).AND. &
            &                (Dest_Tot > -HUGE(Dest_Tot)/2.))
          write(output_string,'("Error_count = ",i8,"  Max error = ",i8," at ",i8)') &
               &                 error_count, &
               &  ABS(Dest_Tot(error_index(1)) - Dest_Expected_Tot(error_index(1))), &
               &  error_index
          call pgslib_error(output_string)
          call pgslib_output(output_string)

          WRITE(output_string,*) 'Comparing Rank arrays'
          call pgslib_error(output_string)
          call pgslib_output(output_string)
          error_index = MAXLOC((ABS(Rank_Tot - Rank_Expected_Tot)), &
            &                MASK=(Rank_Tot > -HUGE(Rank_Tot)/2.))
          error_count = COUNT((ABS(Rank_Tot - Rank_Expected_Tot) > Tolerance).AND. &
            &                (Rank_Tot > -HUGE(Rank_Tot)/2.))
          write(output_string,'("Error_count = ",i8,"  Max error = ",i8," at ",i8)') &
               &                 error_count, &
               &  ABS(Rank_Tot(error_index(1)) - Rank_Expected_Tot(error_index(1))), &
               &  error_index
          call pgslib_error(output_string)
          call pgslib_output(output_string)
       end if
    end if

    warning_error = warning_error .or. error_flag
    fatal_error   = fatal_error   .or. error_flag

    if (.not. (warning_error .or. fatal_error)) then
       call pgslib_output('Passed' // _DATA_TYPE_STRING_ //  _RANK_TEST_STRING_ // ' test')
    endif

    ! Now try a segmented sort.
    ! Set the segment to something random.
    ! Initialize arrays
    Src = (/ (I, I=1,N_Loc) /)
    Src = (MOD((thisPE*INT(Src))*K1 + K2, MaxSize))
    Src = _RAND_OP_(Src)

    ! Need a cutoff to pick T & F
    Temp = (PGSLib_Global_MAXVAL(Src) - PGSLib_Global_MINVAL(Src) ) / 3
    Seg = PGSLib_PARITY_PREFIX( (Src <= Temp) )
    
    ! Now call the segmented sort.
    Rank = PGSlib_Grade_Up(Src, Seg)

    ! Sort the data
    Call PGSLib_Permute(Dest, Src, Rank)

    ! Clean up memory
    DEALLOCATE(Src, Dest, Rank, Seg, TMask, Temp_Array)
    DEALLOCATE(Src_Tot, Dest_Tot, Dest_Expected_Tot, &
         &     Sort_Index_Tot, Rank_Tot, Rank_Expected_Tot, &
         &     Seg_Tot,  TMask_Tot)

    return
  end subroutine test_local_rank_int_1D
end MODULE test_sort
