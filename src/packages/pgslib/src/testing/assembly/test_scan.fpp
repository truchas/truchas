!!FPP!! This file is to be included by test_scan.F

!!FPP!! _DATA_TYPE_ and _RAND_OP_ must be defined before including this file


    ! $Id: test_scan.fpp,v 1.1.1.1 2000/10/11 22:44:26 ferrell Exp $

#define _MESSAGE_STRING_(PF,D,OS,S) PF//D//OS//S//" test"


    implicit none
    ! Subroutine arguments
    logical :: fatal_error, warning_error

    ! Arrays sizes, local size for distributed arrays, total size for total arrays
    integer :: N_Loc, N_Tot, Start, Stop, Step
    ! All arrays come in two sizes, distributed and total
    _DATA_TYPE_, pointer, dimension(:) :: Src, Dest
    _DATA_TYPE_, pointer, dimension(:) :: Src_Tot, Dest_Tot, Dest_Expected_Tot
    _DATA_TYPE_, pointer, dimension(:) :: Dest_Error_Tot
    REAL(PGSLib_Double_Type) :: TOLERANCE
    logical, pointer, dimension(:) :: Seg, TMask
    logical, pointer, dimension(:) :: Seg_Tot, TMask_Tot

    ! Misc scalars
    integer, dimension(1) :: error_index
    integer :: K1, K2, MaxSize, error_count
    integer :: thisPE, nPE, I
    logical :: error_flag
    character (LEN=1024) :: output_string


    ! These are used for psuedo random number
    k1 = 1000003
    k2 = 262147
    MaxSize = 1101
    thisPE = PGSLib_Inquire_thisPE()
    nPE    = PGSLib_Inquire_nPE()
    TOLERANCE = MaxSize*EPSILON(1.0)

    ! Array size is random
    N_Loc = MOD(thispe*MOD(K1,MaxSize) + K2, MaxSize)
    ! Test 0 sized arrays
    if (thisPE == 2) N_Loc = 0
    if (thisPE == 3) N_Loc = 0
    N_Tot = PGSLib_Global_SUM(N_Loc)
    write(output_string,'( "N_Loc = ", i8, " N_Tot = ",i8, " Tol = ",e12.4)') &
         N_Loc, N_Tot, TOLERANCE
    call pgslib_output(output_string)
    ! Allocate arrays
    ALLOCATE(  Src(N_Loc),                    &
         &    Dest(N_Loc),                    &
         &     Seg(N_Loc),                    &
         &   TMask(N_Loc))
    if (PGSLib_Inquire_IO_P()) then
       ALLOCATE(           Src_Tot(N_Tot),                    &
            &             Dest_Tot(N_Tot),                    &
            &    Dest_Expected_Tot(N_Tot),                    &
            &       Dest_Error_Tot(N_Tot),                    &
            &              Seg_Tot(N_Tot),                    &
            &            TMask_Tot(N_Tot))
    else
       ALLOCATE(           Src_Tot(1),                    &
            &             Dest_Tot(1),                    &
            &    Dest_Expected_Tot(1),                    &
            &       Dest_Error_Tot(1),                    &
            &              Seg_Tot(1),                    &
            &            TMask_Tot(1))
    end if

    ! Initialize arrays
    Src = (/ (I, I=1,N_Loc) /)
    Src = (MOD((thisPE*INT(Src))*K1 + K2, MaxSize))
    Src = _RAND_OP_(Src)
    Src = ABS(Src)

    call PGSLib_Collate(Src_Tot, Src)

    ! Initialize loop parameters
    Start = _START_
    Stop  = _STOP_
    Step  = _STEP_

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! No Segment, No Mask
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Dest = _PGSLib_SCAN_ROUTINE_(Src)

    ! Test the result
    call PGSLib_Collate(Dest_Tot, Dest)

    Dest_Expected_Tot = 0
    Dest_Expected_Tot(Start) = Src_Tot(Start)
    DO I = Start + Step, Stop, Step
       Dest_Expected_Tot(I) = Dest_Expected_Tot(I - Step) + Src_Tot(I)
    END DO


    IF(PGSLib_Inquire_IO_P()) then
       DO I = 1, SIZE(Dest_Expected_Tot,1)
          if (Dest_Expected_Tot(I) > TOLERANCE) then
             dest_error_tot(I) = ABS(Dest_Expected_Tot(I) - Dest_Tot(I))/ABS(Dest_Expected_Tot(I))
          else
             dest_error_tot(I) = ABS(Dest_Expected_Tot(I) - Dest_Tot(I))
          end if

       end DO
       error_flag = ANY(dest_error_tot > TOLERANCE)
    else
       error_flag = .FALSE.
    end IF

    error_flag = PGSLib_Global_ANY(Error_flag)
    if (error_flag) then
       output_string= _MESSAGE_STRING_("Failed", _DATA_TYPE_STRING_ , ", no segment, no mask ",_SCAN_TEST_STRING_)
       call pgslib_output(output_string)
       call pgslib_error(output_string)
       if (pgslib_inquire_IO_P()) then
          error_index = MAXLOC(dest_error_tot, MASK=(Dest_Tot > -HUGE(Dest_Tot)/2.))
          error_count = COUNT((dest_error_tot > TOLERANCE) .AND. &
            &                (Dest_Tot > -HUGE(Dest_Tot)/2.))
          write(output_string,'("Error_count = ",i8,"  Max error = ",e10.4," at ",i8)') &
               &                 error_count, &
               &  REAL(Dest_Error_Tot(error_index(1))), error_index
          call pgslib_error(output_string)
          call pgslib_output(output_string)
#ifdef DEBUG_ASSEMBLY
          open(UNIT=20, FILE='Scan-out', STATUS='UNKNOWN')
          WRITE(20, '("Data For ",a,a," test")') _DATA_TYPE_STRING_, &
               &                                 _SCAN_TEST_STRING_
          DO I = 1, SIZE(Dest_Expected_Tot,1)
             IF (Dest_Error_Tot(I) .GT. TOLERANCE) then
                WRITE(20, *) I, Src_Tot(I), Dest_Expected_Tot(I), Dest_Tot(I), &
                     &  Dest_error_Tot(I), "DIFF"
             else
                WRITE(20, *) I, Src_Tot(I), Dest_Expected_Tot(I), Dest_Tot(I), &
                     &  Dest_error_Tot(I), "SAME"
             END IF
          END DO
          close(20)
#endif
       END IF

    end if

    warning_error = warning_error .or. error_flag
    fatal_error   = fatal_error   .or. error_flag

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Segment, No Mask
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Setup Segment.  Pick segments of mean length 20
    Seg = PGSLib_PARITY_PREFIX(MOD(INT(MaxSize*Src),MaxSize/20) == 0)
    call pgslib_collate(Seg_tot, Seg)

    Dest = _PGSLib_SCAN_ROUTINE_(Src, SEGMENT=Seg)

    ! Test the result
    call PGSLib_Collate(Dest_Tot, Dest)

    Dest_Expected_Tot = 0
    Dest_Expected_Tot(Start) = Src_Tot(Start)
    DO I = Start + Step, Stop, Step
       if(Seg_Tot(i) .EQV. Seg_Tot(i - Step)) then
          Dest_Expected_Tot(I) = Dest_Expected_Tot(I - Step) + Src_Tot(I)
       else
          Dest_Expected_Tot(I) = Src_Tot(I)
       end if
    END DO

    IF(PGSLib_Inquire_IO_P()) then
       error_flag = ANY(ABS(Dest_Expected_Tot - Dest_Tot) .GT. TOLERANCE)
    else
       error_flag = .FALSE.
    end IF

    error_flag = PGSLib_Global_ANY(Error_flag)
    if (error_flag) then
       output_string = _MESSAGE_STRING_("Failed ", _DATA_TYPE_STRING_ , " segment, no mask ", _SCAN_TEST_STRING_)
       call pgslib_output(output_string)
       call pgslib_error(output_string)
    end if

    warning_error = warning_error .or. error_flag
    fatal_error   = fatal_error   .or. error_flag


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! No Segment, Mask
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Setup Mask.  Pick 20:1 true/false ratio

    TMASK = MOD(INT(MaxSize*Src),MaxSize/20) == 0
    call pgslib_collate(TMASK_tot, TMASK)

    Dest = _PGSLib_SCAN_ROUTINE_(Src, MASK=TMask)

    ! Test the result
    call PGSLib_Collate(Dest_Tot, Dest)

    Dest_Expected_Tot = 0
    Dest_Expected_Tot(Start) = MERGE(Src_Tot(Start),_OP_ID_,TMask_Tot(Start))
    DO I = Start + Step, Stop, Step
       Dest_Expected_Tot(I) = Dest_Expected_Tot(I - Step) + MERGE(Src_Tot(I),_OP_ID_,TMask_Tot(I))
    END DO

    IF(PGSLib_Inquire_IO_P()) then
       DO I = 1, SIZE(Dest_Expected_Tot,1)
          if (Dest_Expected_Tot(I) > TOLERANCE) then
             dest_error_tot(I) = ABS(Dest_Expected_Tot(I) - Dest_Tot(I))/ABS(Dest_Expected_Tot(I))
          else
             dest_error_tot(I) = ABS(Dest_Expected_Tot(I) - Dest_Tot(I))
          end if
       end DO
       error_flag = ANY(dest_error_tot > TOLERANCE)
    else
       error_flag = .FALSE.
    end IF

    error_flag = PGSLib_Global_ANY(Error_flag)
    if (error_flag) then
       output_string= _MESSAGE_STRING_("Failed", _DATA_TYPE_STRING_ , ", no segment, mask ",_SCAN_TEST_STRING_)
       call pgslib_output(output_string)
       call pgslib_error(output_string)
       if (pgslib_inquire_IO_P()) then
          error_index = MAXLOC(dest_error_tot, MASK=(Dest_Tot > -HUGE(Dest_Tot)/2.))
          error_count = COUNT((dest_error_tot > TOLERANCE) .AND. &
            &                (Dest_Tot > -HUGE(Dest_Tot)/2.))
          write(output_string,'("Error_count = ",i8,"  Max error = ",e10.4," at ",i8)') &
               &                 error_count, &
               &  REAL(Dest_Error_Tot(error_index(1))), error_index
          call pgslib_error(output_string)
          call pgslib_output(output_string)
#ifdef DEBUG_ASSEMBLY
          open(UNIT=20, FILE='Scan-out', STATUS='UNKNOWN')
          WRITE(20, '("Data For ",a,a," test")') _DATA_TYPE_STRING_, &
               &                                 _SCAN_TEST_STRING_
          DO I = 1, SIZE(Dest_Expected_Tot,1)
             IF (Dest_Error_Tot(I) .GT. TOLERANCE) then
                WRITE(20, *) I, Src_Tot(I), Dest_Expected_Tot(I), Dest_Tot(I), &
                     &  Dest_error_Tot(I), "DIFF"
             else
                WRITE(20, *) I, Src_Tot(I), Dest_Expected_Tot(I), Dest_Tot(I), &
                     &  Dest_error_Tot(I), "SAME"
             END IF
          END DO
          close(20)
#endif
       END IF
    end if
    warning_error = warning_error .or. error_flag
    fatal_error   = fatal_error   .or. error_flag


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Segment and Mask
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Setup Mask.  Pick 20:1 true/false ratio

    TMASK = MOD(INT(MaxSize*Src),MaxSize/20) == 0
    call pgslib_collate(TMASK_tot, TMASK)

    ! Setup Segment.  Pick segments of mean length 50
    Seg = PGSLib_PARITY_PREFIX(MOD(INT(MaxSize*Src),MaxSize/20) == 0)
    call pgslib_collate(Seg_tot, Seg)

    Dest = _PGSLib_SCAN_ROUTINE_(Src, SEGMENT=Seg, MASK=TMask)

    ! Test the result
    call PGSLib_Collate(Dest_Tot, Dest)

    Dest_Expected_Tot = 0
    Dest_Expected_Tot(Start) = MERGE(Src_Tot(Start),_OP_ID_,TMask_Tot(Start))
    DO I = Start + Step, Stop, Step
       if(Seg_Tot(i) .EQV. Seg_Tot(i - Step)) then
          Dest_Expected_Tot(I) = Dest_Expected_Tot(I - Step) + MERGE(Src_Tot(I),_OP_ID_,TMask_Tot(I))
       else
          Dest_Expected_Tot(I) = MERGE(Src_Tot(I),_OP_ID_,TMask_Tot(I))
       end if
    END DO

    IF(PGSLib_Inquire_IO_P()) then
       error_flag = ANY(ABS(Dest_Expected_Tot - Dest_Tot) .GT. TOLERANCE)
    else
       error_flag = .FALSE.
    end IF

    error_flag = PGSLib_Global_ANY(Error_flag)
    if (error_flag) then
       output_string = _MESSAGE_STRING_("Failed ", _DATA_TYPE_STRING_ , " segment, mask ", _SCAN_TEST_STRING_)
       call pgslib_output(output_string)
       call pgslib_error(output_string)
    end if

    warning_error = warning_error .or. error_flag
    fatal_error   = fatal_error   .or. error_flag




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finished all tests
    !

    if (.not. (warning_error .or. fatal_error)) then
       output_string = _MESSAGE_STRING_("Passed", _DATA_TYPE_STRING_ , " all segment/mask combinations ", _SCAN_TEST_STRING_)
       call pgslib_output(output_string)
    endif

    DEALLOCATE(Src,          &
         &          Dest,    &
         &          Seg,&
         &          TMask,&
         &          Src_Tot,&
         &          Dest_Tot,&
         &          Dest_Expected_Tot,&
         &          Dest_Error_Tot,&
         &          Seg_Tot,&
         &          TMask_Tot)


    return

#undef _RAND_OP_
#undef _DATA_TYPE_STRING_
#undef _OP_ID_
#undef _DATA_TYPE_
#undef _START_
#undef _STOP_
#undef _STEP_
#undef _PGSLib_SCAN_ROUTINE_
#undef _SCAN_TEST_STRING_
