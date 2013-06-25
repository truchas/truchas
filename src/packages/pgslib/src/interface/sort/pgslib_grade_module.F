MODULE PGSLIB_GRADE_MODULE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Provide PGSLIB_GRADE_UP ranking routines
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! $Id: pgslib_grade_module.F,v 1.3 2001/04/20 19:36:15 ferrell Exp $

  use pgslib_type_module
  use pgslib_gs_module
  use pgslib_io_module
  use pgslib_misc_module
  use pgslib_reductions_module
  use pgslib_scan_module
  use pgslib_utility_module

  implicit none
  PRIVATE
  PUBLIC :: PGSLIB_GRADE_UP, PGSLib_GRADE_UP_Local

  INTERFACE PGSLIB_GRADE_UP
     MODULE PROCEDURE GRADE_UP_Int_1D
  END INTERFACE

  INTERFACE Determine_Buckets
     MODULE PROCEDURE Determine_Buckets_Int_1D
  END INTERFACE

  INTERFACE BucketNumber
     MODULE PROCEDURE BucketNumber_Int
  END INTERFACE

  INTERFACE PGSLib_GRADE_UP_Local
     MODULE PROCEDURE Rank_Up_Int_1D
  END INTERFACE

  INTERFACE Heap_Insert
     MODULE PROCEDURE Heap_Insert_Int
  END INTERFACE
  

  ! The data types, used by the C routines.
  integer, parameter :: FLOAT_DATA = 1
  integer, parameter :: INTEGER_DATA = 2
  ! The word sizes, in bytes
  integer, parameter :: INTEGER_SIZE = 4
  integer, parameter :: FLOAT_SIZE   = 4

  ! How much sampling to do
  integer, parameter :: SAMPLE_MAX = 20

  ! Type used for the BUCKETs
  TYPE INT_BUCKET_BNDRY
     integer (PGSLib_Int_Type) :: data
     integer (PGSLib_Int_Type) :: SegNum
  END TYPE INT_BUCKET_BNDRY

  TYPE FLOAT_BUCKET_BNDRY
     real    (PGSLib_Real_Type) :: data
     integer (PGSLib_Int_Type) :: SegNum
  END TYPE FLOAT_BUCKET_BNDRY

  TYPE FLOAT_BUCKET_TYPE
     type (FLOAT_BUCKET_BNDRY) :: Low, High
     integer (PGSLib_Int_Type) :: Bucket
  end TYPE FLOAT_BUCKET_TYPE

  TYPE INT_BUCKET_TYPE
     type (INT_BUCKET_BNDRY)   :: Low, High
     integer (PGSLib_Int_Type) :: Bucket
  end TYPE INT_BUCKET_TYPE

  type WorkNode_Int
     integer :: Index, Pos, SegNum
     integer :: data
  end type WorkNode_Int

  

  INTERFACE OPERATOR (.GE.)
     MODULE PROCEDURE GE_Bucket_Int
     MODULE PROCEDURE GE_Nodes_Int
  END INTERFACE
  

  INTERFACE OPERATOR (.LE.)
     MODULE PROCEDURE LE_Bucket_Int
  END INTERFACE
     
  INTERFACE OPERATOR (.LT.)
     MODULE PROCEDURE LT_Bucket_Int
  END INTERFACE
     
  INTERFACE OPERATOR (.EQ.)
     MODULE PROCEDURE EQ_Bucket_Int
  END INTERFACE
     



CONTAINS
  function GRADE_UP_Int_1D(KEY_ARRAY, Segment, DIM)
    implicit none

    integer, intent(IN   ), &
         &   dimension(                :) :: KEY_ARRAY
    integer, intent(IN   ), &
         &   optional                     :: DIM

    logical, intent(IN   ), &
         &   optional     , &
         &   dimension(SIZE(KEY_ARRAY,1)) :: Segment

    ! Function type
    integer,   dimension(SIZE(KEY_ARRAY,1)) :: GRADE_UP_Int_1D

    
    ! Local Variables
    type (INT_BUCKET_TYPE), &
         &   pointer,       &
         &   dimension(                :) :: BUCKETS
    integer (PGSLib_Int_Type), dimension(SIZE(KEY_ARRAY,1)) :: SegmentNumber
    logical (PGSLib_Log_Type), dimension(SIZE(KEY_ARRAY,1)) :: SegmentMask
    integer (PGSLib_Int_Type), dimension(SIZE(KEY_ARRAY,1)) :: Key_Bucket
    integer (PGSLib_Int_Type), dimension(SIZE(KEY_ARRAY,1)) :: Key_Rank
    integer (PGSLib_Int_Type), dimension(SIZE(KEY_ARRAY,1)) :: Key_Local
    integer (PGSLib_Int_Type), dimension(SIZE(KEY_ARRAY,1)) :: Global_Index
    integer (PGSLib_Int_Type), dimension(SIZE(KEY_ARRAY,1)) :: Bucket_Rank
    integer (PGSLib_Int_Type), &
         &   POINTER,          &
         &  dimension(:)       :: Supp_Key, Supp_SegNum, Supp_Bucket, Supp_Rank, Supp_GI
    integer (PGSLib_Int_Type), &
         &   POINTER,          &
         &  dimension(:)       :: Work_Key, Work_SegNum, Work_Rank, Work_Enum, Work_GI, GI_Rank
    

    type (PGSLIB_GS_Trace), POINTER :: Bucket_Trace
    integer (PGSLib_Int_Type) :: N_Supplement, N_Duplicate
    integer (PGSLib_Int_Type) :: N_Local, N_Work, MinIndex
    integer :: k, I, Local_I, Global_I
    character (LEN=256):: out_string
    

    ! If all the arrays are empty (global array size = 0) then there is nothing to do
    ! So we can just return.
    if (PGSLib_Global_All(SIZE(Key_Array,1) == 0)) then
       ! Should not need to set the result, but this keeps compilers from complaining
       GRADE_UP_Int_1D = 0
       return
    end if

    ! Check that we got valid argument combinations
    if (PRESENT(DIM)) then
       if (DIM /= 1) then
          call PGSLib_Fatal_ERROR('In PGSLIB_GRADE_UP, if DIM is present it must == 1')
       end if
    end if
    
    ! If we have a segment we need to compute a segment number.
    if (PRESENT(Segment)) then
! This fails if SIZE(Segment) == 0
!       SegmentMask   = PGSLib_Global_EOSHIFT(Segment, SHIFT=-1, BOUNDARY= (.NOT.Segment(1)))
! Instead, first use a CSHIFT, which gets every item correct exect (possibly) the first.
       SegmentMask   = PGSLib_Global_CSHIFT(Segment, SHIFT=-1)
! Now fix the first item.  We do not know which is the first item though, since some procs may have
! arrays of size 0.       
! What we need to do is find the value of the lowest *global* index.
       Global_Index = PGSLib_SUM_PREFIX( (/ (1, K=1, SIZE(Global_Index,1)) /) )
       MinIndex    = PGSLib_Global_MINVAL(Global_Index)
       WHERE (Global_Index == MinIndex)
          SegmentMask = .NOT. Segment
       END WHERE

       SegmentMask   = .NOT.(Segment .EQV. SegmentMask)
! Now we have flagged the beginning of each segment with a .TRUE.
       SegmentNumber = MERGE(1,0, SegmentMask)
       SegmentNumber = PGSlib_SUM_PREFIX(SegmentNumber)
    else
       SegmentNumber = 1
    end if

    ! We need to determine the buckets to use
    Call Determine_Buckets(Key_array, SegmentNumber, Buckets)

    ! Assign a bucket to each key
    DO k=1,SIZE(Key_array,1)
       Key_Bucket(k) = BucketNumber(Key_array(k), SegmentNumber(k), Buckets)
#ifdef DEBUG_SORT
       write(out_string,*) 'K, Key_Array(k), SegmentNumber(k), Key_Bucket(k) = ', &
            &               K, Key_Array(k), SegmentNumber(k), Key_Bucket(k)
       call pgslib_output(out_string)
#endif
    END DO
#ifdef DEBUG_SORT
    WRITE(out_string,*) 'MINVAL(Key_Bucket), MAXVAL(Key_Bucket) =', MINVAL(Key_Bucket), MAXVAL(Key_Bucket)
    call pgslib_output(out_string)
#endif    


    ! Sort the data in buckets
    Bucket_Rank = PGSLib_GRADE_UP_Local(Key_Bucket)
#ifdef DEBUG_SORT
    WRITE(out_string,*) 'MINVAL(Bucket_rank), MAXVAL(Bucket_rank) =', MINVAL(Bucket_rank), MAXVAL(Bucket_rank)
    call pgslib_output(out_string)
#endif
    
    ! Enumerate keys to assure stability of rank
    Global_Index = PGSLib_SUM_PREFIX( (/ (K, K=1, SIZE(Global_Index,1)) /) )

    Key_Local(Bucket_Rank)     = Key_Array
    SegmentNumber(Bucket_Rank) = SegmentNumber
    Key_Bucket(Bucket_Rank)    = Key_Bucket
    Global_Index(Bucket_Rank)  = Global_Index

    ! The off-pe sort is: send buckets, local sort, enumerate, return
    ! How much data goes off-pe?
    N_Supplement = COUNT(Key_Bucket /= ThisBucket())
    ALLOCATE(Supp_Key(N_Supplement),   &
         &   Supp_SegNum(N_Supplement),&
         &   Supp_Bucket(N_Supplement),&
         &   Supp_Rank(N_Supplement),  &
         &   Supp_GI(N_Supplement)     )
    Supp_Key    = PACK(Key_Local,     MASK=(Key_Bucket /= ThisBucket()))
    Supp_SegNum = PACK(SegmentNumber, MASK=(Key_Bucket /= ThisBucket()))
    Supp_Bucket = PACK(Key_Bucket,    MASK=(Key_Bucket /= ThisBucket()))
    Supp_GI     = PACK(Global_Index,  MASK=(Key_Bucket /= ThisBucket()))
    Supp_Rank = 0

    ! Setup the trace for the send-to-queue
    Bucket_Trace => PGSlib_Setup_Basic_Trace(Supp_Bucket)
    N_Duplicate = Bucket_Trace%N_Duplicate

    ! Make work space, now that we know how much data we will get
    N_Local = COUNT(Key_Bucket == ThisBucket())
    N_Work  = N_Local + N_Duplicate
#ifdef DEBUG_SORT
    WRITE(out_string,*) 'N_Local, N_Work = ', N_Local, N_Work
    call pgslib_output(out_string)
    call pgslib_flush_output()
#endif

    ALLOCATE(Work_Key(N_Work),     &
         &   Work_SegNum(N_Work),  &
         &   Work_Rank(N_Work),    &
         &   Work_Enum(N_Work),    &
         &   Work_GI(N_Work),      &       
         &   GI_Rank(N_Work)       )
         

    ! Put local stuff into work arrays through local move
    Work_Key(1:N_Local)    = PACK(Key_Local,     MASK=(Key_Bucket == ThisBucket()))
    Work_SegNum(1:N_Local) = PACK(SegmentNumber, MASK=(Key_Bucket == ThisBucket()))
    Work_GI(1:N_Local)     = PACK(Global_Index,  MASK=(Key_Bucket == ThisBucket()))

    ! Stuff coming in from other PEs uses send-to-queue
    
    Work_Key(N_Local+1:)    = pgslib_scatter_buffer(Supp_Key,    Bucket_Trace)
    Work_SegNum(N_Local+1:) = pgslib_scatter_buffer(Supp_SegNum, Bucket_Trace)
    Work_GI(N_Local+1:)     = pgslib_scatter_buffer(Supp_GI,     Bucket_Trace)

    ! Local Rank is the easiest part
#ifdef DEBUG_SORT
    WRITE(out_string,*) 'Calling PGSLib_GRADE_UP_Local'
    call pgslib_output(out_string)
    call pgslib_flush_output()
#endif

    ! First rank and sort according to global index
    GI_Rank   = PGSLib_GRADE_UP_LOCAL(Work_GI, Work_SegNum)
    Work_Key   (GI_Rank) = Work_Key
    Work_SegNum(GI_Rank) = Work_SegNum

    ! Now rank, and rank is stable so will preserve GI ordering
    Work_Rank = PGSlib_GRADE_UP_Local(Work_Key, Work_SegNum)
    ! This rank has to be re-ordered so that it represents total rank
    Work_Rank = Work_Rank(GI_Rank)

    ! Global Rank is requires global numbering
#ifdef DEBUG_SORT
    WRITE(out_string,*) 'Returned from PGSLib_GRADE_UP_Local'
    call pgslib_output(out_string)
    call pgslib_flush_output()
#endif

    Work_SegNum = 1
#ifdef DEBUG_SORT
    WRITE(out_string,*) 'Calling SUM_PREFIX'
    call pgslib_output(out_string)
#endif
    Work_Enum   = PGSLib_SUM_PREFIX(Work_SegNum)
#ifdef DEBUG_SORT
    WRITE(out_string,*) 'MINVAL(Work_enum), MAXVAL(Work_enum) =', MINVAL(Work_enum), MAXVAL(Work_enum)
    call pgslib_output(out_string)
#endif
    
    ! Global rank is just (local) permuation of global counting
    Work_Rank = Work_Enum(Work_Rank)
#ifdef DEBUG_SORT
    WRITE(out_string,*) 'MINVAL(Work_Rank), MAXVAL(Work_Rank) =', MINVAL(Work_Rank), MAXVAL(Work_Rank)
    call pgslib_output(out_string)
#endif

    ! First step towards getting the stuff home is to bring home the work buffers
    
    Supp_Rank =  pgslib_gather_buffer(Work_Rank(N_Local+1:), Bucket_Trace)

    ! Unpack the rank 
    Local_I  = 1
    Global_I = 1
    DO I = 1, SIZE(Key_Rank,1)
       if (Key_Bucket(I) == ThisBucket()) then
          Key_Rank(I) = Work_Rank(Local_I)
          Local_I = Local_I+1
       else
          Key_Rank(I) = Supp_Rank(Global_I)
          Global_I = Global_I + 1
       endif
    ENDDO

#ifdef DEBUG_SORT
    WRITE(out_string,*) 'MINVAL(Key_Rank), MAXVAL(Key_Rank) =', MINVAL(Key_Rank), MAXVAL(Key_Rank)
    call pgslib_output(out_string)
    call pgslib_flush_output()
#endif
    ! Finally, unscramble  Key_Rank to get funciton result
    GRADE_UP_Int_1D = Key_Rank(Bucket_Rank)

#ifdef DEBUG_SORT
    WRITE(out_string,*) 'Releasing Bucket_Trace'
    call pgslib_output(out_string)
    call pgslib_flush_output()
#endif
    call PGSLib_GS_Release_Trace(Bucket_Trace)

#ifdef DEBUG_SORT
    WRITE(out_string,*) 'Deallocating Arrays'
    call pgslib_output(out_string)
    call pgslib_flush_output()
#endif

    DEALLOCATE(Buckets)
    IF (ASSOCIATED(Supp_Rank)) then
       DEALLOCATE(Supp_Rank, Supp_Bucket, Supp_SegNum, Supp_Key, Supp_GI)
    end IF
    if (ASSOCIATED(Work_Key)) then
       DEALLOCATE(Work_Key, Work_SegNum, Work_Rank, Work_Enum, Work_GI, GI_Rank)
    end if
  
    RETURN

  END function GRADE_UP_Int_1D

  function ThisBucket()
    implicit none
    integer (PGSLib_Int_Type) :: ThisBucket
    ThisBucket = PGSLib_Inquire_thisPE()
    return
  end function ThisBucket

  subroutine Determine_Buckets_Int_1D(KEY_Array, SEGMENTNUMBER, Buckets)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This function picks Samples random items from Key_Array and returns 
      ! Bucket ranges.
    USE PGSLib_Utility_MODULE
    implicit none
      
      ! Arguments
      integer, intent(IN   )              , &
           &   dimension(               :) :: Key_Array
      integer, intent(IN   ),               &
           &   optional,                    &
           &   dimension(SIZE(KEY_ARRAY,1)):: SEGMENTNUMBER
      type (INT_BUCKET_TYPE),               &
           &   POINTER,                     &
           &   dimension(               :) :: Buckets

      ! Local variables
      integer :: N_Keys, N_Buckets, N_Samples, N_Samples_Tot
      integer :: BucketSize, b, High_index, Low_index
      real    :: SampleCutOff
      integer, pointer, dimension(:) :: Samples, Samples_Tot
      integer, pointer, dimension(:) :: Sample_SegNums, Sample_SegNums_Tot
      real   ,          dimension(SIZE(Key_array,1)) :: Sample_Seed
      logical,          dimension(SIZE(Key_array,1)) :: Sample_Mask
      
      character (LEN=256):: out_string

      ! Choose some number of keys for the sample
      N_Keys    = SIZE(KEY_ARRAY,1)
      N_Samples = MIN(N_Keys,MAX(N_Keys/( MAX(50, PGSLib_Inquire_nPE()) ),50) )
      

      ! Randomly choose about that many samples.  Count after selecting to 
      ! determine exact number chosen.
      Call RANDOM_SEED()
      Call RANDOM_NUMBER(HARVEST = Sample_Seed)
      if( N_Keys > 0) then 
         SampleCutOff = REAL(N_Samples)/REAL(N_Keys)
      else
         SampleCutOff = 0.0
      end if
      Sample_Mask = Sample_Seed < SampleCutOff
      N_Samples   = COUNT(Sample_Mask)

      N_Samples_Tot = PGSLib_Global_Sum(N_Samples)

      ALLOCATE(Samples       (N_Samples),     &
           &   Sample_SegNums(N_Samples))
      if (PGSLib_Inquire_IO_P()) then
         ALLOCATE(Samples_Tot(N_Samples_Tot), &
              &   Sample_SegNums_Tot(N_Samples_Tot))
      else
         ALLOCATE(Samples_Tot(1),             &
              &   Sample_SegNums_Tot(1))
      end if
      
      Samples        = PACK(Key_Array, MASK=Sample_Mask)
      Sample_SegNums = PACK(SegmentNumber, MASK=Sample_Mask)
      call pgslib_collate(Samples_Tot, Samples)
      call pgslib_collate(Sample_SegNums_Tot, Sample_SegNums)

      ! Sort the samples
      Samples_Tot(       PGSLib_GRADE_UP_LOCAL(Samples_Tot, SEGNUM=Sample_SegNums_Tot)) = Samples_Tot
      Sample_SegNums_Tot(PGSLib_GRADE_UP_LOCAL(Samples_Tot, SEGNUM=Sample_SegNums_Tot)) = Sample_SegNums_Tot

      ! One bucket per processor is reasonable.
      N_Buckets = PGSLib_Inquire_nPE()
      ALLOCATE(Buckets(N_Buckets))

      Buckets(1)%Low%Data    = PGSLib_Global_MINVAL(Key_Array)
      Buckets(1)%Low%SegNum  = PGSLib_Global_MINVAL(SegmentNumber)
      b = N_Buckets
      Buckets(b)%High%Data   = PGSLib_Global_MAXVAL(Key_Array) + 1
      Buckets(b)%High%SegNum = PGSLib_Global_MAXVAL(SegmentNumber) + 1


      if (PGSLib_Inquire_IO_P()) then

                                ! If multiple PEs then have to define buckets for each PE
         if (N_Buckets > 1) then
            BucketSize = N_Samples_Tot/N_Buckets
            High_index = MIN(BucketSize, N_Samples_Tot)         
            Buckets(1)%High%Data   = MAX(Samples_Tot(High_Index),Buckets(1)%Low%Data+1)
            Buckets(1)%High%SegNum = Sample_SegNums_Tot(High_Index) 


            b = 1
#ifdef DEBUG_SORT
            write(out_string,*) 'b, Low-D, High-D, Low-S, High-S = ', b, &
            & Buckets(b)%Low%Data, Buckets(b)%High%Data, &
            & Buckets(b)%Low%SegNum, Buckets(b)%High%SegNum
            call pgslib_output(out_string)
#endif
            DO b = 2, N_Buckets - 1
               Low_index  = Min(High_Index + 1, N_Samples_Tot)
               High_Index = Min(Low_index + BucketSize - 1, N_Samples_Tot)
               Buckets(b)%Low%Data    = Buckets(b-1)%High%Data
               Buckets(b)%Low%SegNum  = Buckets(b-1)%High%SegNum
               Buckets(b)%High%Data   = MAX(Samples_Tot(High_index),Buckets(b)%Low%Data+1)
               Buckets(b)%High%SegNum = Sample_SegNums_Tot(High_index)
#ifdef DEBUG_SORT
               write(out_string,*) 'b, Low-D, High-D, Low-S, High-S = ', b, &
               & Buckets(b)%Low%Data, Buckets(b)%High%Data, &
               & Buckets(b)%Low%SegNum, Buckets(b)%High%SegNum
               call pgslib_output(out_string)
#endif
            END DO

            b = N_Buckets
            Buckets(b)%Low%Data    = Buckets(b-1)%High%Data
            Buckets(b)%Low%SegNum  = Buckets(b-1)%High%SegNum
         endif
         b = N_Buckets
#ifdef DEBUG_SORT
         write(out_string,*) 'b, Low-D, High-D, Low-S, High-S = ', b, &
         & Buckets(b)%Low%Data, Buckets(b)%High%Data, &
         & Buckets(b)%Low%SegNum, Buckets(b)%High%SegNum
         call pgslib_output(out_string)
#endif
      endif
      call pgslib_bcast(Buckets%Low%Data   )
      call pgslib_bcast(Buckets%High%Data  )
      call pgslib_bcast(Buckets%Low%SegNum )
      call pgslib_bcast(Buckets%High%SegNum)

      DEALLOCATE(Samples, Sample_SegNums, Samples_Tot, Sample_SegNums_Tot)
      RETURN
    END subroutine Determine_Buckets_Int_1D

    function GE_Bucket_Int(A, B)
      implicit none
      type (INT_BUCKET_BNDRY),   intent(IN) :: A, B
      logical (PGSlib_Log_Type)             :: GE_Bucket_Int

      GE_Bucket_Int =  (A%SegNum >  B%SegNum) .OR.   &
           &          ((A%SegNum == B%SegNum) .AND. (A%Data >= B%Data) )

      RETURN
    END function GE_Bucket_Int

    function LE_Bucket_Int(A, B)
      implicit none
      type (INT_BUCKET_BNDRY),   intent(IN) :: A, B
      logical (PGSlib_Log_Type)             :: LE_Bucket_Int

      LE_Bucket_Int =  (A%SegNum <  B%SegNum) .OR.   &
           &          ((A%SegNum == B%SegNum) .AND. (A%Data <= B%Data) )

      RETURN
    END function LE_Bucket_Int

    function LT_Bucket_Int(A, B)
      implicit none
      type (INT_BUCKET_BNDRY),   intent(IN) :: A, B
      logical (PGSlib_Log_Type)             :: LT_Bucket_Int

      LT_Bucket_Int =  (A%SegNum <  B%SegNum) .OR.   &
           &          ((A%SegNum == B%SegNum) .AND.(A%Data < B%Data) )

      RETURN
    END function LT_Bucket_Int

    function EQ_Bucket_Int(A, B)
      implicit none
      type (INT_BUCKET_BNDRY),   intent(IN) :: A, B
      logical (PGSlib_Log_Type)             :: EQ_Bucket_Int

      EQ_Bucket_Int = (A%SegNum == B%SegNum) .AND.   &
           &          (A%Data   == B%Data)

      RETURN
    END function EQ_Bucket_Int

    function BucketNumber_Int(Data, SegNum, Buckets)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Search through Buckets to find which bucket contains (Data,SegNum)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none
      integer (PGSLib_Int_Type), intent(IN   ) :: Data
      integer (PGSLib_Int_Type), intent(IN   ) :: SegNum
      type    (INT_BUCKET_TYPE), intent(IN   ), &
           &                      dimension(:) :: Buckets
      integer (PGSLib_Int_Type)                :: BucketNumber_Int

      !Local Variables
      integer (PGSLib_Int_Type), SAVE :: Last_Bucket = 1
      integer (PGSLib_Int_Type)       :: LowB, HighB, MidB
      integer (PGSLib_Int_Type)       :: b
      type (INT_BUCKET_BNDRY)         :: ThisBucket

      ThisBucket = INT_BUCKET_BNDRY(Data, SegNum)

      ! Use binary search of array buckets
      ! But first check last bucket used
      if (  (ThisBucket >= Buckets(Last_Bucket)%Low ) .AND. &
           &(ThisBucket <  Buckets(Last_Bucket)%High)) then
         BucketNumber_Int = Last_Bucket
         RETURN
      end if

      ! If not in last bucket, then perform search, starting from that 
      ! bucket.
      LowB  = 1
      HighB = SIZE(Buckets,1)
      MidB  = Last_Bucket

      DO b = 1, SIZE(Buckets,1)  
         if (ThisBucket <  Buckets(MidB)%High) then
            HighB = MidB
            if (ThisBucket <  Buckets(LowB)%High) then
               Last_Bucket = LowB
               exit
            endif
         else
            LowB  = MidB
            if (ThisBucket >= Buckets(HighB)%Low) then
               Last_Bucket = HighB
               exit
            endif
         end if

         MidB = LowB + (HighB - LowB)/2
         if ( (ThisBucket >= Buckets(MidB)%Low ) .AND. &
            & (ThisBucket <  Buckets(MidB)%High) ) then
            Last_Bucket = MidB
            exit
         end if
      END DO
      BucketNumber_Int = Last_Bucket
      RETURN
    END function BucketNumber_Int
              
         

      


    function GE_Nodes_Int(A, B)
      implicit none
      type (WorkNode_Int), intent(IN   ) :: A, B
      logical :: GE_Nodes_Int

      GE_Nodes_Int =  ( A%SegNum >  B%SegNum) .OR.  &  ! Different Segments
           &         (( A%SegNum == B%SegNum) .AND. &  ! Same Segment
           &          ( A%Data   >  B%Data))  .OR.  &  ! Different Data
           &         (( A%SegNum == B%SegNum) .AND. &  ! Same Segment
           &          ( A%Data   == B%Data)   .AND. &  ! Same value
           &          ( A%Pos    >  B%Pos))            ! Force stability

      return
    end function GE_Nodes_Int

    function Rank_Up_Int_1D(Key, SegNum)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Given array Key, return in Index the index of each item
      ! If SegNum is present it is the associated segment number for each item
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      integer, intent(IN   ), dimension(          :) :: Key
      integer, intent(IN   ), OPTIONAL,    &
           &                  TARGET,      &
           &                  dimension(SIZE(Key,1)) :: SegNum

      integer,                dimension(SIZE(Key,1)) :: Rank_Up_Int_1D

      ! Local variables
      integer :: N, node, Temp_Index
      integer, pointer,       dimension(          :) :: L_SegNum
      integer,                dimension(size(Key,1)) :: L_Index

      N = SIZE(Key,1)
      L_Index = (/ (node, node=1,N) /)

      if (PRESENT(SegNum)) then
         L_SegNum => SegNum
      else
         ALLOCATE(L_SegNum(N))
         L_SegNum = 1
      end if
      
      
      ! First we build the heap
      DO node = N/2,1,-1
         call Heap_Insert(L_Index, node, Key, L_SegNum)
      enddo

      ! Now finish the indexing
      DO node = N,2,-1
         Temp_Index    = L_Index(1)
         L_Index(1)    = L_Index(node)
         L_Index(node) = Temp_Index
         call Heap_Insert(L_Index(1:node-1), 1, Key, L_SegNum)
      end DO


      if (N > 1) then
         Rank_Up_Int_1D(L_Index) = (/ (node, node=1,N) /)
      elseif (N == 1) then
         Rank_Up_Int_1D(1) = 1
      end if


      return
    end function Rank_Up_Int_1D
    

      

    subroutine Heap_Insert_Int(Position, Node, Heap, SegNum)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This routine takes Node and modifies Position so that Heap(Position(:)) is a heap
      ! 
      ! Input:
      !    Heap(1:N) - array of data
      !   Position(  :) - array, of size heap.  May be smaller than array Heap.
      !                This is the permuatation array
      !                which tells how to re-order Heap to get a heap.
      !                Position is the size of the heap
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      integer, intent(IN   )               :: Node
      integer, intent(INOUT), dimension(:) :: Position
      integer, intent(IN   ), dimension(:) :: Heap
      integer, intent(IN   ), dimension(:) :: SegNum

      ! Local variables
      integer :: HeapSize, TreeHeight, Temp, l
      type (WorkNode_Int), TARGET  :: Current, Left, Right
      type (WorkNode_Int), Pointer :: Child, Largest

      HeapSize      = SIZE(Position,1)
      Current%Index = Node
      TreeHeight    = Log(REAL(SIZE(Heap)))/Log(2.0) + 1

      DO l = 1, TreeHeight  ! At most Log_2(N) steps (N=Size(HEAP) = total number data elements)

         Current%Pos   = Position(Current%Index)
         Current%Data  = Heap(Current%Pos)
         Current%SegNum= SegNum(Current%Pos)

         Left%Index    = 2*Current%Index

         Right%Index   = Left%Index + 1
         
         if (Right%Index <= HeapSize) then ! If there is a right child,
            Right%Pos     = Position(Right%Index)
            Right%Data    = Heap(Right%Pos)
            Right%SegNum  = SegNum(Right%Pos)
            Left%Pos      = Position(Left%Index)
            Left%Data     = Heap(Left%Pos)
            Left%SegNum   = SegNum(Left%Pos)
            if (Right .GE. Left) then      ! then there is a left child, and pick the larger
               Child => Right
            else
               Child => Left
            endif
         else
            if (Left%Index <= HeapSize) then ! If there is only a left child, compare against it only
               Left%Pos      = Position(Left%Index)
               Left%Data     = Heap(Left%Pos)
               Left%SegNum   = SegNum(Left%Pos)
               Child => Left
            else
               exit                          ! If there are no children, then we are done.
            endif
         end if

         if (Current .GE. Child) exit ! If Current is .GE. both children, done
         
                                      ! Otherwise, swap current and child
         Temp                      = Position(Current%Index)
         Position(Current%Index) = Position(Child%Index)
         Position(Child%Index)   = Temp
         Current%Index           = Child%Index
      END DO

      RETURN
    END subroutine Heap_Insert_Int
      


  END MODULE PGSLIB_GRADE_MODULE
