MODULE PARTITION_DATA_TYPES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE
  !   Provide the data types and support routines for paritioned
  !   sets and graphs.
  !
  ! $Id: partition_data_types.F,v 1.6 2002/01/11 04:23:06 lally Exp $
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use pgslib_Globals_Module
  use pgslib_io_module,   only: PGSLib_Collate,  PGSLib_BCast
  use pgslib_reductions_module, only: PGSLib_Global_SUM
  use pgslib_scan_module,       only: PGSLib_SUM_PREFIX
  use pgslib_utility_module, only : PGSLib_Fatal_Error, PGSLib_Inquire_thisPE
  use Partition_Constants

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: A_SET
  PUBLIC :: A_SET_PARTITIONING
  PUBLIC :: PARTITION_TO_PROCESSOR_MAP
  PUBLIC :: INITIALIZE
  PUBLIC :: ALLOC
  PUBLIC :: FREE
  PUBLIC :: SIZE
  PUBLIC :: SET
  PUBLIC :: LOOKUP
  PUBLIC :: TRANSLATE
  PUBLIC :: Get_Partition_ID
  PUBLIC :: Get_Num_Partitions
  PUBLIC :: Get_Num_Partitions_Available
  PUBLIC :: Get_Start
  PUBLIC :: Get_End
  PUBLIC :: Get_Sizes
  PUBLIC :: Get_Start_Available
  PUBLIC :: Get_End_Available
  PUBLIC :: Get_Sizes_Available
  PUBLIC :: Partition_Is_Available
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Type definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type A_SET
     PRIVATE
     ! The sets allowed here are assumed to be "in order". 
     ! Permutations will come later, but and always be distinct from 
     ! partitioning.
     ! At the simplest a set is isomorphic to the integers [1..Size_Of_Set]
     integer(PGSLib_INT_TYPE) :: Size_Of_Set

     ! Our sets are often distributed.  We need to know if they are since
     ! cross-set operations on distributed sets need special treatment.
     logical(PGSLib_Log_Type) :: Set_Is_Distributed

     ! If the sets are distributed, then they must be mapped to processor.
     ! For now we assume that the processor mapping is [1..Number_Of_Processors]
     ! A physical paritioning has a 1-to-1 and onto mapping from partitons to processors.
     ! The partition numbering is also assumed to go [1..Number_Of_Processors]
     ! If the set is NOT distributed then the PHYSICAL_PARTITION may be NULL,
     ! and in any case may not be accessed.

     type (A_SET_PARTITIONING),&
             POINTER          :: PHYSICAL_PARTITION
  end type A_SET

  type PARTITION_TO_PROCESSOR_MAP
     PRIVATE
     ! This is an array of size total-number-of-partitions
     ! For any partition p, MAP(p) returns the processor which
     ! "owns" that partition.  Or, the object for which that
     ! partition is "available".  The map is identical for all objects.
     integer(PGSLib_INT_TYPE), &
          dimension(:)       , &
          POINTER              :: MAP
  end type PARTITION_TO_PROCESSOR_MAP
  
  type A_SET_PARTITIONING
     PRIVATE
     ! First of all, we expect a distributed model, so give each copy
     ! of the partitioning object a unique identifier.  In our simple
     ! distributed computing model we can just use processor number.
     integer(PGSLib_Int_Type) :: This_Processor

     ! This is useful since it allows us to get info such as total set size
     type(A_SET),              &
             POINTER          :: Set_That_Is_Partitioned

     ! How many partitions (this is *total*, of course, since partitionings don''t 
     ! actually know about local vs global, except in special circumstances
     integer(PGSLib_INT_TYPE) :: Number_Of_Partitions

     ! How many partitions available to the calling object.
     integer(PGSLib_INT_TYPE) :: Number_Of_Partitions_Available

     ! This array has length Number_Of_Partitions (once the data structure is loaded).
     ! Each item is the number of elements in that partition.  Note that
     ! SUM(Partition_Sizes) == Set_That_Is_Partitioned%Size_Of_Set, by definition.
     integer(PGSLib_INT_TYPE), &
             dimension(:)    , &
             POINTER          :: Partition_Sizes
     
     ! These arrays contain the first and last item in each of the partitions.  Since the items
     ! must be in order, these arrays are redundant.  But it is convenient having them both.
     integer(PGSLib_INT_TYPE), &
             dimension(:)    , &
             POINTER          :: Partition_Start
     integer(PGSLib_INT_TYPE), &
             dimension(:)    , &
             POINTER          :: Partition_End

     ! This tells how the partitions are mapped to the processors.  That is,
     ! which partition is available to which partitioning object.
     ! This array has length Number_Of_Partitions_Available (once the data structure is loaded).
     ! Each item is the number of elements in that partition.  
     integer(PGSLib_INT_TYPE), &
             dimension(:)    , &
             POINTER          :: Available_Partition_Sizes

     ! These arrays contain the first and last item in each of the partitions.  Since the items
     ! must be in order, these arrays are redundant.  But it is convenient having them both.
     integer(PGSLib_INT_TYPE), &
             dimension(:)    , &
             POINTER          :: Available_Partition_Start
     integer(PGSLib_INT_TYPE), &
             dimension(:)    , &
             POINTER          :: Available_Partition_End

     ! We often need to translate from an available partition number to a global partition number.
     ! This is a simple translation, but requires global communication.  It is a small amount of
     ! data, so store it here.
     integer(PGSLib_INT_TYPE), &
             dimension(:)    , &
             POINTER          :: Available_Global_Part_Num

     ! This indicates whether the data is in place for fast partition lookups.
     logical(PGSLib_Log_Type) :: USE_FAST_LOOKUP
     ! This is the indexing array to use for fast lookups
     integer(PGSLib_Int_Type), &
             dimension(:)    , &
             POINTER          :: FAST_LOOKUP_TABLE

     ! This tells how the partitions are mapped to the processors.  That is,
     ! which partition is available to which partitioning object.
     type(PARTITION_TO_PROCESSOR_MAP) :: MAP

!!$     ! The partitioning is a collection of Number_Of_Partions.  Each item in the collection
!!$     ! is a collection of some number of items from the set that is partitioned.  It is often useful
!!$     ! to think of the union of all the partitions as a set itself.  That is one good way
!!$     ! of describing how partitions are distributed, since that Set_Of_Partitions may be a
!!$     ! distributed set.
!!$     type (A_SET)             :: Set_Of_Partitions
!!$
!!$     ! If the partition is physical then there is a 1-to-1, onto mapping from elements in Set_Of_Partitions
!!$     ! to the number of processors.
!!$     logical(PGSLib_Log_Type) :: Partition_Is_Physical
  end type A_SET_PARTITIONING

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! INTERFACEs for generic procedure names           !!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE INITIALIZE
     MODULE PROCEDURE A_Set_INIT
     MODULE PROCEDURE A_Set_Part_INIT
     MODULE PROCEDURE Map_INIT
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE A_Set_ALLOC
     MODULE PROCEDURE A_Set_Part_ALLOC
     MODULE PROCEDURE Map_ALLOC
  END INTERFACE

  INTERFACE FREE
     MODULE PROCEDURE A_Set_FREE
     MODULE PROCEDURE A_Set_Part_FREE
     MODULE PROCEDURE Map_FREE
  END INTERFACE

  INTERFACE SIZE
     MODULE PROCEDURE A_Set_Size
  END INTERFACE

  INTERFACE SET
     MODULE PROCEDURE A_Set_Part_SET
     MODULE PROCEDURE Map_SET_From_Input
     MODULE PROCEDURE Map_SET_By_Computing
  END INTERFACE
  
  INTERFACE LOOKUP
     MODULE PROCEDURE Map_Lookup
     MODULE PROCEDURE PARTITION_LOOKUP
  END INTERFACE

  INTERFACE TRANSLATE
     MODULE PROCEDURE Translate_Local_To_Global_Part
  END INTERFACE

  INTERFACE Get_Partition_ID
     MODULE PROCEDURE Get_Partition_ID
  END INTERFACE

  INTERFACE Get_Num_Partitions
     MODULE PROCEDURE Number_Of_Partitions
  END INTERFACE

  INTERFACE Get_Num_Partitions_Available
     MODULE PROCEDURE Number_Of_Partitions_Available
  END INTERFACE
  
  INTERFACE Get_Start
     MODULE PROCEDURE Get_Start
  END INTERFACE

  INTERFACE Get_End
     MODULE PROCEDURE Get_End
  END INTERFACE

  INTERFACE Get_Sizes
     MODULE PROCEDURE Get_Sizes
  END INTERFACE

  INTERFACE Get_Start_Available
     MODULE PROCEDURE Get_Start_Available
  END INTERFACE

  INTERFACE Get_End_Available
     MODULE PROCEDURE Get_End_Available
  END INTERFACE

  INTERFACE Get_Sizes_Available
     MODULE PROCEDURE Get_Sizes_Available
  END INTERFACE

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routines for maintaining A_SET data structures                     !
  ! Routines in this section (Generic name = specific name)
  !          INITIALIZE = A_Set_INIT
  !          ALLOC      = A_Set_ALLOC
  !          FREE       = A_Set_FREE
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine A_Set_INIT(The_Set)
    implicit none
    type(A_Set), intent(OUT) :: The_Set
      
    The_Set%Size_Of_Set = PARTITION_INVALID_SIZE
    The_Set%Set_Is_Distributed = .FALSE.
    NULLIFY(The_Set%Physical_Partition)

    RETURN
  END subroutine A_Set_INIT

  subroutine A_Set_ALLOC(The_Set, SIZE)
    implicit none
    type(A_Set), intent(INOUT) :: The_Set
    integer,     intent(IN   ) :: Size

    The_Set%Size_Of_Set = SIZE

    RETURN
  END subroutine A_Set_ALLOC

  subroutine A_Set_FREE(The_Set)
    implicit none
    type(A_Set), intent(INOUT) :: The_Set
    ! There isn''t any memory stored in here, so just initialize
    call INITIALIZE(The_Set)

    RETURN
  END subroutine A_Set_FREE

  function A_Set_Size(The_Set) RESULT(SIZE)
    implicit none
    type(A_Set), intent(IN   ) :: The_Set
    integer(PGSLib_Int_Type)   :: SIZE
    SIZE = The_Set%Size_Of_Set
    return
  end function A_Set_Size
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routines for maintaining PARTITION_TO_PROCESSOR_MAP
  ! Routines in this section (Generic name = specific name)
  !          INITIALIZE = Map_INIT
  !          ALLOC      = Map_ALLOC
  !          SET        = Map_SET
  !          FREE       = Map_FREE
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Map_INIT(MAP)
    implicit none
    type(PARTITION_TO_PROCESSOR_MAP), intent(INOUT) :: MAP

    ! The only thing to do is to nullify the pointer
    NULLIFY(MAP%MAP)
    return
  end subroutine Map_INIT
  

  subroutine Map_ALLOC(MAP, SIZE)
    implicit none
    type(PARTITION_TO_PROCESSOR_MAP), intent(INOUT) :: MAP
    integer(PGSLib_INT_Type),         intent(IN   ) :: SIZE

    ALLOCATE(MAP%MAP(SIZE))
    MAP%MAP = PARTITION_INVALID_PROCESSOR

    return
  end subroutine Map_ALLOC

  subroutine Map_SET_From_Input(MAP, PROCESSORS)
    implicit none
    type(PARTITION_TO_PROCESSOR_MAP), intent(INOUT) :: MAP
    integer(PGSLib_INT_Type),         dimension(:),  &
                                      intent(IN   ) :: PROCESSORS
    ! the map must have been initialized.
    ! If it has not been allocated, allocate it now.

    if (.NOT. ASSOCIATED(MAP%MAP)) then
       ALLOCATE(MAP%MAP(SIZE(PROCESSORS)))
    end if

    ! This form of the map_set routine takes a map as input
    MAP%MAP = PROCESSORS
    return
  end subroutine Map_SET_From_Input

  subroutine Map_SET_By_Computing(MAP, Available_Number_Of_Partitions)
    implicit none
    type(PARTITION_TO_PROCESSOR_MAP), intent(INOUT) :: MAP
    integer(PGSLib_INT_Type),         intent(IN   ) :: Available_Number_Of_Partitions

    ! This form of the map_set routine assumes that the
    ! partitions are in order on each processor, but we
    ! only have local information on each processor.  So we
    ! have to do some global operations to figure out the
    ! global map.

    ! Local variables
    integer :: Number_Of_Partitions
    integer, dimension(Available_Number_Of_Partitions) :: Local_Processor_Map

    ! If the map has not been allocated, allocate it now.  That
    ! will require a global operation.
    if (.NOT. ASSOCIATED(MAP%MAP)) then
       ! First we have to find the total number of partitons
       Number_Of_Partitions = PGSLib_Global_SUM(Available_Number_Of_Partitions)
       ALLOCATE(MAP%MAP(Number_Of_Partitions))
    else
       Number_Of_Partitions = SIZE(Map%Map)
    end if

    ! Now we assume that the partitions are in processor order
    Local_Processor_Map = PGSLib_Inquire_thisPE()

    ! Now cluster all the local maps into a single map
    call PGSLib_Collate(Map%Map, Local_Processor_Map)
    call PGSLib_BCast(Map%Map)
    
    return
  end subroutine Map_SET_By_Computing

  subroutine Map_FREE(MAP)
    implicit none
    type(PARTITION_TO_PROCESSOR_MAP), intent(INOUT) :: MAP

    DEALLOCATE(MAP%MAP)

    return
  end subroutine Map_FREE

  function Map_Lookup(Map, PARTITION)  RESULT(PROC)
    implicit none
    type(PARTITION_TO_PROCESSOR_MAP), intent(IN   ) :: MAP
    integer(PGSLib_INT_Type),         intent(IN   ) :: PARTITION
    integer(PGSLib_INT_Type)                        :: PROC

    PROC = Map%Map(PARTITION)
    return
  end function Map_Lookup
    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routines for maintaining A_SET_PARTITIONING data structures                     !
  ! Routines in this section (Generic name = specific name)
  !          INITIALIZE = A_Set_Part_INIT
  !          ALLOC      = A_Set_Part_ALLOC
  !          FREE       = A_Set_Part_FREE
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine A_Set_Part_INIT(The_Set_Part)
    implicit none
    type(A_SET_PARTITIONING), intent(OUT) :: The_Set_Part
      
    The_Set_Part%This_Processor = PARTITION_INVALID_PROCESSOR

    ! An initialized, but not allocated, partition does not point to any set
    NULLIFY(The_Set_Part%Set_That_Is_Partitioned)

    ! The number of partitions is invalid
    The_Set_Part%Number_Of_Partitions           = PARTITION_INVALID_PARTITIONS
    The_Set_Part%Number_Of_Partitions_Available = PARTITION_INVALID_PARTITIONS
    The_Set_Part%USE_FAST_LOOKUP = .FALSE.

    ! The size arrays do not exist
    NULLIFY(The_Set_Part%Partition_Sizes)
    NULLIFY(The_Set_Part%Partition_Start)
    NULLIFY(The_Set_Part%Partition_End)
    NULLIFY(The_Set_Part%Available_Partition_Sizes)
    NULLIFY(The_Set_Part%Available_Partition_Start)
    NULLIFY(The_Set_Part%Available_Partition_End)

!!$    ! The set of partitions is a non-existant set
!!$    call INITIALIZE(The_Set_Part%Set_Of_Partitions)
!!$
!!$    ! The partition is not physical, of course
!!$    The_Set_Part%Partition_Is_Physical = .FALSE.

    RETURN
  END subroutine A_Set_Part_INIT

  subroutine A_Set_Part_ALLOC(The_Set_Part, AVAILABLE_PARTITIONS, TOTAL_PARTITIONS)
    implicit none
    type(A_SET_PARTITIONING), intent(INOUT) :: The_Set_Part
    integer(PGSLib_INT_Type), intent(IN   ) :: AVAILABLE_PARTITIONS
    integer(PGSLib_INT_Type), OPTIONAL,      &
                              intent(IN   ) :: TOTAL_PARTITIONS

    ! Local variables
    integer :: Number_Of_Partitions
    
    ! Allocate space for the available partitions
    The_Set_Part%Number_Of_Partitions_Available = AVAILABLE_PARTITIONS
    ALLOCATE(The_Set_Part%Available_Partition_Sizes(Available_Partitions))
    ALLOCATE(The_Set_Part%Available_Partition_Start(Available_Partitions))
    ALLOCATE(The_Set_Part%Available_Partition_End  (Available_Partitions))
    ALLOCATE(The_Set_Part%Available_Global_Part_Num(Available_Partitions))

    ! If the total number of partitions was input, then do not need to do any further
    ! computation.  If not, need to calculate it, using GLOBAL operations.

    if (PRESENT(TOTAL_PARTITIONS)) then
       Number_Of_Partitions = TOTAL_PARTITIONS
    else
       Number_Of_Partitions = PGSLib_Global_SUM(AVAILABLE_PARTITIONS)
    end if
    The_Set_Part%Number_Of_Partitions = Number_Of_Partitions

    ! Finally, allocate the things which require knowing the total number of partitions
    ALLOCATE(The_Set_Part%Partition_Sizes(Number_Of_Partitions))
    ALLOCATE(The_Set_Part%Partition_Start(Number_Of_Partitions))
    ALLOCATE(The_Set_Part%Partition_End  (Number_Of_Partitions))
    
    ! The map needs to know about all the partitions
    call ALLOC(The_Set_Part%MAP, Number_Of_Partitions)
    return
  end subroutine A_Set_Part_ALLOC


  subroutine A_Set_Part_FREE(The_Set_Part)
    implicit none
    type(A_SET_PARTITIONING), intent(INOUT) :: The_Set_Part
      
    ! Deallocate any arrays which exist
    if (ASSOCIATED(The_Set_Part%Partition_Sizes)) NULLIFY(The_Set_Part%Partition_Sizes)
    if (ASSOCIATED(The_Set_Part%Partition_Start)) NULLIFY(The_Set_Part%Partition_Start)
    if (ASSOCIATED(The_Set_Part%Partition_End)  ) NULLIFY(The_Set_Part%Partition_End)
    if (ASSOCIATED(The_Set_Part%Available_Partition_Sizes)) NULLIFY(The_Set_Part%Available_Partition_Sizes)
    if (ASSOCIATED(The_Set_Part%Available_Partition_Start)) NULLIFY(The_Set_Part%Available_Partition_Start)
    if (ASSOCIATED(The_Set_Part%Available_Partition_End)  ) NULLIFY(The_Set_Part%Available_Partition_End)
    if (ASSOCIATED(The_Set_Part%Available_Global_Part_Num)) NULLIFY(The_Set_Part%Available_Global_Part_Num)

    ! And re-initialize
    call INITIALIZE(The_Set_Part)

    RETURN
  END subroutine A_Set_Part_FREE

  subroutine A_Set_Part_SET(The_Set_Part, PARTITIONED_SET, PARTITION_COLORS, SCOPE)
    implicit none
    type(A_SET_PARTITIONING), intent(INOUT) :: The_Set_Part
    type(A_SET), TARGET,      intent(IN   ) :: PARTITIONED_SET
    integer,     dimension(:),intent(IN   ) :: PARTITION_COLORS
    type(PGSLib_SCOPE),       intent(IN   ) :: SCOPE

    ! Local variables
    integer :: N_Available_Partitions, Total_Color_Size, ThisColor, e, i
    integer, dimension(:), POINTER :: TOT_PARTITION_COLORS


    ! For the moment, only LOCAL scope is allowed.  That means that all
    ! the colors that are input are only for available partitions.
    if (.NOT. (SCOPE == PGSLib_LOCAL)) then
       call PGSLib_Fatal_Error("SCOPE must be PGSLib_LOCAL in SET of A_SET_PARTITIONING")
       return ! never will get here, since fatal error abors
    end if

    ! Identify the object/processor
    The_Set_Part%This_Processor = PGSLib_Inquire_thisPE()

    ! Set the pinter to the set which is partitioned by this partitioning
    The_Set_Part%Set_That_Is_Partitioned => PARTITIONED_SET
    
    ! The color array must already be ordered.  Will relax that requirement eventually
    ! when I have a better setup procedure.  Not going to check that that condition
    ! is satisified.  User beware:)  Actually, only require that the colors
    ! be grouped.  Every transition from one color to another I assume is the boundary
    ! between partitions.

    ! Since colors are only for available partitions, only need for this object to
    ! parse its own color array.

    ! First off we need to count how many colors there are.
    N_Available_Partitions = 0
    ! We will look for color transitions.  The line below makes sure that
    ! the initial color is distinct from all valid colors.
    ThisColor = MINVAL(PARTITION_COLORS) - 1
    do e = 1, SIZE(PARTITION_COLORS)
       ! Is this the same as the previous color
       if (ThisColor /= PARTITION_COLORS(e)) then
          N_Available_Partitions = N_Available_Partitions + 1
          ThisColor = PARTITION_COLORS(e)
       end if
    end do

    The_Set_Part%Number_Of_Partitions_Available = N_Available_Partitions

    ! Global computation of total number of partitions is done internally
    call ALLOC(The_Set_Part, AVAILABLE_PARTITIONS = N_Available_Partitions)


    ! Now set up the start and end indices.
    call Setup_Start_End_Arrays(The_Set_Part%Available_Partition_Start, &
                                The_Set_Part%Available_Partition_End  , &
                                The_Set_Part%Available_Partition_Sizes, &
                                PARTITION_COLORS)

    ! We need to set up the same information for the global partition information.
    ! For that we need the total color array on each processor.
    Total_Color_Size = PGSLib_Global_SUM(SIZE(PARTITION_COLORS))
    ALLOCATE(TOT_PARTITION_COLORS(Total_Color_Size))
    call PGSLib_COLLATE(TOT_PARTITION_COLORS, PARTITION_COLORS)
    call PGSLib_BCAST(TOT_PARTITION_COLORS)

    call Setup_Start_End_Arrays(The_Set_Part%Partition_Start, &
                                The_Set_Part%Partition_End  , &
                                The_Set_Part%Partition_Sizes, &
                                TOT_PARTITION_COLORS)
    DEALLOCATE(TOT_PARTITION_COLORS)

    ! We also need a partition to processor map.  We only have local info, so use that form
    call SET(The_Set_Part%Map, N_Available_Partitions)

    ! Finally, we need to setup the local partition number to global partition number map.
    ! For our limited support of partitions, we know that partitions are numbered starting on the
    ! first processor.  
    The_Set_Part%Available_Global_Part_Num = PGSLib_SUM_PREFIX( (/ (1, i=1, N_Available_Partitions) /) )

    ! And now we are done
    return
  end subroutine A_Set_Part_SET


  subroutine Setup_Start_End_Arrays(Start, End, Sizes, Colors)
    implicit none
    ! Based on the input color array, which is assumed to be in order
    ! find the start and end of each color transition.  Also return the
    ! size.
    integer(PGSLib_Int_Type), dimension(:), intent(  OUT) :: Start
    integer(PGSLib_Int_Type), dimension(:), intent(  OUT) :: End
    integer(PGSLib_Int_Type), dimension(:), intent(  OUT) :: Sizes
    integer(PGSLib_Int_Type), dimension(:), intent(IN   ) :: Colors

    ! Local variables
    integer :: ThisColor, Partition, e

    ThisColor = COLORS(1)
    Partition = 1
    Start(Partition) = 1
    do e = 2, SIZE(COLORS)
       if (ThisColor /= COLORS(e)) then
          ! The end of the last partition was one iteration ago
          End(Partition) = e - 1
          ! The size of the last partition 
          Sizes(Partition) = End(Partition) - Start(Partition) + 1
          ! We''ve found another partition
          Partition = Partition + 1
          Start(Partition) = e
          ThisColor = COLORS(e)
       endif
    end do
    ! We didn''t get the final information for the last partition
    ! The end of the last partition was the end of the array
    e = SIZE(COLORS)
    End(Partition) = e
    ! The size of the last partition 
    Sizes(Partition) = End(Partition) - Start(Partition) + 1

    return
  end subroutine Setup_Start_End_Arrays

  function Get_Partition_ID(The_Set_Part) RESULT(P)
    implicit none
    type (A_SET_PARTITIONING), intent(IN) :: The_Set_Part
    integer(PGSLib_Int_Type)              :: P

    P = The_Set_Part%This_Processor
    
    return
  end function Get_Partition_ID
    

  function Number_Of_Partitions(PARTITION) RESULT(N)
    implicit none
    type (A_SET_PARTITIONING), intent(IN) :: PARTITION
    integer(PGSLib_Int_Type)             :: N

    N = PARTITION%Number_Of_Partitions
    RETURN
  END function Number_Of_Partitions

  function Number_Of_Partitions_Available(PARTITION) RESULT(N)
    implicit none
    type (A_SET_PARTITIONING), intent(IN) :: PARTITION
    integer(PGSLib_Int_Type)             :: N

    N = PARTITION%Number_Of_Partitions_Available
    RETURN
  END function Number_Of_Partitions_Available

  function Get_Start_Available(PARTITION) RESULT(Start)
    ! Return the array with the first item of each partition
    implicit none
    type(A_SET_PARTITIONING), intent(IN)  :: PARTITION
    integer(PGSLib_Int_Type), dimension(:),&
                              POINTER     :: Start

    Start => PARTITION%Available_Partition_Start
    RETURN
  END function Get_Start_Available

  function Get_End_Available(PARTITION) RESULT(End)
    ! Return the array with the last item of each partition
    implicit none
    type(A_SET_PARTITIONING), intent(IN)  :: PARTITION
    integer(PGSLib_Int_Type), dimension(:),&
                              POINTER     :: End

    End => PARTITION%Available_Partition_End
    RETURN
  END function Get_End_Available

  function Get_Start(PARTITION) RESULT(Start)
    ! Return the array with the first item of each partition
    implicit none
    type(A_SET_PARTITIONING), intent(IN)  :: PARTITION
    integer(PGSLib_Int_Type), dimension(:),&
                              POINTER     :: Start

    Start => PARTITION%Partition_Start
    RETURN
  END function Get_Start

  function Get_End(PARTITION) RESULT(End)
    ! Return the array with the last item of each partition
    implicit none
    type(A_SET_PARTITIONING), intent(IN)  :: PARTITION
    integer(PGSLib_Int_Type), dimension(:),&
                              POINTER     :: End

    End => PARTITION%Partition_End
    RETURN
  END function Get_End

  function Get_Sizes(PARTITION) RESULT(Size)
    ! Return the array with the size of each available partition
    implicit none
    type(A_SET_PARTITIONING), intent(IN)  :: PARTITION
    integer(PGSLib_Int_Type), dimension(:),&
                              POINTER     :: Size

    Size => PARTITION%Partition_Sizes
    RETURN
  END function Get_Sizes

  function Get_Sizes_Available(PARTITION) RESULT(Size)
    ! Return the array with the size of each available partition
    implicit none
    type(A_SET_PARTITIONING), intent(IN)  :: PARTITION
    integer(PGSLib_Int_Type), dimension(:),&
                              POINTER     :: Size

    Size => PARTITION%Available_Partition_Sizes
    RETURN
  END function Get_Sizes_Available

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Routines for extracting info from the partition  !!!!!!!!!!

  function Translate_Local_To_Global_Part(The_Set_Part, AvailablePartition) RESULT(GlobalPart)
    ! Returns the global partition number of the available partition number.
    ! (That is, partition 3 on process 4 may be global partition 17.)
    ! The input must be in the range [1..Number_Available_Partitions].  That is
    ! it must be a local partition number of the calling object.
    implicit none
    type(A_SET_PARTITIONING), intent(IN) :: The_Set_Part
    integer(PGSLib_INT_TYPE), intent(IN) :: AvailablePartition
    integer(PGSLib_Int_Type)             :: GlobalPart

    GlobalPart = The_Set_Part%Available_Global_Part_Num(AvailablePartition)
    RETURN
  END function Translate_Local_To_Global_Part

  function PARTITION_LOOKUP(The_Set_Partition, item) RESULT(PART_NUM)
    use pgslib_utility_module
    ! Returns the partition number of the partition which contains the item
    implicit none
    type(A_SET_PARTITIONING), intent(IN) :: The_Set_Partition
    integer(PGSLib_INT_TYPE), intent(IN) :: Item
    integer(PGSLib_Int_Type)             :: PART_NUM

    ! This routine has two options.  If we have pre-setup for a fast lookup, that
    ! means each object has a whole lookup array ready to dereference.  If we don''t have that,
    ! then we have to do the slow search.

    if (The_Set_Partition%USE_FAST_LOOKUP) then
       PART_NUM = FAST_PARTITION_LOOKUP(The_Set_Partition, item)
    else
       PART_NUM = SLOW_PARTITION_LOOKUP(The_Set_Partition, item)
    end if
       

    return
  end function PARTITION_LOOKUP

  function SLOW_PARTITION_LOOKUP(The_Set_Partition, item) RESULT(PART_NUM)
    use pgslib_utility_module
    ! Returns the partition number of the partition which contains the item
    ! This version uses a search.
    ! If the partition cannot be found, PARTITION_INVALID_PARTITIONS is returned
    implicit none
    type(A_SET_PARTITIONING), intent(IN) :: The_Set_Partition
    integer(PGSLib_INT_TYPE), intent(IN) :: Item
    integer(PGSLib_Int_Type)             :: PART_NUM

    ! Local variables
    integer :: p, Number_Of_Partions
    integer(PGSLib_INT_TYPE), dimension(:),&
                              pointer     :: Start
    integer(PGSLib_INT_TYPE), dimension(:),&
                              pointer     :: End

    ! The simplest technique is by linear search.

    ! If the item is not in the set, then no point in searching
    ! Just return invalid partition.

    if ( (item <= 0) .OR. &
         (SIZE(The_Set_Partition%Set_That_Is_Partitioned) < item)) then
       PART_NUM = PARTITION_INVALID_PARTITIONS
       RETURN
    end if

    ! These arrays contain the first and last item of each partition
    Start => Get_START(The_Set_Partition)
    End   => Get_END(The_Set_Partition)
    Number_Of_Partions = Get_Num_Partitions(The_Set_Partition)

    PART_NUM = PARTITION_INVALID_PARTITIONS

    do p = 1, Number_Of_Partions
       if ( (Start(p) <= item) .AND. (item <= End(p) ) ) then
          PART_NUM = p
          exit
       end if
    end do

    ! If we didn''t find a partition, something went wrong
    if (PART_NUM == PARTITION_INVALID_PARTITIONS) then
       call PGSLib_Fatal_Error("Could not find partition which contains item in FIND_PARTITION")
    end if

    RETURN
  END function SLOW_PARTITION_LOOKUP

  function FAST_PARTITION_LOOKUP(The_Set_Partition, item) RESULT(PART_NUM)
    use pgslib_utility_module
    ! Returns the partition number of the partition which contains the item
    ! This version uses indexing into a large array.
    implicit none
    type(A_SET_PARTITIONING), intent(IN) :: The_Set_Partition
    integer(PGSLib_INT_TYPE), intent(IN) :: Item
    integer(PGSLib_Int_Type)             :: PART_NUM

    ! This routine does not work yet, so if we get here it is an error.

    ! First make sure that we can do fast lookup.
    ! If we are not initialized for it, that is an error.
    if (.NOT. The_Set_Partition%USE_FAST_LOOKUP) then
       call PGSLib_Fatal_Error("Must prepare for fast lookup by calling SETUP_FAST_LOOKUP before FAST_PARTITION_LOOKUP")
       return ! never do get here
    end if

    ! Now we can just do the indexing
!    PART_NUM = The_Set_Partition%FAST_LOOKUP_TABLE(Item)

    ! This routine does not work yet, so if we get here it is an error.
    ! this line makes the compiler happy
    PART_NUM =item
    PART_NUM = PARTITION_INVALID_PARTITIONS
    RETURN
  END function FAST_PARTITION_LOOKUP


  function PARTITION_IS_AVAILABLE(The_Set_Part, Partition) RESULT(Is_Avail)
    ! Returns TRUE is the partition is available to The_Set_Part object.
    implicit none
    type(A_SET_PARTITIONING),   intent(IN) :: The_Set_Part
    integer(PGSLib_Int_Type),   intent(IN) :: Partition
    logical(PGSLib_Log_Type)               :: Is_Avail
    

    ! Our test is whether the "processor" which holds the partition is
    ! the same one which holds the quereying object.

    Is_Avail = Get_Partition_ID(The_Set_Part) == LOOKUP(The_Set_Part%Map, Partition)

    RETURN
  END function PARTITION_IS_AVAILABLE

END MODULE PARTITION_DATA_TYPES
