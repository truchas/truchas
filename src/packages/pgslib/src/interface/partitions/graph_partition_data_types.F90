MODULE GRAPH_PARTITION_DATA_TYPES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE
  !   Provide the data types and support routines for paritioned
  !   graphs.
  !
  ! $Id: 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use Pgslib_Globals_Module
  use Partition_Constants
  use Partition_Data_Types

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: A_PARTITIONED_GRAPH
  PUBLIC :: INITIALIZE
  PUBLIC :: ALLOC
  PUBLIC :: FREE
  PUBLIC :: SET
  PUBLIC :: Get_Num_Edges_Available
  PUBLIC :: Get_Set_Partitioning
  PUBLIC :: Get_Head_Partitions
  PUBLIC :: Get_Tail_Partitions
  PUBLIC :: Get_Head_Local
  PUBLIC :: Get_Head_Available
  PUBLIC :: LOOKUP_HEAD
  PUBLIC :: LOOKUP_TAIL
  PUBLIC :: Head_Is_Local
  PUBLIC :: Head_Is_Available
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Constant declarations

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Type definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type A_PARTITIONED_GRAPH
     PRIVATE
     integer(PGSLib_Int_Type)                 :: Number_Available_Edges
     ! These are "global" edge numbers.  
     integer(PGSLib_Int_Type), dimension(:,:), &
                               POINTER        :: Edges
     ! These are the partitions for the heads and tails of the edges.
     ! I use two different arrays (rather than a 2-d array) since I often
     ! want to get at only the head or the tail.  In particular, I provide
     ! pointer access to each component individually, and having two arrays
     ! gives higher performance than pointing at an array section
     integer(PGSLib_Int_Type), dimension(:),   &
                               POINTER        :: Head_Partitions
     integer(PGSLib_Int_Type), dimension(:),   &
                               POINTER        :: Tail_Partitions
     ! This flag indicates whether this edge connects to vertices in the same partition.
     ! TRUE means the head is in the same partition as the tail
     logical(PGSLib_Log_Type), dimension(:),   &
                               POINTER        :: Head_Is_Local
     ! This flag indicates whether the edge connects to another partition which
     ! is available to the partitioned_graph object.  TRUE means that
     ! The head is available to this partitoned_graph object.
     logical(PGSLib_Log_Type), dimension(:),   &
                               POINTER        :: Head_Is_Available
     type(A_SET_PARTITIONING), POINTER        :: Head_Partition
     type(A_SET_PARTITIONING), POINTER        :: Tail_Partition
  end type A_PARTITIONED_GRAPH
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! INTERFACEs for generic procedure names           !!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE INITIALIZE
     MODULE PROCEDURE A_Part_Graph_INIT
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE A_Part_Graph_ALLOC
  END INTERFACE

  INTERFACE FREE
     MODULE PROCEDURE A_Part_Graph_FREE
  END INTERFACE

  INTERFACE SET
     MODULE PROCEDURE A_Part_Graph_SET
  END INTERFACE


  INTERFACE Get_Num_Edges_Available
     MODULE PROCEDURE Part_Graph_Get_Avail_Edges
  END INTERFACE
  
  INTERFACE Get_Set_Partitioning
     MODULE PROCEDURE Get_Set_Partitioning
  END INTERFACE

  INTERFACE Get_Head_Partitions
     MODULE PROCEDURE Get_Head_Partitions
  END INTERFACE

  INTERFACE Get_Tail_Partitions
     MODULE PROCEDURE Get_Tail_Partitions
  END INTERFACE

  INTERFACE Get_Head_Local
     MODULE PROCEDURE Get_Head_Local
  END INTERFACE

  INTERFACE Get_Head_Available
     MODULE PROCEDURE Get_Head_Available
  END INTERFACE

  INTERFACE LOOKUP_HEAD
     MODULE PROCEDURE LOOKUP_Head_Partition
  END INTERFACE

  INTERFACE LOOKUP_Tail
     MODULE PROCEDURE LOOKUP_Tail_Partition
  END INTERFACE
  
  INTERFACE Head_Is_Local
     MODULE PROCEDURE Head_Is_Local
  END INTERFACE

  INTERFACE Head_Is_Available
     MODULE PROCEDURE Head_Is_Available
  END INTERFACE

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routines for maintaining A_PARTITIONED_GRAPH data structures                     !
  ! Routines in this section (Generic name = specific name)
  !          INITIALIZE = A_Part_Graph_INIT
  !          ALLOC      = A_Part_Graph_ALLOC
  !          FREE       = A_Part_Graph_FREE
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine A_Part_Graph_INIT(The_Part_Graph)
    implicit none
    type(A_PARTITIONED_GRAPH), intent(OUT) :: The_Part_Graph
      
    The_Part_Graph%Number_Available_Edges = PARTITION_INVALID_SIZE

    NULLIFY(The_Part_Graph%Edges)
    NULLIFY(The_Part_Graph%Head_Partitions)
    NULLIFY(The_Part_Graph%Tail_Partitions)
    NULLIFY(The_Part_Graph%Head_Is_Local)
    NULLIFY(The_Part_Graph%Head_Is_Available)
    NULLIFY(The_Part_Graph%Head_Partition)
    NULLIFY(The_Part_Graph%Tail_Partition)

    RETURN
  END subroutine A_Part_Graph_INIT

  subroutine A_Part_Graph_Alloc(The_Part_Graph, Number_Available_Edges)
    implicit none
    type(A_PARTITIONED_GRAPH), intent(INOUT) :: The_Part_Graph
    integer,                   intent(IN   ) :: Number_Available_Edges

    The_Part_Graph%Number_Available_Edges = Number_Available_Edges
    ALLOCATE(The_Part_Graph%Edges(2, Number_Available_Edges))
    ALLOCATE(The_Part_Graph%Head_Partitions(Number_Available_Edges))
    ALLOCATE(The_Part_Graph%Tail_Partitions(Number_Available_Edges))
    ALLOCATE(The_Part_Graph%Head_Is_Local(Number_Available_Edges))
    ALLOCATE(The_Part_Graph%Head_Is_Available(Number_Available_Edges))

    RETURN
  END subroutine A_Part_Graph_ALLOC

  subroutine A_Part_Graph_Free(The_Part_Graph)
    implicit none
    type(A_PARTITIONED_GRAPH), intent(INOUT) :: The_Part_Graph

    DEALLOCATE(The_Part_Graph%Edges)
    DEALLOCATE(The_Part_Graph%Tail_Partitions)
    DEALLOCATE(The_Part_Graph%Head_Partitions)
    DEALLOCATE(The_Part_Graph%Head_Is_Local)
    DEALLOCATE(The_Part_Graph%Head_Is_Available)
    call INITIALIZE(The_Part_Graph)

    RETURN
  END subroutine A_Part_Graph_Free

  function Part_Graph_Get_Avail_Edges(The_Part_Graph) RESULT(Num_Edges)
    implicit none
    type(A_PARTITIONED_GRAPH), intent(IN   ) :: The_Part_Graph
    integer(PGSLib_Int_Type)                 :: Num_Edges

    Num_Edges = The_Part_Graph%Number_Available_Edges
    RETURN
  END function Part_Graph_Get_Avail_Edges    

  function Get_Set_Partitioning(The_Part_Graph, HeadOrTail) RESULT(A_Set_Part)
    implicit none
    type(A_PARTITIONED_GRAPH), TARGET   :: The_Part_Graph
    integer(PGSLib_Int_Type)            :: HeadOrTail
    type(A_SET_PARTITIONING),  POINTER  :: A_Set_Part

    if (HeadOrTail == GRAPH_HEAD) then
       A_Set_Part => The_Part_Graph%Head_Partition
    end if

    if (HeadOrTail == GRAPH_TAIL) then
       A_Set_Part => The_Part_Graph%Tail_Partition
    end if
    
    RETURN
  END function Get_Set_Partitioning

  function Get_Head_Partitions(The_Part_Graph) RESULT(Parts)
    implicit none
    type(A_PARTITIONED_GRAPH), TARGET       :: The_Part_Graph
    integer(PGSLib_Int_Type),  dimension(:), &
                               POINTER      :: Parts
    Parts => The_Part_Graph%Head_Partitions
    RETURN
  END function Get_Head_Partitions

  function Get_Tail_Partitions(The_Part_Graph) RESULT(Parts)
    implicit none
    type(A_PARTITIONED_GRAPH), TARGET       :: The_Part_Graph
    integer(PGSLib_Int_Type),  dimension(:), &
                               POINTER      :: Parts
    Parts => The_Part_Graph%Tail_Partitions
    RETURN
  END function Get_Tail_Partitions

  function Get_Head_Local(The_Part_Graph) RESULT(P)
    implicit none
    type(A_PARTITIONED_GRAPH), TARGET       :: The_Part_Graph
    logical(PGSLib_LOG_Type),  dimension(:), &
                               POINTER      :: P
    P => The_Part_Graph%Head_Is_Local
    RETURN
  END function Get_Head_Local

  function Get_Head_Available(The_Part_Graph) RESULT(P)
    implicit none
    type(A_PARTITIONED_GRAPH), TARGET       :: The_Part_Graph
    logical(PGSLib_LOG_Type),  dimension(:), &
                               POINTER      :: P
    P => The_Part_Graph%Head_Is_Available
    RETURN
  END function Get_Head_Available

  subroutine A_Part_Graph_SET(The_Part_Graph, Number_Available_Edges, &
                              Edges, Head_Partitioned_Set, Tail_Partitioned_Set)
    ! Set the fields in a partitioned graph
    implicit none
    type(A_PARTITIONED_GRAPH), intent(INOUT) :: The_Part_Graph
    integer(PGSLib_Int_Type),  intent(IN   ) :: Number_Available_Edges
    integer(PGSLib_Int_Type),  dimension(:,:),&
                               intent(IN   ) :: Edges
    type(A_SET_PARTITIONING),  TARGET        :: Head_Partitioned_Set
    type(A_SET_PARTITIONING),  TARGET        :: Tail_Partitioned_Set

    ! Local variables.
    integer :: e
    integer(PGSLib_Int_Type) :: Head_Partition
    integer(PGSLib_Int_Type) :: Tail_Partition

    ! First we have to allocate the partitioned graph
    call ALLOC(The_Part_Graph, Number_Available_Edges)

    ! Set the partitioning of the head and tail vertices
    The_Part_Graph%Head_Partition => Head_Partitioned_Set
    The_Part_Graph%Tail_Partition => Tail_Partitioned_Set

    ! Now put the edges in.  The edges must be with global vertex numbering.
    The_Part_Graph%Edges = Edges

    ! Finally, setup the partition numbers for the vertices of the edges
    do e = 1, Number_Available_Edges
       Head_Partition = LOOKUP(Head_Partitioned_Set, Edges(Graph_HEAD, e))
       Tail_Partition = LOOKUP(Tail_Partitioned_Set, Edges(Graph_TAIL, e))

       The_Part_Graph%Head_Partitions(e) = Head_Partition
       The_Part_Graph%Tail_Partitions(e) = Tail_Partition
       ! If the head is in a valid parition, stash some more info about that partition
       if (.NOT. (Head_Partition == PARTITION_INVALID_PARTITIONS)) then
          ! This flag is true if the head is in the same partition as the tail.
          The_Part_Graph%Head_Is_Local(e)                = (Head_Partition == Tail_Partition)
          ! This flag is true if the head is available to this partition
          The_Part_Graph%Head_Is_Available(e)            = Partition_Is_Available(Head_Partitioned_Set, &
                                                                                  Head_Partition)
       else
          The_Part_Graph%Head_Is_Local(e)                = .FALSE.
          The_Part_Graph%Head_Is_Available(e)            = .FALSE.
       end if
    end do

    return
  end subroutine A_Part_Graph_SET
    
  function LOOKUP_Tail_Partition(The_Part_Graph, Edge) RESULT(Part)
    implicit none
    type(A_PARTITIONED_GRAPH),intent(IN) :: The_Part_Graph
    integer(PGSLib_Int_Type), intent(IN) :: Edge
    integer(PGSLib_Int_Type)             :: Part

    ! Return the Tail partition of the given edge.

    Part = The_Part_Graph%Tail_Partitions(edge)
    return
  end function LOOKUP_Tail_Partition

  function LOOKUP_Head_Partition(The_Part_Graph, Edge) RESULT(Part)
    implicit none
    type(A_PARTITIONED_GRAPH),intent(IN) :: The_Part_Graph
    integer(PGSLib_Int_Type), intent(IN) :: Edge
    integer(PGSLib_Int_Type)             :: Part

    ! Return the Head partition of the given edge.

    Part = The_Part_Graph%Head_Partitions(edge)
    return
  end function LOOKUP_Head_Partition

  function Head_Is_Local(The_Part_Graph, Edge) RESULT(Is_Local)
    implicit none
    type(A_PARTITIONED_GRAPH),intent(IN) :: The_Part_Graph
    integer(PGSLib_Int_Type), intent(IN) :: Edge
    logical(PGSLib_Log_Type)             :: Is_Local

    ! Return .TRUE. if the head is in the same partition as the tail

    Is_Local = The_Part_Graph%Head_Is_Local(edge)
    return
  end function Head_Is_Local    

  function Head_Is_Available(The_Part_Graph, Edge) RESULT(Is_Available)
    implicit none
    type(A_PARTITIONED_GRAPH),intent(IN) :: The_Part_Graph
    integer(PGSLib_Int_Type), intent(IN) :: Edge
    logical(PGSLib_Log_Type)             :: Is_Available

    ! Return .TRUE. if the head is in the same partition as the tail

    Is_Available = The_Part_Graph%Head_Is_Available(edge)
    return
  end function Head_Is_Available    

END MODULE GRAPH_PARTITION_DATA_TYPES
