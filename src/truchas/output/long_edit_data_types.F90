#include "f90_assert.fpp"

MODULE LONG_EDIT_DATA_TYPES
  !=======================================================================
  ! Purpose(s):
  !
  !   Define and support the data types used for long edits.
  !   
  !   Data types provided
  !      LONG_EDIT_DATA
  !      LONG_EDIT_ITEM
  !
  !   Routines Provided
  !      COLLATE - for Long_Edit_Data types
  !
  ! Author(s): Robert C Ferrell (dbk@cpca.gov)
  !
  !  
  !=======================================================================
  use kinds, only: r8
  use matl_module,      only: MATERIAL, MATL_SLOT, COLLATE, PERMUTE_MATL
  use mesh_module,      only: CELL_GEOMETRY, COLLATE, PERMUTE_CELL
  use parameter_module, only: max_slots, mat_slot, nvc, ndim, nfc
  use zone_module,      only: CELL_AVG, Zone, COLLATE, PERMUTE_ZONE
  use solid_mechanics_data, only: CELL_MECH_INVARIANT
  use solid_mechanics_module, only: STRESS_STRAIN_INVARIANTS
  use truchas_logging_services
  implicit none
  public

  PUBLIC :: LONG_EDIT_LIST
  PUBLIC :: CREATE, DESTROY
  PUBLIC :: SET
  PUBLIC :: CELLNUMBER
  PUBLIC :: ZONE_DATA
  PUBLIC :: MATL_DATA
  PUBLIC :: MECH_EDIT_LIST
  PUBLIC :: CELL_GEO_DATA
  PUBLIC :: VERTEX_GEO_DATA
  PUBLIC :: VERTEX_DATA
  PUBLIC :: FACE_GEO_DATA
  PUBLIC :: FACE_DATA
  PUBLIC :: COLLATE
  PUBLIC :: SORT

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Data type for vertex information about a cell
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type VERTEX_DATA
     ! This is the vertex info for a particular neighbor vertex
     integer :: Ngbr_Vrtx
     real(r8), dimension(ndim) :: Vertex_Coord
  end type VERTEX_DATA

  type CELL_VERTEX_DATA
     ! This is the all the vertex info for a particular cell
     PRIVATE
     type(VERTEX_DATA), dimension(nvc) :: Vertices
  end type CELL_VERTEX_DATA

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Data type for face/cell neighbor information about a cell
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type FACE_DATA
     ! This is the face neigbhor information for a particular face
     integer :: Ngbr_Cell
     integer :: Ngbr_Face
     real(r8), dimension(ndim) :: Face_Normal
     real(r8), dimension(ndim) :: Face_Centroid
     real(r8) :: Face_Area
  end type FACE_DATA

  type CELL_NEIGHBOR_DATA
     ! This is the face/cell neighbor information for a cell
     type (FACE_DATA), dimension(nfc) :: Faces
  end type CELL_NEIGHBOR_DATA

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Data types to contain all information we want in a long edit
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type LONG_EDIT_PHYS_DATA
     ! This is the data which varies from cycle to cycle.
     PRIVATE
     type (CELL_AVG), POINTER, dimension(:) :: Zone => NULL()
     type (MATL_SLOT), dimension(max_slots) :: Matl
  end type LONG_EDIT_PHYS_DATA

  type LONG_EDIT_GEO_DATA
     ! This is the data which is geometric, and is constant from cycle to cycle
     PRIVATE
     type (CELL_GEOMETRY),     POINTER, dimension(:) :: Cell_Geo => NULL()
     type (CELL_VERTEX_DATA),  POINTER, dimension(:) :: Cell_Vrtx => NULL()
     type (CELL_NEIGHBOR_DATA),POINTER, dimension(:) :: Cell_Ngbr => NULL()
  end type LONG_EDIT_GEO_DATA

  type LONG_EDIT_LIST
     PRIVATE
     integer, dimension(:), pointer :: CellNumber  => NULL()! This is global cell number
     type (LONG_EDIT_PHYS_DATA)     :: Phys_Data            ! This is the data of the physical fields for that cell
     type (LONG_EDIT_GEO_DATA)      :: Geo_Data             ! This is the geometry of the cell
  end type LONG_EDIT_LIST
  
  type MECH_EDIT_LIST
     Private
     integer, dimension(:), pointer :: CellNumber => NULL() ! This is global cell number
     type (CELL_MECH_INVARIANT), dimension(:), pointer :: Mech_Data => NULL()
  end type MECH_EDIT_LIST

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  ! Interface statements
  INTERFACE CREATE
     MODULE PROCEDURE CREATE_LIST
     MODULE PROCEDURE CREATE_PHYS_DATA
     MODULE PROCEDURE CREATE_GEO_DATA
     MODULE PROCEDURE CREATE_MECH_LIST
     MODULE PROCEDURE CREATE_MECH_DATA
  END INTERFACE

  INTERFACE CLONE
     MODULE PROCEDURE CLONE_LIST
     MODULE PROCEDURE CLONE_PHYS_DATA
     MODULE PROCEDURE CLONE_GEO_DATA
     MODULE PROCEDURE CLONE_MECH_LIST
     MODULE PROCEDURE CLONE_MECH_DATA
  END INTERFACE

  INTERFACE DESTROY
     MODULE PROCEDURE DESTROY_LIST
     MODULE PROCEDURE DESTROY_PHYS_DATA
     MODULE PROCEDURE DESTROY_GEO_DATA
     MODULE PROCEDURE DESTROY_MECH_LIST
     MODULE PROCEDURE DESTROY_MECH_DATA
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE ASSIGN_LIST
     MODULE PROCEDURE ASSIGN_PHYS_DATA
     MODULE PROCEDURE ASSIGN_GEO_DATA
!     MODULE PROCEDURE ASSIGN_MECH_LIST
!     MODULE PROCEDURE ASSIGN_MECH_DATA
  END INTERFACE

  INTERFACE SIZE
     MODULE PROCEDURE SIZE_LIST
     MODULE PROCEDURE SIZE_PHYS_DATA
     MODULE PROCEDURE SIZE_GEO_DATA
     MODULE PROCEDURE SIZE_MECH_LIST
  END INTERFACE

  INTERFACE SET
     MODULE PROCEDURE SET_LIST
     MODULE PROCEDURE SET_PHYS_DATA
     MODULE PROCEDURE SET_GEO_DATA
     MODULE PROCEDURE SET_VERTEX_DATA
     MODULE PROCEDURE SET_FACE_NEIGHBOR_DATA
     MODULE PROCEDURE SET_MECH_LIST
     MODULE PROCEDURE SET_MECH_DATA
  END INTERFACE

  INTERFACE CELLNUMBER
     MODULE PROCEDURE LIST_CELLNUMBER
     MODULE PROCEDURE MECH_LIST_CELLNUMBER
  END INTERFACE

  INTERFACE PHYS_DATA
     MODULE PROCEDURE LIST_PHYS_DATA
  END INTERFACE

  INTERFACE VERTEX_GEO_DATA
     MODULE PROCEDURE vertex_data_vertex_geo_number
     MODULE PROCEDURE geo_data_item_vnumber_vertex
     MODULE PROCEDURE geo_data_cell_vertex_data
     MODULE PROCEDURE geo_data_item_cell_vertex_data
     MODULE PROCEDURE list_item_cell_vertex_data
     MODULE PROCEDURE list_item_cell_vnumber_vertex
  END INTERFACE

  INTERFACE NGBR_GEO_DATA
     MODULE PROCEDURE geo_data_item_cell_ngbr_data
     MODULE PROCEDURE geo_data_cell_ngbr_data
     MODULE PROCEDURE list_item_cell_ngbr_data
  END INTERFACE

  INTERFACE FACE_GEO_DATA
     MODULE PROCEDURE get_ngbr_data_face_number_face
     MODULE PROCEDURE geo_data_item_fnumber_face
     MODULE PROCEDURE list_item_cell_fnumber_face
  END INTERFACE

  INTERFACE CELL_GEO_DATA
     MODULE PROCEDURE geo_data_item_cell_geo
     MODULE PROCEDURE geo_data_cell_geo
     MODULE PROCEDURE LIST_ITEM_CELL_GEO
  END INTERFACE

  INTERFACE GEO_DATA
     MODULE PROCEDURE list_geo_data
  END INTERFACE

  INTERFACE ZONE_DATA
     MODULE PROCEDURE phys_data_item_zone
     MODULE PROCEDURE phys_data_zone
     MODULE PROCEDURE LIST_ZONE
     MODULE PROCEDURE LIST_ITEM_ZONE
  END INTERFACE

  INTERFACE MATL_DATA
     MODULE PROCEDURE phys_data_item_matl
     MODULE PROCEDURE phys_data_slot_matl
     MODULE PROCEDURE phys_data_matl
     MODULE PROCEDURE LIST_MATL_SLOT_MATL
     MODULE PROCEDURE LIST_ITEM_MATL
     MODULE PROCEDURE LIST_MATL
  END INTERFACE

  INTERFACE COLLATE
     MODULE PROCEDURE LIST_COLLATE
     MODULE PROCEDURE PHYS_DATA_COLLATE
     MODULE PROCEDURE GEO_DATA_COLLATE
  END INTERFACE

  INTERFACE COLLATE
     MODULE PROCEDURE collate_vertex_list
     MODULE PROCEDURE collate_ngbr_list
     MODULE PROCEDURE mech_list_collate
  END INTERFACE

  INTERFACE SORT
     MODULE PROCEDURE LIST_SORT
     MODULE PROCEDURE PHYS_DATA_SORT
     MODULE PROCEDURE GEO_DATA_SORT
     MODULE PROCEDURE MECH_LIST_SORT
  END INTERFACE

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Support routines for CELL_VERTEX_DATA and VERTEX_DATA
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SET_VERTEX_DATA(CELL_VERTEX, VERTEX_NUMBER, VERTEX)
    type (CELL_VERTEX_DATA), target, intent(INOUT) :: CELL_VERTEX
    integer, intent(IN) :: VERTEX_NUMBER
    type (VERTEX_DATA),target, intent(IN) :: VERTEX

    CELL_VERTEX%Vertices(VERTEX_NUMBER) = VERTEX
  end subroutine SET_VERTEX_DATA

  function vertex_data_vertex_geo_number(CELL_VERTEX, VERTEX_NUMBER) RESULT(VERTEX)
    type (CELL_VERTEX_DATA), intent(IN), &
                             TARGET        :: CELL_VERTEX
    integer, intent(IN) :: VERTEX_NUMBER
    type(VERTEX_DATA), POINTER :: VERTEX

    VERTEX => CELL_VERTEX%Vertices(VERTEX_NUMBER)
  end function vertex_data_vertex_geo_number

  subroutine collate_vertex_list(collated_list, local_list)
    use pgslib_module
    type(CELL_VERTEX_DATA), intent(INOUT), &
                           TARGET,        &
                           dimension(:)  :: collated_list
    type(CELL_VERTEX_DATA), intent(IN), &
                           TARGET,        &
                           dimension(:)  :: local_list

    ! Local variables
    integer :: v, d
    type(VERTEX_DATA), POINTER, dimension(:) :: collated_vrtx_data => NULL()
    type(VERTEX_DATA), POINTER, dimension(:) :: local_vrtx_data => NULL()
    integer,      POINTER,      &
                            dimension(:) :: collated_int_temp => NULL()
    integer,      POINTER,      &
                            dimension(:) :: local_int_temp => NULL()
    real(r8),        POINTER,      &
                            dimension(:) :: collated_real_temp => NULL()
    real(r8),        POINTER,      &
                            dimension(:) :: local_real_temp => NULL()

    do v = 1, nvc
       local_vrtx_data => local_list%Vertices(v)
       collated_vrtx_data => collated_list%Vertices(v)
       
       ! Collate the ngbr vertex
       local_int_temp    => local_vrtx_data%Ngbr_Vrtx
       collated_int_temp => collated_vrtx_data%Ngbr_Vrtx
       call pgslib_collate(collated_int_temp, local_int_temp)

       ! Collate the coordinates
       do d = 1, ndim
          local_real_temp => local_vrtx_data%Vertex_Coord(d)
          collated_real_temp => collated_vrtx_data%Vertex_Coord(d)
          call pgslib_collate(collated_real_temp, local_real_temp)
       end do

    end do
  end subroutine collate_vertex_list
       
  subroutine permute_vertex_list(permuted_list, orig_list, PERMUTER, SCOPE)
    use parallel_scope
    use pgslib_module
    type(CELL_VERTEX_DATA), intent(INOUT), &
                           TARGET,        &
                           dimension(:)  :: permuted_list
    type(CELL_VERTEX_DATA), intent(IN), &
                           TARGET,        &
                           dimension(:)  :: orig_list
    integer,     intent(IN), &
                           dimension(:)  :: PERMUTER
    type (PL_SCOPE),       intent(IN) :: SCOPE

    ! Local variables
    integer :: c

    ! This routine works only for local scope so far
    if (SCOPE /= LOCAL_SCOPE) then
       call TLS_fatal ('PERMUTE_VERTEX_LIST: only works with local scope')
    end if

    do c = 1, SIZE(Permuter)
       permuted_list(permuter(c)) = orig_list(c)
    end do

  end subroutine permute_vertex_list
       

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Support routines for FACE_NEIGHBOR_DATA and CELL_NEIGHBOR_DATA
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  subroutine SET_FACE_NEIGHBOR_DATA(CELL_NGBR, FACE_NUMBER, FACE)
    type (CELL_NEIGHBOR_DATA), intent(INOUT) :: CELL_NGBR
    integer,        intent(IN) :: FACE_NUMBER
    type (FACE_DATA),          intent(IN) :: FACE
                            

    CELL_NGBR%Faces(FACE_NUMBER) = FACE
  end subroutine SET_FACE_NEIGHBOR_DATA

  function get_ngbr_data_face_number_face(CELL_NGBR, FACE_NUMBER) RESULT(FACE)
    type (CELL_NEIGHBOR_DATA), intent(IN), &
                               TARGET        :: CELL_NGBR
    integer,         intent(IN) :: FACE_NUMBER
    type(FACE_DATA),           POINTER       :: FACE

    FACE => CELL_NGBR%Faces(Face_NUMBER)
  end function get_ngbr_data_face_number_face

  subroutine collate_ngbr_list(collated_list, local_list)
    use pgslib_module
    type(CELL_NEIGHBOR_DATA), intent(INOUT), &
                           TARGET,        &
                           dimension(:)  :: collated_list
    type(CELL_NEIGHBOR_DATA), intent(IN), &
                           TARGET,        &
                           dimension(:)  :: local_list

    ! Local variables
    integer :: f, d
    type(FACE_DATA), POINTER,      &
                            dimension(:) :: collated_face_data => NULL()
    type(FACE_DATA), POINTER,      &
                            dimension(:) :: local_face_data => NULL()
    integer,      POINTER,      &
                            dimension(:) :: collated_int_temp => NULL()
    integer,      POINTER,      &
                            dimension(:) :: local_int_temp => NULL()
    real(r8),        POINTER,      &
                            dimension(:) :: collated_real_temp => NULL()
    real(r8),        POINTER,      &
                            dimension(:) :: local_real_temp => NULL()

    do f = 1, nfc
       local_face_data => local_list%Faces(f)
       collated_face_data => collated_list%Faces(f)
       
       ! Collate the ngbr cell
       local_int_temp    => local_face_data%Ngbr_Cell
       collated_int_temp => collated_face_data%Ngbr_Cell
       call pgslib_collate(collated_int_temp, local_int_temp)

       ! Collate the ngbr face
       local_int_temp    => local_face_data%Ngbr_face
       collated_int_temp => collated_face_data%Ngbr_face
       call pgslib_collate(collated_int_temp, local_int_temp)
       
       ! Collate the face normal coordinates
       do d = 1, ndim
          local_real_temp => local_face_data%Face_Normal(d)
          collated_real_temp => collated_face_data%Face_Normal(d)
          call pgslib_collate(collated_real_temp, local_real_temp)
       end do

       ! Collate the face centroid coordinates
       do d = 1, ndim
          local_real_temp => local_face_data%Face_Centroid(d)
          collated_real_temp => collated_face_data%Face_Centroid(d)
          call pgslib_collate(collated_real_temp, local_real_temp)
       end do

       ! Collate the face area
       local_real_temp => local_face_data%Face_area
       collated_real_temp => collated_face_data%Face_area
       call pgslib_collate(collated_real_temp, local_real_temp)

    end do
  end subroutine collate_ngbr_list
       
  subroutine permute_ngbr_list(permuted_list, orig_list, PERMUTER, SCOPE)
    use parallel_scope
    use pgslib_module
    type(CELL_NEIGHBOR_DATA), intent(INOUT), &
                           TARGET,        &
                           dimension(:)  :: permuted_list
    type(CELL_NEIGHBOR_DATA), intent(IN), &
                           TARGET,        &
                           dimension(:)  :: orig_list
    integer,     intent(IN), &
                           dimension(:)  :: PERMUTER
    type (PL_SCOPE),       intent(IN) :: SCOPE

    ! Local variables
    integer :: c

    ! This routine works only for local scope so far
    if (SCOPE /= LOCAL_SCOPE) then
       call TLS_fatal ('PERMUTE_VERTEX_LIST: only works with local scope')
    end if

    do c = 1, SIZE(Permuter)
       permuted_list(permuter(c)) = orig_list(c)
    end do

  end subroutine permute_ngbr_list
       

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Support routines for LONG_EDIT_PHYS_DATA type !!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine CREATE_PHYS_DATA(DATA, ITEMS)
    ! Allocate the internal storage for the data
    type (LONG_EDIT_PHYS_DATA), intent(INOUT) :: DATA
    integer,     intent(IN) :: ITEMS

    ! local variables
    integer :: s

    ! we could save a little memory by only allocating the needed slot pointers,
    ! but then we might have to allocate some slot pointers during the run
    do s = 1, SIZE(Data%Matl)
       ! Allocate space for all the slots we will use
       ALLOCATE(Data%Matl(s)%Cell(ITEMS))
    end do
    
    ALLOCATE(Data%Zone(ITEMS))
  end subroutine CREATE_PHYS_DATA

  subroutine CLONE_PHYS_DATA(NEW, OLD)
    ! Clone a long edit item, the new one is empty
    type (LONG_EDIT_PHYS_DATA), intent(OUT) :: NEW
    type (LONG_EDIT_PHYS_DATA), intent(IN) :: OLD

    call CREATE(DATA=NEW, ITEMS=SIZE(OLD))
  end subroutine CLONE_PHYS_DATA
  
  subroutine DESTROY_PHYS_DATA(DATA)
    ! DeAllocate the internal storage for the data
    type (LONG_EDIT_PHYS_DATA), intent(INOUT) :: DATA

    ! local variables
    integer :: s

    do s = 1, SIZE(Data%Matl)
       ! DeAllocate space for all the slots we used
       DEALLOCATE(Data%Matl(s)%Cell)
    end do
    
    DEALLOCATE(Data%Zone)
  end subroutine DESTROY_PHYS_DATA

  subroutine ASSIGN_PHYS_DATA(LDATA, RDATA)
    type (LONG_EDIT_PHYS_DATA), intent(INOUT) :: LDATA
    type (LONG_EDIT_PHYS_DATA), intent(IN) :: RDATA

    ! Local variables
    integer :: s

    LDATA%Zone = RData%Zone

    do s = 1, mat_slot
       LData%Matl(s)%Cell = RData%Matl(s)%Cell
    end do
  end subroutine ASSIGN_PHYS_DATA

  function SIZE_PHYS_DATA(DATA) RESULT(S)
    type (LONG_EDIT_PHYS_DATA), intent(IN) :: DATA
    integer                :: S
    ! We assume that the create worked properly so we need only look
    ! at size of one of the elements of DATA
    S = SIZE(DATA%Zone)
  end function SIZE_PHYS_DATA
    
  subroutine SET_PHYS_DATA(DATA, ITEM, ZONE, MATL_SLOT, MATL)
    ! Put the values into the long edit item
    type (LONG_EDIT_PHYS_DATA), intent(INOUT) :: DATA
    integer,     intent(IN) :: ITEM
    type (CELL_AVG),       OPTIONAL,      &
                           intent(IN) :: ZONE
    integer,     OPTIONAL,      &
                           intent(IN) :: MATL_SLOT
    type (MATERIAL),       OPTIONAL,      &
                           intent(IN) :: MATL

    ! Local variables
    
    if (PRESENT(ZONE)) then
       Data%Zone(Item) = Zone
    end if
    
    if (PRESENT(MATL)) then
       ASSERT(PRESENT(MATL_SLOT))
       Data%Matl(MATL_SLOT)%Cell(Item) = MATL
    end if

  end subroutine SET_PHYS_DATA


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Support routines for Mech data (Solid Mechanics) !!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine CREATE_MECH_DATA(DATA, ITEMS)
    ! Allocate the internal storage for the data
    type (CELL_MECH_INVARIANT), pointer, dimension(:) :: DATA
    integer,     intent(IN)              :: ITEMS

    ! local variables

    ALLOCATE(Data(ITEMS))
  end subroutine CREATE_MECH_DATA

  subroutine CLONE_MECH_DATA(NEW, OLD)
    ! Clone a long edit mech item, the new one is empty
    type (CELL_MECH_INVARIANT), pointer, dimension(:) :: NEW
    type (CELL_MECH_INVARIANT), pointer, dimension(:) :: OLD

    call CREATE(DATA=NEW, ITEMS=SIZE(OLD))
  end subroutine CLONE_MECH_DATA
  
  subroutine DESTROY_MECH_DATA(DATA)
    ! DeAllocate the internal storage for the data
    type (CELL_MECH_INVARIANT), pointer, dimension(:) :: DATA

    DEALLOCATE(Data)
  end subroutine DESTROY_MECH_DATA

!  subroutine ASSIGN_MECH_DATA(LDATA, RDATA)
!    implicit none
!    type (CELL_MECH_INVARIANT), pointer, intent(IN)  :: LDATA
!    type (CELL_MECH_INVARIANT), pointer, intent(OUT) :: RDATA
!
!    LDATA%mises_stress = RData%mises_stress
!    LDATA%eff_plastic_strain = RDATA%eff_plastic_strain
!    LDATA%mean_stress = RDATA%mean_stress
!    LDATA%volumetric_strain = RDATA%volumetric_strain
!
!    return
!  end subroutine ASSIGN_MECH_DATA

  function SIZE_MECH_DATA(DATA) RESULT(S)
    type (CELL_MECH_INVARIANT), pointer, dimension(:) :: DATA
    integer                :: S
    ! We assume that the create worked properly so we need only look
    ! at size of one of the elements of DATA
    S = SIZE(DATA)
  end function SIZE_MECH_DATA
    
  subroutine SET_MECH_DATA(DATA, ITEM, MECH)
    ! Put the values into the long edit item
    type (CELL_MECH_INVARIANT), pointer, dimension(:) :: DATA
    type (CELL_MECH_INVARIANT), pointer               :: Mech
    integer,     intent(IN) :: ITEM

    Data(ITEM) = Mech

  end subroutine SET_MECH_DATA

  !!!!!!!!!! Support routines for LONG_EDIT_GEO_DATA type !!!!!!!!!!

  subroutine CREATE_GEO_DATA(DATA, ITEMS)
    ! Allocate the internal storage for the data
    type (LONG_EDIT_GEO_DATA), intent(INOUT) :: DATA
    integer,         intent(IN) :: ITEMS

    ALLOCATE(Data%Cell_Geo(ITEMS))
    ALLOCATE(Data%Cell_Vrtx(ITEMS))
    ALLOCATE(Data%Cell_Ngbr(ITEMS))
  end subroutine CREATE_GEO_DATA

  subroutine CLONE_GEO_DATA(NEW, OLD)
    ! Clone a long edit item, the new one is empty
    type (LONG_EDIT_GEO_DATA), intent(OUT) :: NEW
    type (LONG_EDIT_GEO_DATA), intent(IN) :: OLD

    call CREATE(DATA=NEW, ITEMS=SIZE(OLD))
  end subroutine CLONE_GEO_DATA
  
  subroutine DESTROY_GEO_DATA(DATA)
    ! DeAllocate the internal storage for the data
    type (LONG_EDIT_GEO_DATA), intent(INOUT) :: DATA

    DEALLOCATE(Data%Cell_Geo)
    DEALLOCATE(Data%Cell_Vrtx)
    DEALLOCATE(Data%Cell_Ngbr)
  end subroutine DESTROY_GEO_DATA

  subroutine ASSIGN_GEO_DATA(LDATA, RDATA)
    type (LONG_EDIT_GEO_DATA), intent(INOUT) :: LDATA
    type (LONG_EDIT_GEO_DATA), intent(IN) :: RDATA

    LDATA%Cell_Geo  = RData%Cell_Geo
    LDATA%Cell_Vrtx = RData%Cell_Vrtx
    LDATA%Cell_Ngbr = RData%Cell_Ngbr
  end subroutine ASSIGN_GEO_DATA

  function SIZE_GEO_DATA(DATA) RESULT(S)
    type (LONG_EDIT_GEO_DATA), intent(IN) :: DATA
    integer                     :: S

    S = SIZE(DATA%Cell_Geo)

  end function SIZE_GEO_DATA
    
  subroutine SET_GEO_DATA(DATA, ITEM, CELL_GEO, &
                          VERTEX_NUMBER, VERTEX,&
                          FACE_NUMBER, FACE)
    ! Put the values into the long edit item
    type (LONG_EDIT_GEO_DATA), intent(INOUT) :: DATA
    integer,         intent(IN) :: ITEM
    type (CELL_GEOMETRY),      OPTIONAL,      &
                               intent(IN) :: CELL_GEO
    integer,         OPTIONAL,      &
                               intent(IN) :: VERTEX_NUMBER
    type (VERTEX_DATA), OPTIONAL, intent(IN) :: VERTEX

    integer,         OPTIONAL,      &
                               intent(IN) :: FACE_NUMBER
    type (FACE_DATA),   OPTIONAL,      &
                               intent(IN) :: FACE
    ! Local variables

    IF (PRESENT(CELL_GEO))    Data%Cell_Geo(Item)  = CELL_GEO
    
    ! If either VERTEX_NUMBER of CELL_VERTEX are provided, both must be provided
    ASSERT(PRESENT(VERTEX_NUMBER) .EQV. PRESENT(VERTEX))
    if (PRESENT(VERTEX_NUMBER)) then
       call SET(Data%Cell_Vrtx(Item), VERTEX_NUMBER = VERTEX_NUMBER, VERTEX = VERTEX)
    end IF

    ! If either FACE_NUMBER of CELL_FACE are provided, both must be provided
    ASSERT(PRESENT(FACE_NUMBER) .EQV. PRESENT(FACE))
    if (PRESENT(FACE_NUMBER)) then
       call SET(Data%Cell_Ngbr(Item), FACE_NUMBER = FACE_NUMBER, FACE = FACE)
    end IF

  end subroutine SET_GEO_DATA

  !!!!!!!!!! Support routines for LONG_EDIT_LIST type !!!!!!!!!!

  subroutine CREATE_LIST(LIST, ITEMS)
    ! Allocate the internal storage for the item list
    type (LONG_EDIT_LIST), intent(INOUT) :: LIST
    integer,     intent(IN) :: ITEMS

    ! Allocate space for the cell number
    ALLOCATE(List%CellNumber(ITEMS))

    ! Now allocate the data
    call CREATE(LIST%Phys_Data, ITEMS)
    call CREATE(LIST%Geo_Data, ITEMS)
  end subroutine CREATE_LIST

  subroutine CLONE_LIST(NEW, OLD)
    ! Clone a long edit item, the new one is empty
    type (LONG_EDIT_LIST), intent(OUT) :: NEW
    type (LONG_EDIT_LIST), intent(IN) :: OLD

    call CREATE(LIST=NEW, ITEMS=SIZE(OLD))
  end subroutine CLONE_LIST
  
  subroutine DESTROY_LIST(LIST)
    ! DeAllocate the internal storage for the item list
    type (LONG_EDIT_LIST), intent(INOUT) :: LIST

    ! Allocate space for the cell number
    DEALLOCATE(List%CellNumber)

    ! Now deallocate the data
    call DESTROY(LIST%Phys_Data)
    call DESTROY(LIST%Geo_Data)
  end subroutine DESTROY_LIST

  subroutine ASSIGN_LIST(LITEM, RITEM)
    type (LONG_EDIT_LIST), intent(INOUT) :: LITEM
    type (LONG_EDIT_LIST), intent(IN) :: RITEM

    LITEM%CellNumber = RITEM%CellNumber
    LITEM%Phys_Data  = RITEM%Phys_Data
    LITEM%Geo_Data   = RITEM%Geo_Data

  end subroutine ASSIGN_LIST

  function SIZE_LIST(LIST) RESULT(S)
    type (LONG_EDIT_LIST), intent(IN) :: LIST
    integer                 :: S
    ! We assume that the create worked properly so we need only look
    ! at size of one of the elements of LIST
    S = SIZE(List%CellNumber)
  end function SIZE_LIST
    
  subroutine SET_LIST(LIST, ITEM, CELLNUMBER, ZONE, MATL_SLOT, MATL, &
                      CELL_GEO,                                      &
                      VERTEX_NUMBER, VERTEX,                         &
                      FACE_NUMBER,   FACE)
    ! Put the values into the long edit item
    type (LONG_EDIT_LIST), intent(INOUT) :: LIST
    integer,     intent(IN) :: ITEM
    integer,     OPTIONAL,      &
                           intent(IN) :: CELLNUMBER
    type (CELL_AVG),       OPTIONAL,      &
                           intent(IN) :: ZONE
    integer,     OPTIONAL,      &
                           intent(IN) :: MATL_SLOT
    type (MATERIAL),       OPTIONAL,      &
                           intent(IN) :: MATL
    type (CELL_GEOMETRY),  OPTIONAL,      &
                           intent(IN) :: CELL_GEO
    integer,     OPTIONAL,      &
                           intent(IN) :: VERTEX_NUMBER
    type(VERTEX_DATA),OPTIONAL, intent(IN) :: VERTEX
    integer,     OPTIONAL,      &
                           intent(IN) :: FACE_NUMBER
    type(FACE_DATA),OPTIONAL,      &
                           intent(IN) :: FACE
    ! Local variables

    if (PRESENT(CellNumber)) then
       LIST%CellNumber(Item) = CELLNUMBER
    end if
    if (PRESENT(ZONE) .OR. PRESENT(MATL)) then
       call SET(LIST%Phys_DATA, ITEM, ZONE, MATL_SLOT, MATL)
    end if
    ! If either VERTEX_NUMBER of VERTEX are provided, both must be provided
    ASSERT(PRESENT(VERTEX_NUMBER) .EQV. PRESENT(VERTEX))

    if (PRESENT(VERTEX_NUMBER)) then
       call SET(LIST%Geo_Data, ITEM, &
                VERTEX_NUMBER = VERTEX_NUMBER, VERTEX = VERTEX)
    endif

    ! If either FACE_NUMBER of FACE are provided, both must be provided
    ASSERT(PRESENT(FACE_NUMBER) .EQV. PRESENT(FACE))

    if (PRESENT(FACE_NUMBER)) then
       call SET(LIST%Geo_Data, ITEM, &
                FACE_NUMBER = FACE_NUMBER, FACE = FACE)
    endif

    if (PRESENT(Cell_GEO)) then
       call SET(LIST%Geo_Data, ITEM, CELL_GEO = CELL_GEO)
    end if

  end subroutine SET_LIST
  

  !!!!!!!!!! Support routines for MECH_EDIT_LIST type !!!!!!!!!!

  subroutine CREATE_MECH_LIST(LIST, ITEMS)
    ! Allocate the internal storage for the item list
    type (MECH_EDIT_LIST), intent(INOUT) :: LIST
    integer,     intent(IN) :: ITEMS

    ! Allocate space for the cell number
    ALLOCATE(List%CellNumber(ITEMS))

    ! Now allocate the data
    call CREATE(LIST%Mech_Data, ITEMS)
  end subroutine CREATE_MECH_LIST

  subroutine CLONE_MECH_LIST(NEW, OLD)
    ! Clone a long edit item, the new one is empty
    type (MECH_EDIT_LIST), intent(OUT) :: NEW
    type (MECH_EDIT_LIST), intent(IN)  :: OLD

    call CREATE(LIST=NEW, ITEMS=SIZE(OLD))
  end subroutine CLONE_MECH_LIST
  
  subroutine DESTROY_MECH_LIST(LIST)
    ! DeAllocate the internal storage for the item list
    type (MECH_EDIT_LIST), intent(INOUT) :: LIST

    ! Deallocate space for the cell number
    DEALLOCATE(List%CellNumber)

    ! Now deallocate the data
    call DESTROY(LIST%MECH_Data)
  end subroutine DESTROY_MECH_LIST

!  subroutine ASSIGN_MECH_LIST(LITEM, RITEM)
!    implicit none
!    type (MECH_EDIT_LIST), intent(OUT) :: LITEM
!    type (MECH_EDIT_LIST), intent(IN) :: RITEM
!
!    ! Local variables
!    integer :: s
!
!    LITEM%CellNumber = RITEM%CellNumber
!    LITEM%Mech_Data  = RITEM%Phys_Data
!
!    return
!  end subroutine ASSIGN_MECH_LIST

  function SIZE_MECH_LIST(LIST) RESULT(S)
    type (MECH_EDIT_LIST), intent(IN) :: LIST
    integer                 :: S
    ! We assume that the create worked properly so we need only look
    ! at size of one of the elements of LIST
    S = SIZE(List%CellNumber)
  end function SIZE_MECH_LIST
    
  subroutine SET_MECH_LIST(LIST, ITEM, CELLNUMBER, DATA)
    ! Put the values into the long edit item
    type (MECH_EDIT_LIST), intent(INOUT) :: LIST
    integer,     intent(IN) :: ITEM
    integer,     intent(IN) :: CELLNUMBER
    type (CELL_MECH_INVARIANT), pointer :: DATA
    ! Local variables

    LIST%CellNumber(Item) = CELLNUMBER

    call SET(LIST%Mech_DATA, ITEM, DATA)

  end subroutine SET_MECH_LIST

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  ACCESSOR FUNCTIONS AND ROUTINES                                 !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function phys_data_item_zone(data, item) RESULT(zone)
    type (LONG_EDIT_PHYS_DATA), TARGET,     &
                           intent(IN) :: DATA
    integer,     intent(IN) :: item
    type (CELL_AVG), POINTER          :: Zone

    Zone => data%Zone(item)
  end function phys_data_item_zone
  
  function phys_data_zone(data) RESULT(zone)
    type (LONG_EDIT_PHYS_DATA), TARGET,     &
                           intent(IN) :: DATA
    type (CELL_AVG),       dimension(:),&
                           POINTER    :: Zone

    Zone => data%Zone
  end function phys_data_zone
  
  function phys_data_item_matl(data, item, matl_slot) RESULT(matl)
    type (LONG_EDIT_PHYS_DATA), TARGET,     &
                           intent(IN) :: DATA
    integer,     intent(IN) :: item
    integer,     intent(IN) :: matl_slot
    type (MATERIAL),       POINTER    :: Matl

    Matl => data%matl(matl_slot)%Cell(item)
  end function phys_data_item_matl

  function phys_data_slot_matl(data, matl_slot) RESULT(matl)
    type (LONG_EDIT_PHYS_DATA), TARGET,     &
                           intent(IN) :: DATA
    integer,     intent(IN) :: matl_slot
    type (MATERIAL),       dimension(:),&
                           POINTER    :: Matl

    Matl => data%matl(matl_slot)%Cell
  end function phys_data_slot_matl

  function phys_data_matl(data) RESULT(matl)
    type (LONG_EDIT_PHYS_DATA), TARGET,     &
                           intent(IN) :: DATA
    type (MATL_SLOT),      dimension(:),&
                           POINTER    :: Matl

    Matl => data%matl
  end function phys_data_matl

  !!!!!!!!!! Accessor routines for LONG_EDIT_GEO_DATA type !!!!!!!!!!

  function geo_data_item_cell_geo(data, item) RESULT(cell)
    type (LONG_EDIT_GEO_DATA), TARGET,     &
                             intent(IN) :: DATA
    integer,       intent(IN) :: item
    type (CELL_GEOMETRY),    POINTER    :: Cell

    Cell => data%cell_geo(item)
  end function geo_data_item_cell_geo
  
  function geo_data_cell_geo(data) RESULT(cell)
    type (LONG_EDIT_GEO_DATA), TARGET,     &
                             intent(IN) :: DATA
    type (CELL_GEOMETRY),    dimension(:), &
                             POINTER    :: Cell

    Cell => data%cell_geo
  end function geo_data_cell_geo
  
  function geo_data_cell_vertex_data(geo_data) RESULT(cell_vrtx_geo)
    type (LONG_EDIT_GEO_DATA), TARGET,     &
                             intent(IN) :: GEO_DATA
    type (CELL_VERTEX_DATA),  dimension(:), &
                             POINTER    :: cell_vrtx_geo

    cell_vrtx_geo => geo_data%cell_vrtx
  end function geo_data_cell_vertex_data
  
  function geo_data_item_cell_vertex_data(geo_data, item) RESULT(cell_vrtx_geo)
    type (LONG_EDIT_GEO_DATA), TARGET,     &
                             intent(IN) :: GEO_DATA
    integer,       intent(IN) :: item
    type (CELL_VERTEX_DATA), POINTER    :: cell_vrtx_geo

    cell_vrtx_geo => geo_data%cell_vrtx(item)
  end function geo_data_item_cell_vertex_data

  function geo_data_item_vnumber_vertex(geo_data, item, vertex_number) RESULT(vertex)
    type (LONG_EDIT_GEO_DATA), TARGET,     &
                             intent(IN) :: GEO_DATA
    integer,       intent(IN) :: item
    integer,       intent(IN) :: vertex_number
    type (VERTEX_DATA),      POINTER    :: Vertex

    ! Local variables
    type (CELL_VERTEX_DATA), POINTER :: CELL_VERTEX => NULL()
    
    Cell_Vertex => VERTEX_GEO_DATA(geo_data, item)
    Vertex      => VERTEX_GEO_DATA(cell_vertex, vertex_number)
  end function geo_data_item_vnumber_vertex

  function geo_data_cell_ngbr_data(geo_data) RESULT(cell_ngbr)
    type (LONG_EDIT_GEO_DATA), TARGET,     &
                             intent(IN)   :: GEO_DATA
    type (CELL_NEIGHBOR_DATA),dimension(:),&
                              POINTER     :: cell_ngbr

    cell_ngbr => geo_data%cell_Ngbr
  end function geo_data_cell_ngbr_data
  
  function geo_data_item_cell_ngbr_data(geo_data, item) RESULT(cell_ngbr)
    type (LONG_EDIT_GEO_DATA), TARGET,     &
                               intent(IN) :: GEO_DATA
    integer,         intent(IN) :: item
    type (CELL_NEIGHBOR_DATA),POINTER     :: cell_ngbr

    ! Local variables
    type (CELL_NEIGHBOR_DATA),dimension(:),&
                              POINTER     :: cell_ngbrs => NULL()

    cell_ngbrs => NGBR_GEO_DATA(geo_data)
    cell_ngbr  => cell_ngbrs(item)
  end function geo_data_item_cell_ngbr_data
  
  function geo_data_item_fnumber_face(geo_data, item, face_number) RESULT(face)
    type (LONG_EDIT_GEO_DATA), TARGET,     &
                             intent(IN) :: GEO_DATA
    integer,       intent(IN) :: item
    integer,       intent(IN) :: face_number
    type (FACE_DATA),      POINTER    :: Face

    ! Local variables
    type (CELL_NEIGHBOR_DATA), POINTER :: CELL_NGBR
    
    Cell_Ngbr => NGBR_GEO_DATA(geo_data, item)
    Face      => FACE_GEO_DATA(cell_Ngbr, face_number)
  end function geo_data_item_fnumber_face
  

  !!!!!!!!!! Accessor routines for LONG_EDIT_LIST type and MECH_LIST_TYPE !!!!!!!!!

  function list_cellnumber(LIST, ITEM) RESULT(Cell)
    type (LONG_EDIT_LIST), intent(IN) :: LIST
    integer,     intent(IN) :: ITEM
    integer                 :: Cell

    Cell = LIST%CellNumber(ITEM)
  end function list_cellnumber

  function mech_list_cellnumber(LIST, ITEM) RESULT(Cell)
    type (MECH_EDIT_LIST), intent(IN) :: LIST
    integer,     intent(IN) :: ITEM
    integer                 :: Cell

    Cell = LIST%CellNumber(ITEM)
  end function mech_list_cellnumber

  function list_phys_data(List) RESULT(Data)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: LIST
    type (LONG_EDIT_PHYS_DATA), POINTER    :: Data

    DATA => list%Phys_data
  end function list_phys_data
  
  function list_geo_data(List) RESULT(Data)
    type (LONG_EDIT_LIST),     TARGET,     &
                               intent(IN) :: LIST
    type (LONG_EDIT_GEO_DATA), POINTER    :: Data

    DATA => list%geo_data
  end function list_geo_data
  
  function list_mech_data(List) RESULT(Data)
    type (MECH_EDIT_LIST),     TARGET,     &
                               intent(IN) :: LIST
    type (CELL_MECH_INVARIANT), POINTER, dimension(:)    :: Data

    DATA => list%mech_data
  end function list_mech_data
  
  function list_item_zone(list, item) RESULT(zone)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    integer,     intent(IN) :: item
    type (CELL_AVG), POINTER          :: Zone

    ! Local variables
    type (LONG_EDIT_PHYS_DATA), POINTER    :: LData

    LData => PHYS_DATA(List)
    Zone => ZONE_DATA(LData, ITEM)
  end function list_item_zone
  
  function list_item_mech(list, item) RESULT(mech)
    type (MECH_EDIT_LIST), intent(IN) :: List
    integer,     intent(IN) :: item
    type (CELL_MECH_INVARIANT), POINTER     :: mech

    ! Local variables
    type (CELL_MECH_INVARIANT), POINTER, dimension(:)    :: MData => NULL()

    MData => List_Mech_Data(List)
    mech => MData(ITEM)
  end function list_item_mech
  
  function list_item_matl(list, item, matl_slot) RESULT(matl)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    integer,     intent(IN) :: item
    integer,     intent(IN) :: matl_slot
    type (MATERIAL), POINTER          :: Matl

    ! Local variables
    type (LONG_EDIT_PHYS_DATA), POINTER    :: LData => NULL()

    LData => PHYS_DATA(List)

    Matl => LDATA%matl(matl_slot)%Cell(item)
  end function list_item_matl

  function list_item_cell_geo(list, item) RESULT(geo)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    integer,     intent(IN) :: item
    type (CELL_GEOMETRY), POINTER     :: Geo

    ! Local variables
    type (LONG_EDIT_GEO_DATA), POINTER    :: LData => NULL()

    LData => GEO_DATA(List)
    Geo => CELL_GEO_DATA(LData, ITEM)
  end function list_item_cell_geo
  
  function list_item_cell_vertex_data(list, item) RESULT(cell_vertex)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    integer,     intent(IN) :: item
    type (CELL_VERTEX_DATA),POINTER   :: cell_vertex

    ! Local variables
    type (LONG_EDIT_GEO_DATA), POINTER    :: LData => NULL()

    LData       => GEO_DATA(List)
    cell_vertex => VERTEX_GEO_DATA(LData, ITEM)
  end function list_item_cell_vertex_data
  
  function list_item_cell_vnumber_vertex(list, item, vertex_number) RESULT(vrtx)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    integer,     intent(IN) :: item
    integer,     intent(IN) :: vertex_number
    type (VERTEX_DATA),    POINTER    :: vrtx

    ! Local variables
    type (CELL_VERTEX_DATA), POINTER    :: CELL_VERTEX => NULL()

    CELL_VERTEX => VERTEX_GEO_DATA(List, item)
    vrtx        => VERTEX_GEO_DATA(CELL_VERTEX, VERTEX_NUMBER)
  end function list_item_cell_vnumber_vertex
  
  function list_item_cell_ngbr_data(list, item) RESULT(cell_ngbr)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    integer,     intent(IN) :: item
    type (CELL_NEIGHBOR_DATA),POINTER :: cell_ngbr

    ! Local variables
    type (LONG_EDIT_GEO_DATA), POINTER    :: LData => NULL()

    LData       => GEO_DATA(List)
    cell_ngbr => NGBR_GEO_DATA(LData, ITEM)
  end function list_item_cell_ngbr_data
  
  function list_item_cell_fnumber_face(list, item, face_number) RESULT(face)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    integer,     intent(IN) :: item
    integer,     intent(IN) :: face_number
    type (FACE_DATA),      POINTER    :: face

    ! Local variables
    type (CELL_NEIGHBOR_DATA), POINTER    :: CELL_NGBR => NULL()

    CELL_NGBR => NGBR_GEO_DATA(List, item)
    face      => FACE_GEO_DATA(CELL_NGBR, FACE_NUMBER)
  end function list_item_cell_fnumber_face
  
  function list_zone(list) RESULT(zone)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    type (CELL_AVG),       dimension(:),&
                           POINTER    :: Zone

    ! Local variables
    type (LONG_EDIT_PHYS_DATA), POINTER    :: LData => NULL()

    LData => PHYS_DATA(List)
    Zone  => ZONE_DATA(LData)
  end function list_zone
  
  function list_matl_slot_matl(list, matl_slot) RESULT(matl)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    integer,     intent(IN) :: matl_slot
    type (MATERIAL),       dimension(:),&
                           POINTER    :: Matl

    ! Local variables
    type (LONG_EDIT_PHYS_DATA), POINTER    :: LData => NULL()

    LData => PHYS_DATA(List)

    Matl  => LDATA%matl(matl_slot)%Cell
  end function list_matl_slot_matl

  function list_matl(list) RESULT(matl)
    type (LONG_EDIT_LIST), TARGET,     &
                           intent(IN) :: List
    type (MATL_SLOT),      dimension(:),&
                           POINTER    :: Matl

    ! Local variables
    type (LONG_EDIT_PHYS_DATA), POINTER    :: LData => NULL()

    LData => PHYS_DATA(List)

    Matl => LDATA%matl
  end function list_matl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!           SORT and COLLATE Routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine phys_data_collate(COLLATED_DATA, LOCAL_DATA)
    type (LONG_EDIT_PHYS_DATA), target, intent(INOUT) :: COLLATED_DATA
    type (LONG_EDIT_PHYS_DATA), target, intent(IN) :: LOCAL_DATA

    ! Local variables
    type (CELL_AVG), dimension(:), POINTER :: Local_Zone  => NULL()
    type (CELL_AVG), dimension(:), POINTER :: Collated_Zone  => NULL()
    type (MATL_SLOT),dimension(:), POINTER :: Collated_Matl  => NULL()
    type (MATL_SLOT),dimension(:), POINTER :: Local_Matl  => NULL()
    

    ! Go through each item in the LONG_EDIT_PHYS_DATA type individually
    ! Be sure to edit this code when you add new data to the LONG_EDIT_PHYS_DATA
       
    ! Take care of Zone
    Local_Zone    => ZONE_DATA(LOCAL_DATA)
    Collated_Zone => ZONE_DATA(COLLATED_DATA)

    call COLLATE(Collated_Zone, Local_Zone)

    ! Now take care of Material
    Local_Matl    => MATL_DATA(LOCAL_DATA)
    Collated_Matl => MATL_DATA(COLLATED_DATA)
    
    call COLLATE(Collated_Matl, Local_Matl)

  end subroutine phys_data_collate
  
  subroutine phys_data_sort(DATA, RANK)
    use parallel_scope
    type (LONG_EDIT_PHYS_DATA), intent(INOUT) :: Data
    integer,     dimension(:),  &
                           intent(IN) :: Rank

    ! Local variables
    integer :: s
    type (LONG_EDIT_PHYS_DATA) :: Temp_Data
    type (CELL_AVG), dimension(:), POINTER :: Temp_Zone => NULL()
    type (CELL_AVG), dimension(:), POINTER :: Permuted_Zone => NULL()
    type (MATERIAL), dimension(:), POINTER :: Temp_Matl => NULL()
    type (MATERIAL), dimension(:), POINTER :: Permuted_Matl => NULL()

    ! First copy the original data into the temporary
    call CLONE(NEW=Temp_DATA, OLD = DATA)
    Temp_Data = Data
    
    ! Now permute the temporary data into the original data
    ! Take care of Zone first
    Temp_Zone     => ZONE_DATA(Temp_DATA)
    Permuted_Zone => ZONE_DATA(DATA)
    call PERMUTE_ZONE(Permuted_Zone, Temp_Zone, RANK, SCOPE=LOCAL_SCOPE)

    ! Now permute the matl stuff 
    do s = 1, mat_slot
       Temp_Matl     => MATL_DATA(Temp_DATA, MATL_SLOT=s)
       Permuted_Matl => MATL_DATA(DATA,      MATL_SLOT=s)
       call PERMUTE_MATL(Permuted_Matl, Temp_Matl, RANK, SCOPE=LOCAL_SCOPE)
    end do

    ! Now get rid of temporary and go home
    call DESTROY(Temp_Data)

  end subroutine phys_data_sort

  subroutine geo_data_collate(COLLATED_DATA, LOCAL_DATA)
    use mesh_module,    only : CELL_GEOMETRY, COLLATE
    type (LONG_EDIT_GEO_DATA), intent(INOUT) :: COLLATED_DATA
    type (LONG_EDIT_GEO_DATA), intent(IN) :: LOCAL_DATA

    ! Local variables
    type (CELL_GEOMETRY), dimension(:), POINTER :: Local_cell_geo => NULL()
    type (CELL_GEOMETRY), dimension(:), POINTER :: Collated_cell_geo => NULL()
    type (CELL_VERTEX_DATA), dimension(:), POINTER :: Local_vrtx_geo => NULL()
    type (CELL_VERTEX_DATA), dimension(:), POINTER :: Collated_vrtx_geo => NULL()
    type (CELL_NEIGHBOR_DATA),dimension(:),POINTER :: Local_ngbr_geo => NULL()
    type (CELL_NEIGHBOR_DATA),dimension(:),POINTER :: Collated_ngbr_geo => NULL()


    ! Go through each item in the LONG_EDIT_GEO_DATA type individually
    ! Be sure to edit this code when you add new data to the LONG_EDIT_GEO_DATA
       
    ! Take care of Cell Geometry information
    Local_cell_geo    => CELL_GEO_DATA(LOCAL_DATA)
    Collated_cell_geo => CELL_GEO_DATA(COLLATED_DATA)

    call COLLATE(Collated_cell_geo, Local_cell_geo)

    ! Take care of Vertex Geometry information
    Local_vrtx_geo    => VERTEX_GEO_DATA(LOCAL_DATA)
    Collated_vrtx_geo => VERTEX_GEO_DATA(COLLATED_DATA)

    call COLLATE(Collated_vrtx_geo, Local_vrtx_geo)
    
    ! Take care of Neigbhor Geometry information
    Local_ngbr_geo    => NGBR_GEO_DATA(LOCAL_DATA)
    Collated_ngbr_geo => NGBR_GEO_DATA(COLLATED_DATA)

    call COLLATE(Collated_ngbr_geo, Local_ngbr_geo)

  end subroutine geo_data_collate
  
  subroutine geo_data_sort(DATA, RANK)
    use parallel_scope
    type (LONG_EDIT_GEO_DATA), intent(INOUT) :: Data
    integer,         dimension(:),  &
                               intent(IN) :: Rank

    ! Local variables
    type (LONG_EDIT_GEO_DATA) :: Temp_Data
    type (CELL_GEOMETRY), dimension(:), POINTER :: Temp_cell_Geo => NULL()
    type (CELL_GEOMETRY), dimension(:), POINTER :: Permuted_cell_Geo => NULL()
    type (CELL_VERTEX_DATA), dimension(:), POINTER :: Temp_vertex_Geo => NULL()
    type (CELL_VERTEX_DATA), dimension(:), POINTER :: Permuted_vertex_Geo => NULL()
    type (CELL_NEIGHBOR_DATA), dimension(:), POINTER :: Temp_ngbr_Geo => NULL()
    type (CELL_NEIGHBOR_DATA), dimension(:), POINTER :: Permuted_ngbr_Geo => NULL()


    ! First copy the original data into the temporary
    call CLONE(NEW=Temp_DATA, OLD = DATA)
    Temp_Data = Data
    
    ! Now permute the temporary data into the original data
    ! Take care of Cell geometry first
    Temp_cell_Geo     => CELL_GEO_DATA(Temp_DATA)
    Permuted_cell_Geo => CELL_GEO_DATA(DATA)
    call PERMUTE_Cell(Permuted_cell_Geo, Temp_cell_Geo, RANK, SCOPE=LOCAL_SCOPE)

    ! THen take care of vertex info
    Temp_vertex_Geo     => VERTEX_GEO_DATA(Temp_DATA)
    Permuted_vertex_Geo => VERTEX_GEO_DATA(DATA)
    call PERMUTE_vertex_list(Permuted_vertex_Geo, Temp_vertex_Geo, RANK, SCOPE=LOCAL_SCOPE)

    ! Don't forget the neighbor info
    Temp_ngbr_Geo     => NGBR_GEO_DATA(Temp_DATA)
    Permuted_ngbr_Geo => NGBR_GEO_DATA(DATA)
    call PERMUTE_ngbr_list(Permuted_ngbr_Geo, Temp_ngbr_Geo, RANK, SCOPE=LOCAL_SCOPE)

    ! Now get rid of temporary and go home
    call DESTROY(Temp_Data)

  end subroutine geo_data_sort

  subroutine list_collate(COLLATED_LIST, LOCAL_LIST)
    use pgslib_module, only : pgslib_collate
    type (LONG_EDIT_LIST), intent(INOUT) :: LOCAL_LIST
    type (LONG_EDIT_LIST), intent(INOUT) :: COLLATED_LIST

    ! Go through each item in the LONG_EDIT_LIST type individually
    ! Be sure to edit this code when you add new data to the LONG_EDIT_LIST
    call pgslib_collate(COLLATED_LIST%CellNumber, LOCAL_LIST%CellNumber)
    
    call COLLATE(COLLATED_LIST%Phys_Data, LOCAL_LIST%Phys_Data)    

    call COLLATE(COLLATED_LIST%Geo_Data,  LOCAL_LIST%Geo_Data)    
    
  end subroutine list_collate

  subroutine mech_list_collate(COLLATED_LIST, LOCAL_LIST)
    use pgslib_module, only : pgslib_collate
    type (MECH_EDIT_LIST), intent(INOUT) :: LOCAL_LIST
    type (MECH_EDIT_LIST), intent(INOUT) :: COLLATED_LIST

    call pgslib_collate(COLLATED_LIST%CellNumber, LOCAL_LIST%CellNumber)
    
    call PGSLib_COLLATE(COLLATED_LIST%Mech_Data%mises_stress, LOCAL_LIST%Mech_Data%mises_stress) 
    call PGSLib_COLLATE(COLLATED_LIST%Mech_Data%eff_plastic_strain, LOCAL_LIST%Mech_Data%eff_plastic_strain) 
    call PGSLib_COLLATE(COLLATED_LIST%Mech_Data%mean_stress, LOCAL_LIST%Mech_Data%mean_stress) 
    call PGSLib_COLLATE(COLLATED_LIST%Mech_Data%volumetric_strain, LOCAL_LIST%Mech_Data%volumetric_strain) 
    
  end subroutine mech_list_collate
    
  subroutine list_sort(LIST)
    ! This is a local sort, stuff gets sorted locally on each processor only
    use pgslib_module
    type (LONG_EDIT_LIST), target, intent(INOUT) :: LIST

    ! Local variables
    integer :: i
    integer, allocatable, dimension(:) :: ItemRank
    integer, allocatable, dimension(:) :: Temp_Cell_Number
    type(LONG_EDIT_PHYS_DATA),  POINTER          :: P_Data => NULL()
    type(LONG_EDIT_GEO_DATA),   POINTER          :: L_Geo_Data => NULL()
    integer                                      :: status

    allocate (ItemRank(SIZE(LIST)), STAT=status)
    if (status /= 0) call TLS_panic ('SORT_LONG_EDIT_DATA: allocate failed: ItemRank')

    allocate (Temp_Cell_Number(SIZE(LIST)), STAT=status)
    if (status /= 0) call TLS_panic ('SORT_LONG_EDIT_DATA: allocate failed: Temp_Cell_Number')

    ! Find the rank of the global cell numbers.
    ItemRank = PGSLib_GRADE_UP_LOCAL(List%CellNumber)

    ! Now sort the items.  
    Temp_Cell_Number = List%CellNumber
    do i = 1, SIZE(List%CellNumber)
       List%CellNumber(ItemRank(i)) = Temp_Cell_Number(i)
    end do

    P_Data => PHYS_DATA(List)
    call SORT(DATA=P_Data, RANK = ItemRank)

    L_GEO_Data => GEO_DATA(List)
    call SORT(DATA=L_GEO_Data, RANK = ItemRank)

    deallocate (Temp_Cell_Number)
    deallocate (ItemRank)

  end subroutine list_sort
    
  subroutine mech_list_sort(LIST)
    ! This is a local sort, stuff gets sorted locally on each processor only
    use pgslib_module
    type (MECH_EDIT_LIST), intent(INOUT) :: LIST

    ! Local variables
    integer :: i
    integer, allocatable, dimension(:) :: ItemRank
    integer, allocatable, dimension(:) :: Temp_Cell_Number
    type(CELL_MECH_INVARIANT),allocatable, dimension(:) :: M_Data
    integer                                      :: status

    allocate (ItemRank(SIZE(LIST)), STAT=status)
    if (status /= 0) call TLS_panic ('MECH_LIST_SORT: allocate failed: ItemRank')

    allocate (Temp_Cell_Number(SIZE(LIST)), STAT=status)
    if (status /= 0) call TLS_panic ('MECH_LIST_SORT: allocate failed: Temp_Cell_Number')

    allocate (M_Data(SIZE(LIST)), STAT=status)
    if (status /= 0) call TLS_panic ('MECH_LIST_SORT: allocate failed: M_Data')

    ! Find the rank of the global cell numbers.
    ItemRank = PGSLib_GRADE_UP_LOCAL(List%CellNumber)

    ! Now sort the items.  
    Temp_Cell_Number = List%CellNumber
    M_Data = List%Mech_Data
    do i = 1, SIZE(List%CellNumber)
       List%CellNumber(ItemRank(i)) = Temp_Cell_Number(i)
       List%Mech_Data(ItemRank(i))%mises_stress = M_Data(i)%mises_stress
       List%Mech_Data(ItemRank(i))%eff_plastic_strain = M_Data(i)%eff_plastic_strain
       List%Mech_Data(ItemRank(i))%mean_stress = M_Data(i)%mean_stress
       List%Mech_Data(ItemRank(i))%volumetric_strain = M_Data(i)%volumetric_strain
    end do

    deallocate (Temp_Cell_Number)
    deallocate (ItemRank)
    deallocate (M_Data)

  end subroutine mech_list_sort
    
  
END MODULE LONG_EDIT_DATA_TYPES
