!!#define DEBUG_MESH_MODULE

MODULE MESH_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define structures and procedures associated with mesh connectivity.
  !
  !   Public Interface(s):
  !
  !     * call ASSIGN_VRTX_BITS (Vrtx_Bit)
  !
  !       Assign bit positions in integer array Vrtx_Bit.
  !
  !     * call ASSIGN_CELL_BITS (Cell_Bit)
  !
  !         Assign bit positions in integer array Cell_Bit.
  !
  !     * call ASSIGN_CELL_EDGES (Cell_Edge)
  !
  !         Assign edge positions in integer array Cell_Edge.
  !
  ! Contains: ASSIGN_CELL_BITS
  !           ASSIGN_CELL_EDGES
  !           ASSIGN_VRTX_BITS
  !           VERTEX_PRESET_SCALAR
  !
  ! Author(s): Douglas B. Kothe, LANL (dbk@lanl.gov)
  !            Modifed by ferrell@cpca.com
  !
  !=======================================================================
  use constants_module, only: preset
  use kind_module,      only: int_kind, log_kind, real_kind
  use parameter_module, only: ndim, nfc, nvc, nvf, nec
  use var_vector_types

  implicit none

  ! Private Module
  private

  ! Public Variables and Types
  public :: CELL_BIT, MESH_CONNECTIVITY, VERTEX_BIT, CELL_EDGE, &
            CELL_GEOMETRY, VERTEX_DATA, BOUNDARY, Vrtx_Bdy, &
            Cell, Mesh, Vertex, Mesh_Connectivity_Preset, & 
            CELL_TET, CELL_PYRAMID, CELL_PRISM, CELL_HEX, &
            GAP_ELEMENT_1,GAP_ELEMENT_3,GAP_ELEMENT_5

  public :: Vertex_Ngbr_All, Vertex_Ngbr_All_Orig
  public :: MESH_COLLATE_VERTEX, VERTEX_COLLATE

  ! Public Subroutines
  public :: ASSIGN_CELL_BITS, ASSIGN_VRTX_BITS, ASSIGN_CELL_EDGES, &
            Vertex_Data_Preset
  public :: Initialize_Face_Bit_Mask,        &
            Set_Face_Neighbor,               &
            Clear_Face_Neighbor,             &
            Is_Face_Ngbr

  public :: COLLATE, PERMUTE_CELL
  public :: read_mesh_data

  ! Public Operators
  public :: operator(.DISTRIBUTE.)

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Various mesh flags.
  logical(KIND = log_kind), public, save :: orthogonal_mesh

  !! Parameters to identify the state of cells, faces and nodes.
  ! Degenerate face and node parameters.
  integer(int_kind), parameter, public :: EXTERNAL_FACE   = 0
  integer(int_kind), parameter, public :: DEGENERATE_FACE = -HUGE(int_kind)
  integer(int_kind), parameter, public :: DEGENERATE_NODE = -1
  ! Define unknown, default, vertex (node), cell, face and PE numbers
  integer(int_kind), parameter :: DEFAULT_CELL   = -1
  integer(int_kind), parameter :: DEFAULT_VERTEX = -1
  integer(int_kind), parameter :: DEFAULT_FACE   = -1
  integer(int_kind), parameter :: DEFAULT_PE     = -1

  ! Mesh diagnostics variables.
  integer(KIND = int_kind), public, save :: degenerate_points, degenerate_lines, &
                                            triangle_faces, quad_faces
  real(KIND = real_kind),   public, save :: volume_min, volume_max

  ! Face-Vertex and Vertex-Face mappings.
  integer(KIND = int_kind), public, save, dimension(nfc,nvf)  :: Face_Vrtx
  integer(KIND = int_kind), public, save, dimension(nvc,ndim) :: Vrtx_Face
  
  ! MESH_CONNECTIVITY structure
  Type MESH_CONNECTIVITY
     Integer (KIND=int_kind), Dimension(nfc)  :: Ngbr_Cell      ! cell number across each face
     Integer (KIND=int_kind), Dimension(nfc)  :: Ngbr_Cell_PE   ! holding data for Ngbr_Cell
     Integer (KIND=int_kind), Dimension(nfc)  :: Ngbr_Cell_Orig ! original cell numbers across each face
     Integer (KIND=int_kind), Dimension(nfc)  :: Ngbr_Face      ! face number of cell across each face
     Integer (KIND=int_kind), Dimension(nvc)  :: Ngbr_Vrtx      ! cell vertex numbers
     Integer (KIND=int_kind), Dimension(nvc)  :: Ngbr_Vrtx_PE   ! holding data for Ngbr_Vrtx
     Integer (KIND=int_kind), Dimension(nvc)  :: Ngbr_Vrtx_Orig ! original vertex numbers
     Type (int_var_vector)                    :: Ngbr_Cells_ALL ! all neighbor cells
     Type (int_var_vector)                    :: Ngbr_Cells_ALL_Orig ! Orig, global, cell numbers, 
     Type (int_var_vector)                    :: Ngbr_Cells_Face ! Bits to identify face ngbrs
     Integer                                  :: Cell_shape     ! indicator for hex, tet, ...
     Integer                                  :: CBlockID = 0   ! cell block ID
  End Type MESH_CONNECTIVITY

  ! define instance of MESH_CONNECTIVITY type
  type(MESH_CONNECTIVITY), dimension(:), pointer, SAVE :: Mesh => NULL()
  logical, save, public :: mesh_has_cblockid_data = .false.

  ! define the four cell shapes, to be stored in Mesh%Cell_shape, 
  !   determined in MESH_TYPE
  integer, parameter :: CELL_TET     = 1
  integer, parameter :: CELL_PYRAMID = 2
  integer, parameter :: CELL_PRISM   = 3
  integer, parameter :: CELL_HEX     = 4
  ! define another cell shape for "zero volume" gap elements
  ! The gap element number corresponds to the first non-degenerate cell face
  integer, parameter :: GAP_ELEMENT_1  = 5
  integer, parameter :: GAP_ELEMENT_3  = 6
  integer, parameter :: GAP_ELEMENT_5  = 7

  ! Permutation vectors for mesh and vertices (used by Mesh_Permute)
  integer(int_kind), POINTER,      &
                     PUBLIC,       &
                     SAVE,         &
                     dimension(:) :: UnPermute_Mesh_Vector => NULL()
  integer(int_kind), POINTER,      &
                     PUBLIC,       &
                     SAVE,         &
                     dimension(:) :: Permute_Mesh_Vector => NULL()
  logical(log_kind), SAVE,         &
                     PUBLIC       :: UnPermute_Mesh_Initialized = .FALSE.
  integer(int_kind), POINTER,      &
                     PUBLIC,       &
                     SAVE,         &
                     dimension(:) :: UnPermute_Vertex_Vector => NULL()
  integer(int_kind), POINTER,      &
                     PUBLIC,       &
                     SAVE,         &
                     dimension(:) :: Permute_Vertex_Vector => NULL()
  logical(log_kind), SAVE,         &
                     PUBLIC       :: UnPermute_Vertex_Initialized = .FALSE.

  ! Define VERTEX_BIT structure
  type VERTEX_BIT
     integer(KIND = int_kind), dimension(nvc) :: Bit
  end type VERTEX_BIT

  ! Declare a VERTEX_BIT type
  type(VERTEX_BIT), public, save :: Vrtx

  ! Define CELL_BIT structure
  type CELL_BIT
     integer(KIND = int_kind), dimension(nfc) :: Bit
  end type CELL_BIT

  ! Declare a CELL_BIT instance
  type(CELL_BIT), public, save :: CllNgbr

  ! Define the edges surrounding a cell
  integer(KIND = int_kind), dimension(2,nec) :: Cell_Edge

  ! define CELL_GEOMETRY derived type
  Type CELL_GEOMETRY
     real (real_kind), dimension(ndim)     :: Centroid        ! cell centroid (physical coordinates)
!    real (real_kind), dimension(ndim)     :: Centroid_L      ! cell centroid (logical coordinates) (place holder)
     real (real_kind)                      :: Volume          ! cell volume
     real (real_kind), dimension(nfc)      :: Halfwidth       ! cell halfwidth
     real (real_kind), dimension(ndim,nfc) :: Face_Centroid   ! face centroid (physical coordinates)
     real (real_kind), dimension(ndim,nfc) :: Face_Centroid_L ! face centroid (logical coordinates)
     real (real_kind), dimension(nfc)      :: Face_Area       ! face area
     real (real_kind), dimension(ndim,nfc) :: Face_Normal     ! face unit normal
  End Type CELL_GEOMETRY

  ! Declare instance of CELL_GEOMETRY type.
  type(CELL_GEOMETRY), dimension(:), SAVE, pointer :: Cell => NULL()

  ! Vertex connectivity
  Type(int_var_vector), dimension(:), POINTER :: Vertex_Ngbr_All => NULL()
  Type(int_var_vector), dimension(:), POINTER :: Vertex_Ngbr_All_Orig => NULL()

  ! VERTEX_DATA structure
  Type VERTEX_DATA
     Real (KIND=real_kind) :: Coord(ndim)   ! vertex coordinates
     Real (KIND=real_kind) :: Rsum_rvol     ! reciprocal sum of inverse volumes around vertices
  End Type VERTEX_DATA

  ! Ghost cells for the vertex data
  ! **DANGER** If coord data changes these need to be nullified.
  ! Define type for boundary data
  type BOUNDARY
    real(KIND = real_kind), dimension(:), pointer :: Data => NULL()
  end type BOUNDARY

  ! This is needed because F90 may not put pointers into known state
  ! When we have gone fully to F95 this will no longer be necessary.
  ! Only place this is used is in Vertex_Preset_Scalar%
  logical(KIND = log_kind), dimension(ndim), SAVE :: Bdy_Initialized = .false.
  
  type(BOUNDARY), dimension(ndim), SAVE :: Vrtx_Bdy

  ! Default Mesh Connectivity Data
  interface Mesh_Connectivity_Preset
     module procedure Mesh_Preset_Scalar
  end interface

  ! Default Vertex Data
  interface VERTEX_DATA_PRESET
     module procedure VERTEX_PRESET_SCALAR
  end interface

  ! Define instance of VERTEX_DATA type.
  type(VERTEX_DATA), dimension(:), pointer, save :: Vertex => NULL()

  ! Bit mask for identifying face neighbors
  integer(int_kind), dimension(:), pointer :: Face_Bit_Mask => NULL()
  logical(log_kind), SAVE                  :: Face_Bit_Mask_Initialized = .FALSE.

  interface COLLATE
     module procedure COLLATE_CELL
  end interface

  interface PERMUTE_CELL
     module procedure PERMUTE_CELL
  end interface

  interface operator(.DISTRIBUTE.)
     module procedure MESH_DISTRIBUTE
     module procedure VERTEX_DISTRIBUTE
     module procedure CELL_DISTRIBUTE
  end interface
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  FUNCTION MESH_PRESET_SCALAR () RESULT(MESH_PRESET)
    !====================================================================
    ! Purpose(s):
    !
    !   Default scalar for mesh connectivity
    !   This is a funciton so we can nullify the Ngbr_Cells_All field
    !====================================================================

    ! Return value
    type(MESH_CONNECTIVITY)        :: Mesh_Preset

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize.
    Mesh_Preset%Ngbr_Cell      = DEFAULT_CELL
    Mesh_Preset%Ngbr_Cell_PE   = DEFAULT_PE

    Mesh_Preset%Ngbr_Face      = DEFAULT_FACE

    Mesh_Preset%Ngbr_Vrtx      = DEFAULT_VERTEX
    Mesh_Preset%Ngbr_Vrtx_PE   = DEFAULT_PE

    Mesh_Preset%Ngbr_Cell_Orig = DEFAULT_CELL
    Mesh_Preset%Ngbr_Vrtx_Orig = DEFAULT_VERTEX
    NULLIFY(Mesh_Preset%Ngbr_Cells_All%v)

    return

  END FUNCTION MESH_PRESET_SCALAR

  FUNCTION VERTEX_PRESET_SCALAR ()
    !====================================================================
    ! Purpose(s):
    !
    !   Assign Vertex_Data a default value
    !   This is a function so we can nullify the boundary stuff
    !
    !====================================================================
    use kind_module, only: int_kind

    implicit none

    ! Local Variables
    integer(KIND = int_kind) :: n
    type(VERTEX_DATA)        :: Vertex_Preset_Scalar

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize.
    Vertex_Preset_Scalar%Coord     = preset
    Vertex_Preset_Scalar%Rsum_rvol = preset

    do n = 1, ndim

       if (.not.Bdy_Initialized(n)) then
          NULLIFY(Vrtx_Bdy(n)%Data)
          Bdy_Initialized(n) = .true.
       else
          if (ASSOCIATED(Vrtx_Bdy(n)%Data)) then
             DEALLOCATE(Vrtx_Bdy(n)%Data)
          else
             NULLIFY(Vrtx_Bdy(n)%Data)
          end if
       end if

    end do

  END FUNCTION VERTEX_PRESET_SCALAR

  SUBROUTINE ASSIGN_CELL_BITS (Cell_Bit)
    !=======================================================================
    ! Purpose(s):
    !
    !   Assign bit positions for the integer Ngbr_CELL_pe_flag which
    !   is part of the Mesh structure. Bits 8-13 of 
    !   Ngbr_pe_flag correspond to the 6 cells surrounding each cell
    !   A bit set to 1 means the cell is on-processor.
    !
    !=======================================================================
    use kind_module,      only: int_kind
    use parameter_module, only: nfc

    implicit none

    ! Arguments
    integer(KIND = int_kind), dimension(nfc), intent(OUT) :: Cell_Bit

    ! Local Variables
    integer(KIND = int_kind) :: c

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! The following convention holds for the bits in Mesh%Ngbr_pe_flag:
    !
    !      Bit     Cell
    !      ---     ------
    !      8-14      1-6   (1 bit per cell)

    ! Assign bit locations for the CELL type
    do c = 1,nfc
       Cell_Bit(c) = c - 1 + SIZE(Vrtx%Bit,1)
    end do

    return

  END SUBROUTINE ASSIGN_CELL_BITS

  SUBROUTINE ASSIGN_VRTX_BITS (Vrtx_bit)
    !=======================================================================
    ! Purpose(s):
    !
    !   Assign bit positions for the integer Ngbr_vrtx_pe_flag which
    !   is part of the Mesh structure. The first 8 bits in
    !   Ngbr_pe_flag correspond to the 8 vertices of each cell;
    !   A bit set to 1 means the vertex is on-processor.
    !
    !=======================================================================
    use kind_module,      only: int_kind
    use parameter_module, only: nvc

    implicit none

    ! Argument List
    integer(KIND = int_kind), dimension(nvc), intent(OUT) :: Vrtx_bit

    ! Local Variables
    integer(KIND = int_kind) :: v

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! The following convention holds for the bits in Mesh%Vrtx_pe_flag:
    !
    !      Bit     Vertex
    !      ---     ------
    !      0-7      1-8   (1 bit per vertex)

    ! Assign bit locations for the Vrtx type
    do v = 1,nvc
       Vrtx_bit(v) = v - 1
    end do

    return

  END SUBROUTINE ASSIGN_VRTX_BITS

  SUBROUTINE ASSIGN_CELL_EDGES (Cell_Edge)
    !=======================================================================
    ! Purpose(s):
    !
    !    Assign cell edges.
    !
    !=======================================================================
    use kind_module,      only: int_kind
    use parameter_module, only: nec

    implicit none

    ! Argument List
    integer(KIND = int_kind), dimension(2,nec), intent(OUT) :: Cell_Edge

    ! Local Variables
    integer(KIND = int_kind) :: edge

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do edge = 1,nec
       select case (edge)
          case (1)
             Cell_Edge(:,edge) = (/1,2/)
          case (2)
             Cell_Edge(:,edge) = (/2,3/)
          case (3)
             Cell_Edge(:,edge) = (/3,4/)
          case (4)
             Cell_Edge(:,edge) = (/4,1/)
          case (5)
             Cell_Edge(:,edge) = (/2,6/)
          case (6)
             Cell_Edge(:,edge) = (/3,7/)
          case (7)
             Cell_Edge(:,edge) = (/4,8/)
          case (8)
             Cell_Edge(:,edge) = (/1,5/)
          case (9)
             Cell_Edge(:,edge) = (/5,6/)
          case (10)
             Cell_Edge(:,edge) = (/6,7/)
          case (11)
             Cell_Edge(:,edge) = (/7,8/)
          case (12)
             Cell_Edge(:,edge) = (/8,5/)
       end select
    end do

    return

  END SUBROUTINE ASSIGN_CELL_EDGES

  SUBROUTINE Initialize_Face_Bit_Mask()
    !==================================================================
    ! Purpose(s):
    !   Initialize the array of bit masks which are used to flag
    !   face neighbors.  Also, set the initialization flag to .TRUE.
    !==================================================================
    
    ! Local variables
    integer :: face


    if (.NOT. Face_Bit_Mask_Initialized) then
       ! Allocate the array, and clear all bits
       ALLOCATE(Face_Bit_Mask(nfc))
       Face_Bit_Mask = 0
       
       do face = 1, nfc
          Face_Bit_Mask(face) = IBSET(Face_Bit_Mask(face), face)
       end do
    
       Face_Bit_Mask_Initialized = .TRUE.
    end if

    RETURN
  END SUBROUTINE Initialize_Face_Bit_Mask
       

  SUBROUTINE Set_Face_Neighbor(Face_Bits, Face)
    !==================================================================
    ! Purpose(s):
    !   Set the neighbor flag to indicate that it is a neighbor
    !   of Face.
    !==================================================================
    
    integer(int_kind), intent(INOUT) :: Face_Bits
    integer(int_kind), intent(IN   ) :: Face

    Face_Bits = IOR(Face_Bits, Face_Bit_Mask(Face))
      
    RETURN
  END SUBROUTINE Set_Face_Neighbor
    
  SUBROUTINE Clear_Face_Neighbor(Face_Bits, Face)
    !==================================================================
    ! Purpose(s):
    !   Clear the neighbor flag which indicates that it is a neighbor
    !   of Face.
    !==================================================================
    
    integer(int_kind), intent(INOUT) :: Face_Bits
    integer(int_kind), intent(IN   ) :: Face

    Face_Bits = IAND(Face_Bits, NOT(Face_Bit_Mask(Face)))
      
    RETURN
  END SUBROUTINE Clear_Face_Neighbor
      
  FUNCTION Is_Face_Ngbr(Face_Bits, Face)
    !==================================================================
    ! Purpose(s):
    !   Test the flag of an integer to see if that neighbor is also
    !   a neighbor of Face.
    
    !   Set the neighbor flag to indicate that it is a neighbor
    !   of Face.
    !==================================================================
    use kind_module, only: int_kind
    implicit none
        
    integer(int_kind), intent(IN   ) :: Face_Bits
    integer(int_kind), intent(IN   ) :: Face
    logical(log_kind)                :: Is_Face_Ngbr

    Is_Face_Ngbr = (IAND(Face_Bits, Face_Bit_Mask(Face)) /= 0)
      
    RETURN
  END FUNCTION Is_Face_Ngbr
      
  FUNCTION MESH_COLLATE_VERTEX (Mesh)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed mesh into a single large mesh on IO PE
    !==================================================================
    use parallel_info_module, only: p_info
    use parameter_module,     only: ncells_tot, nvc, ncells
    use pgslib_module,        only: PGSLib_COLLATE

    implicit none

    ! Arguments
    type(MESH_CONNECTIVITY),  dimension(ncells), intent(IN) :: Mesh
    integer(KIND = int_kind), dimension(:,:), pointer       :: Mesh_Collate_Vertex
    
    ! Local variables
    integer(KIND = int_kind) :: v

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (p_info%IOP) then
       ALLOCATE (Mesh_Collate_Vertex(nvc,ncells_tot))
    else
       ALLOCATE (Mesh_Collate_Vertex(nvc,0))
    end if
    
    do v = 1,nvc
       call PGSLib_COLLATE (Mesh_Collate_Vertex(v,:), Mesh%Ngbr_Vrtx_Orig(v))
    end do

    return

  END FUNCTION MESH_COLLATE_VERTEX

  FUNCTION VERTEX_COLLATE (Vertex)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed vertex into a single large vertex on IO PE
    !==================================================================
    use parallel_info_module, only: p_info
    use parameter_module,     only: nnodes_tot, nnodes
    use pgslib_module,        only: PGSLib_COLLATE

    implicit none

    ! Arguments
    type(VERTEX_DATA), dimension(nnodes), intent(IN) :: Vertex
    type(VERTEX_DATA), dimension(:), pointer         :: VERTEX_COLLATE
    
    ! Local variables
    integer(KIND = int_kind) :: n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (p_info%IOP) then
       ALLOCATE (Vertex_Collate(nnodes_tot))
    else
       ALLOCATE (Vertex_Collate(0))
    end if
    
    do n = 1,ndim
       call PGSLib_COLLATE (Vertex_Collate%Coord(n), Vertex%Coord(n))
    end do
    call PGSLib_COLLATE (Vertex_Collate%Rsum_Rvol, Vertex%Rsum_Rvol)

    return

  END FUNCTION VERTEX_COLLATE

  FUNCTION CELL_COLLATE (Cell)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed cell into a single large cell on IO PE
    !==================================================================
    use parallel_info_module, only: p_info
    use parameter_module,     only: ncells_tot, ncells

    implicit none

    ! Arguments
    type(CELL_GEOMETRY), dimension(ncells), intent(IN) :: Cell
    type(CELL_GEOMETRY), dimension(:), pointer         :: Cell_Collate

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (p_info%IOP) then
       ALLOCATE(Cell_Collate(ncells_tot))
    else
       ALLOCATE(Cell_Collate(0))
    end if
    
    call COLLATE(Cell_Collate, Cell)

    return

  END FUNCTION CELL_COLLATE

  SUBROUTINE COLLATE_CELL (Collated_Cell, Local_Cell)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed cell into a single large cell on IO PE
    !==================================================================
    use parameter_module,     only: ndim, nfc
    use pgslib_module,        only: PGSLib_COLLATE

    implicit none

    ! Arguments
    type(CELL_GEOMETRY), dimension(:), intent(IN   ) :: Local_Cell
    type(CELL_GEOMETRY), dimension(:), intent(  OUT) :: Collated_Cell
    
    ! Local variables
    integer(KIND = int_kind) :: f, n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call PGSLib_COLLATE(Collated_Cell%Volume, Local_Cell%Volume)
    do n = 1,ndim
       call PGSLib_COLLATE(Collated_Cell%Centroid(n), Local_Cell%Centroid(n))
    end do

    do f = 1, nfc
       call PGSLib_COLLATE(Collated_Cell%Face_Area(f), Local_Cell%Face_Area(f))
       call PGSLib_COLLATE(Collated_Cell%Halfwidth(f), Local_Cell%Halfwidth(f))
       do n = 1,ndim
          call PGSLib_COLLATE(Collated_Cell%Face_Normal(n,f), Local_Cell%Face_Normal(n,f))
          call PGSLib_COLLATE(Collated_Cell%Face_Centroid_L(n,f), Local_Cell%Face_Centroid_L(n,f))
       end do
    end do

    return

  END SUBROUTINE COLLATE_CELL

  SUBROUTINE PERMUTE_CELL (Permuted_Cell, Orig_Cell, Permuter, SCOPE)
    !==================================================================
    ! Purpose(s):
    !   Permute cell according to Permuter vector
    !==================================================================
    use parallel_scope

    ! Arguments
    type(CELL_GEOMETRY), dimension(:), intent(IN   ) :: Orig_Cell
    type(CELL_GEOMETRY), dimension(:), intent(  OUT) :: Permuted_Cell
    integer(int_kind), dimension(:), intent(IN) :: Permuter
    type (PL_SCOPE), OPTIONAL,    intent(IN   ) :: SCOPE

    ! Local variables
    type (PL_SCOPE) :: Desired_Scope

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Default scope is global
    if (PRESENT(SCOPE)) then
       Desired_Scope = SCOPE
    else
       Desired_Scope = GLOBAL_SCOPE
    end if

    if (DESIRED_SCOPE == GLOBAL_SCOPE) then
       call PERMUTE_CELL_GLOBAL(Permuted_Cell, Orig_Cell, Permuter)
    end if

    if (DESIRED_SCOPE == LOCAL_SCOPE) then
       call PERMUTE_CELL_LOCAL(Permuted_Cell, Orig_Cell, Permuter)
    end if

    return
  END SUBROUTINE PERMUTE_CELL

  SUBROUTINE PERMUTE_CELL_GLOBAL (Permuted_Cell, Orig_Cell, Permuter)
    !==================================================================
    ! Purpose(s):
    !   Permute cell according to Permuter vector
    !   Global version, so Permuter must be global indices
    !==================================================================
    use parameter_module,     only: ndim, nfc
    use pgslib_module,        only: PGSLib_Permute,    &
                                    PGSLIB_Deallocate_Trace, &
                                    PGSLib_GS_Trace

    implicit none

    ! Arguments
    type(CELL_GEOMETRY), dimension(:), intent(IN   ) :: Orig_Cell
    type(CELL_GEOMETRY), dimension(:), intent(  OUT) :: Permuted_Cell
    integer(int_kind), dimension(:), intent(IN) :: Permuter

    ! Local variables
    integer(KIND = int_kind) :: f, n
    type (PGSLib_GS_Trace), POINTER :: Cell_Trace => null()

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    NULLIFY(Cell_Trace)

    call PGSLib_PERMUTE(DEST   = Permuted_Cell%Volume,&
                        SOURCE = Orig_Cell%Volume,    &
                        INDEX  = Permuter,            &
                        TRACE  = Cell_Trace)
    do n = 1,ndim
       call PGSLib_PERMUTE(DEST   = Permuted_Cell%Centroid(n),&
                           SOURCE = Orig_Cell%Centroid(n),    &
                           INDEX  = Permuter,                 &
                           TRACE  = Cell_Trace)
    end do

    do f = 1, nfc
       call PGSLib_PERMUTE(DEST   = Permuted_Cell%Face_Area(f),&
                           SOURCE = Orig_Cell%Face_Area(f),    &
                           INDEX  = Permuter,                  &
                           TRACE  = Cell_Trace)
       call PGSLib_PERMUTE(DEST   = Permuted_Cell%Halfwidth(f),&
                           SOURCE = Orig_Cell%Halfwidth(f),    &
                           INDEX  = Permuter,                  &
                           TRACE  = Cell_Trace)
       do n = 1,ndim
          call PGSLib_PERMUTE(DEST   = Permuted_Cell%Face_Normal(n,f),&
                              SOURCE = Orig_Cell%Face_Normal(n,f),    &
                              INDEX  = Permuter,                      &
                              TRACE  = Cell_Trace)
          call PGSLib_PERMUTE(DEST   = Permuted_Cell%Face_Centroid_L(n,f),&
                              SOURCE = Orig_Cell%Face_Centroid_L(n,f),    &
                              INDEX  = Permuter,                        &
                              TRACE  = Cell_Trace)
       end do
    end do

    ! Done with the trace
    call PGSLib_DEALLOCATE_TRACE (Cell_Trace)

    return

  END SUBROUTINE PERMUTE_CELL_GLOBAL

  SUBROUTINE PERMUTE_CELL_LOCAL (Permuted_Cell, Orig_Cell, Permuter)
    !==================================================================
    ! Purpose(s):
    !   Permute cell according to Permuter vector
    !   Local version, so Permuter must be local indices
    !==================================================================

    ! Arguments
    type(CELL_GEOMETRY), dimension(:), intent(IN   ) :: Orig_Cell
    type(CELL_GEOMETRY), dimension(:), intent(  OUT) :: Permuted_Cell
    integer(int_kind), dimension(:), intent(IN) :: Permuter

    ! Local variables
    integer(int_kind) :: cell

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do cell = 1, SIZE(Permuter)
       Permuted_Cell(Permuter(cell)) = Orig_Cell(Cell)
    end do

    return

  END SUBROUTINE PERMUTE_CELL_LOCAL

  FUNCTION MESH_DISTRIBUTE (Mesh_Tot)
    !====================================================================
    ! Purpose(s):
    !
    !   Distribute the components of the mesh data structure. If the
    !   mesh data structure changes, this routine must be updated too.
    !
    !====================================================================
    use kind_module,      only: int_kind
    use parameter_module, only: ncells, nfc, nvc
    use pgslib_module,    only: PGSLib_DIST

    implicit none

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(:),     intent(IN) :: Mesh_Tot
    type(MESH_CONNECTIVITY), dimension(ncells)            :: Mesh_Distribute

    ! Local Variables
    integer(KIND = int_kind) :: f, v

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Distribute the ngbr_cell structure
    DIST_NGBR_CELL: do f = 1, nfc

#ifdef DEBUG_MESH_MODULE
       print *, "DIST_NGBR_CELL: ", p_info%thisPE, f
#endif
       call PGSLib_DIST (Mesh_Distribute%Ngbr_Cell(f), Mesh_Tot%Ngbr_Cell(f))

    end do DIST_NGBR_CELL

    ! Distribute the ngbr_face structure
    DIST_NGBR_FACE: do f = 1, nfc

       call PGSLib_DIST (Mesh_Distribute%Ngbr_Face(f), Mesh_Tot%Ngbr_Face(f))

    end do DIST_NGBR_FACE

    ! Distribute the ngbr_vrtx structure
    DIST_NGBR_VRTX: do v = 1, nvc

       call PGSLib_DIST (Mesh_Distribute%Ngbr_Vrtx(v), Mesh_Tot%Ngbr_Vrtx(v))

    end do DIST_NGBR_VRTX

    ! Distribute the cell shape integer
    call PGSLib_DIST (Mesh_Distribute%Cell_shape, Mesh_Tot%Cell_shape)
    
    ! Distribute the cell block ID data
    call PGSLib_DIST (Mesh_Distribute%CBlockID, Mesh_Tot%CBlockID)

    return

  END FUNCTION MESH_DISTRIBUTE

  FUNCTION VERTEX_DISTRIBUTE (Vertex_Tot)
    !====================================================================
    ! Purpose(s):
    !
    !   Distribute the components of the vertex data structure. If that
    !   data structure changes this routine must be updated too.
    !
    !====================================================================
    use kind_module,          only: int_kind
    use parameter_module,     only: ndim, nnodes
    use pgslib_module,        only: PGSLib_DIST

    implicit none

    ! Arguments
    type(VERTEX_DATA), dimension(:),     intent(IN) :: Vertex_Tot
    type(VERTEX_DATA), dimension(nnodes)            :: Vertex_Distribute

    ! Local Variables
    integer(KIND = int_kind) :: n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Distribute the vertex coordinates
    DIST_VERTEX_COORD: do n = 1,ndim

       call PGSLib_DIST (Vertex_Distribute%Coord(n), Vertex_Tot%Coord(n))

    end do DIST_VERTEX_COORD

    call PGSLib_DIST (Vertex_Distribute%Rsum_Rvol, Vertex_Tot%Rsum_Rvol)

    return

  END FUNCTION VERTEX_DISTRIBUTE
  
  FUNCTION CELL_DISTRIBUTE (Cell_Tot)
    !====================================================================
    ! Purpose(s):
    !
    !   Distribute the components of the cell data structure. If that
    !   data structure changes this routine must be updated too.
    !
    !====================================================================
    use kind_module,          only: int_kind
    use parameter_module,     only: ncells
    use pgslib_module,        only: PGSLib_DIST

    implicit none

    ! Arguments
    type(CELL_GEOMETRY), dimension(:), intent(IN) :: Cell_Tot
    type(CELL_GEOMETRY), dimension(ncells)        :: Cell_Distribute

    ! Local Variables
    integer(KIND = int_kind) :: n, f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Distribute the cell coordinates
    DIST_CELL_FACE_NORMAL: do f = 1,nfc

       do n = 1,ndim
          call PGSLib_DIST (Cell_Distribute%Face_Normal(n,f), Cell_Tot%Face_Normal(n,f))
       end do

    end do DIST_CELL_FACE_NORMAL

    DIST_CELL_FACE_AREA: do f = 1,nfc

       call PGSLib_DIST (Cell_Distribute%Face_Area(f), Cell_Tot%Face_Area(f))

    end do DIST_CELL_FACE_AREA

    DIST_CELL_FACE_CENTROID: do f = 1,nfc

       do n = 1,ndim
          call PGSLib_DIST (Cell_Distribute%Face_Centroid_L(n,f), Cell_Tot%Face_Centroid_L(n,f))
       end do

    end do DIST_CELL_FACE_CENTROID


    DIST_CELL_CENTROID: do n = 1,ndim

       call PGSLib_DIST (Cell_Distribute%Centroid(n), Cell_Tot%Centroid(n))

    end do DIST_CELL_CENTROID

    DIST_CELL_HALFWIDTH: do f = 1,nfc

       call PGSLib_DIST (Cell_Distribute%Halfwidth(f), Cell_Tot%Halfwidth(f))

    end do DIST_CELL_HALFWIDTH

    call PGSLib_DIST (Cell_Distribute%Volume, Cell_Tot%Volume)

    return

  END FUNCTION CELL_DISTRIBUTE

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_MESH_DATA
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 23 Apr 2005
 !!
 !! This subroutine reads the mesh data from a restart file opened (and pre-
 !! positioned) on UNIT, and initializes certain components of the module
 !! structures MESH and VERTEX with this data (properly distributed).  Only
 !! the NGBR_VRTX and CBLOCKID components of MESH, and the COORD component
 !! of VERTEX are initialized.  VERSION is the version number of the restart
 !! file format.
 !!
 !! It is assumed that the structures MESH and VERTEX have already been
 !! suitably allocated.  Note that the data read here is *not* permuted
 !! unlike the subsequent cell and node-based data read from the restart
 !! file.
 !!

  subroutine read_mesh_data (unit, version)

    use restart_utilities, only: read_var, read_dist_array

    integer, intent(in) :: unit, version

    integer :: n

    do n = 1, nvc
      call read_dist_array (unit, mesh%ngbr_vrtx(n), errmsg='READ_MESH_DATA: error reading VERTEX records')
    end do

    call read_var (unit, n, 'READ_MESH_DATA: error reading NCBLOCK record.')
    if (n /= 0) then
      mesh_has_cblockid_data = .true.
      call read_dist_array (unit, mesh%cblockid, errmsg='READ_MESH_DATA: error reading CBLOCKID record')
    else
      mesh_has_cblockid_data = .false.
    end if

    do n = 1, ndim
      call read_dist_array (unit, vertex%coord(n), errmsg='READ_MESH_DATA: error reading COORD records')
    end do

  end subroutine read_mesh_data

END MODULE MESH_MODULE
