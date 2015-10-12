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
  use kinds, only: r8
  use parameter_module, only: ndim, nfc, nvc, nvf, nec
  use var_vector_types
  implicit none
  private

  ! Public Variables and Types
  public :: MESH_CONNECTIVITY, CELL_EDGE, &
            CELL_GEOMETRY, VERTEX_DATA, BOUNDARY, Vrtx_Bdy, &
            Cell, Mesh, Vertex, & 
            CELL_TET, CELL_PYRAMID, CELL_PRISM, CELL_HEX, &
            GAP_ELEMENT_1,GAP_ELEMENT_3,GAP_ELEMENT_5

  public :: Vertex_Ngbr_All, Vertex_Ngbr_All_Orig
  public :: MESH_COLLATE_VERTEX, VERTEX_COLLATE

  ! Public Subroutines
  public :: Initialize_Face_Bit_Mask,        &
            Set_Face_Neighbor,               &
            Clear_Face_Neighbor,             &
            Is_Face_Ngbr

  public :: COLLATE_CELL, PERMUTE_CELL
  public :: read_mesh_data

  ! Public Operators
  public :: operator(.DISTRIBUTE.)

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Various mesh flags.
  logical, public, save :: orthogonal_mesh

  !! Parameters to identify the state of cells, faces and nodes.
  ! Degenerate face and node parameters.
  integer, parameter, public :: EXTERNAL_FACE   = 0
  integer, parameter, public :: DEGENERATE_FACE = -HUGE(1)
  integer, parameter, public :: DEGENERATE_NODE = -1
  ! Define unknown, default, vertex (node), cell, face and PE numbers
  integer, parameter :: DEFAULT_CELL   = -1
  integer, parameter :: DEFAULT_VERTEX = -1
  integer, parameter :: DEFAULT_FACE   = -1
  integer, parameter :: DEFAULT_PE     = -1

  ! Mesh diagnostics variables.
  integer, public, save :: degenerate_points, degenerate_lines, &
                                            triangle_faces, quad_faces
  real(r8),   public, save :: volume_min, volume_max

  ! Face-Vertex and Vertex-Face mappings.
  integer, public, save, dimension(nfc,nvf)  :: Face_Vrtx
  data Face_Vrtx /3,1,4,2,3,8,  4,2,1,3,2,5,  8,6,5,7,1,6,  7,5,8,6,4,7/
  integer, public, save, dimension(nvc,ndim) :: Vrtx_Face
  data Vrtx_Face /2,2,1,1,2,2,1,1,  3,4,4,3,3,4,4,3,  5,5,5,5,6,6,6,6/
  
  ! MESH_CONNECTIVITY structure
  Type MESH_CONNECTIVITY
     integer :: Ngbr_Cell(nfc)      = DEFAULT_CELL   ! cell number across each face
     integer :: Ngbr_Cell_PE(nfc)   = DEFAULT_PE     ! holding data for Ngbr_Cell
     integer :: Ngbr_Cell_Orig(nfc) = DEFAULT_CELL   ! original cell numbers across each face
     integer :: Ngbr_Face(nfc)      = DEFAULT_FACE   ! face number of cell across each face
     integer :: Ngbr_Vrtx(nvc)      = DEFAULT_VERTEX ! cell vertex numbers
     integer :: Ngbr_Vrtx_PE(nvc)   = DEFAULT_PE     ! holding data for Ngbr_Vrtx
     integer :: Ngbr_Vrtx_Orig(nvc) = DEFAULT_VERTEX ! original vertex numbers
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
  integer, POINTER,      &
                     PUBLIC,       &
                     SAVE,         &
                     dimension(:) :: UnPermute_Mesh_Vector => NULL()
  integer, POINTER,      &
                     PUBLIC,       &
                     SAVE,         &
                     dimension(:) :: Permute_Mesh_Vector => NULL()
  logical, SAVE,         &
                     PUBLIC       :: UnPermute_Mesh_Initialized = .FALSE.
  integer, POINTER,      &
                     PUBLIC,       &
                     SAVE,         &
                     dimension(:) :: UnPermute_Vertex_Vector => NULL()
  integer, POINTER,      &
                     PUBLIC,       &
                     SAVE,         &
                     dimension(:) :: Permute_Vertex_Vector => NULL()
  logical, SAVE,         &
                     PUBLIC       :: UnPermute_Vertex_Initialized = .FALSE.

  ! Define the edges surrounding a hex cell
  integer :: Cell_Edge(2,12)
  data Cell_Edge/1,2, 2,3, 3,4, 4,1, 2,6, 3,7, 4,8, 1,5, 5,6, 6,7, 7,8, 8,5/

  ! define CELL_GEOMETRY derived type
  Type CELL_GEOMETRY
     real(r8), dimension(ndim)     :: Centroid        ! cell centroid (physical coordinates)
!    real(r8), dimension(ndim)     :: Centroid_L      ! cell centroid (logical coordinates) (place holder)
     real(r8)                      :: Volume          ! cell volume
     real(r8), dimension(nfc)      :: Halfwidth       ! cell halfwidth
     real(r8), dimension(ndim,nfc) :: Face_Centroid   ! face centroid (physical coordinates)
     real(r8), dimension(ndim,nfc) :: Face_Centroid_L ! face centroid (logical coordinates)
     real(r8), dimension(nfc)      :: Face_Area       ! face area
     real(r8), dimension(ndim,nfc) :: Face_Normal     ! face unit normal
  End Type CELL_GEOMETRY

  ! Declare instance of CELL_GEOMETRY type.
  type(CELL_GEOMETRY), dimension(:), SAVE, pointer :: Cell => NULL()

  ! Vertex connectivity
  Type(int_var_vector), dimension(:), POINTER :: Vertex_Ngbr_All => NULL()
  Type(int_var_vector), dimension(:), POINTER :: Vertex_Ngbr_All_Orig => NULL()

  ! VERTEX_DATA structure
  Type VERTEX_DATA
     real(r8) :: Coord(ndim) = 0.0_r8  ! vertex coordinates
     real(r8) :: Rsum_rvol   = 0.0_r8  ! reciprocal sum of inverse volumes around vertices
  End Type VERTEX_DATA

  ! Ghost cells for the vertex data
  ! **DANGER** If coord data changes these need to be nullified.
  ! Define type for boundary data
  type BOUNDARY
    real(r8), dimension(:), pointer :: Data => NULL()
  end type BOUNDARY

  type(BOUNDARY), dimension(ndim), SAVE :: Vrtx_Bdy

  ! Define instance of VERTEX_DATA type.
  type(VERTEX_DATA), dimension(:), pointer, save :: Vertex => NULL()

  ! Bit mask for identifying face neighbors
  integer, dimension(:), pointer :: Face_Bit_Mask => NULL()
  logical, SAVE                  :: Face_Bit_Mask_Initialized = .FALSE.

  interface operator(.DISTRIBUTE.)
     module procedure MESH_DISTRIBUTE
     module procedure VERTEX_DISTRIBUTE
     module procedure CELL_DISTRIBUTE
  end interface
  
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

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

  END SUBROUTINE Initialize_Face_Bit_Mask
       
  !==================================================================
  ! Purpose(s):
  !   Set the neighbor flag to indicate that it is a neighbor
  !   of Face.
  !==================================================================

  SUBROUTINE Set_Face_Neighbor(Face_Bits, Face)
    integer, intent(INOUT) :: Face_Bits
    integer, intent(IN   ) :: Face
    Face_Bits = IOR(Face_Bits, Face_Bit_Mask(Face))
  END SUBROUTINE Set_Face_Neighbor
    
    !==================================================================
    ! Purpose(s):
    !   Clear the neighbor flag which indicates that it is a neighbor
    !   of Face.
    !==================================================================

  SUBROUTINE Clear_Face_Neighbor(Face_Bits, Face)
    integer, intent(INOUT) :: Face_Bits
    integer, intent(IN   ) :: Face
    Face_Bits = IAND(Face_Bits, NOT(Face_Bit_Mask(Face)))
  END SUBROUTINE Clear_Face_Neighbor
      
    !==================================================================
    ! Purpose(s):
    !   Test the flag of an integer to see if that neighbor is also
    !   a neighbor of Face.
    !   Set the neighbor flag to indicate that it is a neighbor
    !   of Face.
    !==================================================================

  FUNCTION Is_Face_Ngbr(Face_Bits, Face)
    integer, intent(IN   ) :: Face_Bits
    integer, intent(IN   ) :: Face
    logical                :: Is_Face_Ngbr
    Is_Face_Ngbr = (IAND(Face_Bits, Face_Bit_Mask(Face)) /= 0)
  END FUNCTION Is_Face_Ngbr
      
  FUNCTION MESH_COLLATE_VERTEX (Mesh)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed mesh into a single large mesh on IO PE
    !==================================================================
    use parallel_info_module, only: p_info
    use parameter_module,     only: ncells_tot, nvc, ncells
    use pgslib_module,        only: PGSLib_COLLATE

    ! Arguments
    type(MESH_CONNECTIVITY),  dimension(ncells), intent(IN) :: Mesh
    integer, dimension(:,:), pointer       :: Mesh_Collate_Vertex
    
    ! Local variables
    integer :: v

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (p_info%IOP) then
       ALLOCATE (Mesh_Collate_Vertex(nvc,ncells_tot))
    else
       ALLOCATE (Mesh_Collate_Vertex(nvc,0))
    end if
    
    do v = 1,nvc
       call PGSLib_COLLATE (Mesh_Collate_Vertex(v,:), Mesh%Ngbr_Vrtx_Orig(v))
    end do

  END FUNCTION MESH_COLLATE_VERTEX

  FUNCTION VERTEX_COLLATE (Vertex)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed vertex into a single large vertex on IO PE
    !==================================================================
    use parallel_info_module, only: p_info
    use parameter_module,     only: nnodes_tot, nnodes
    use pgslib_module,        only: PGSLib_COLLATE

    ! Arguments
    type(VERTEX_DATA), dimension(nnodes), intent(IN) :: Vertex
    type(VERTEX_DATA), dimension(:), pointer         :: VERTEX_COLLATE
    
    ! Local variables
    integer :: n

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

  END FUNCTION VERTEX_COLLATE

  FUNCTION CELL_COLLATE (Cell)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed cell into a single large cell on IO PE
    !==================================================================
    use parallel_info_module, only: p_info
    use parameter_module,     only: ncells_tot, ncells

    ! Arguments
    type(CELL_GEOMETRY), dimension(ncells), intent(IN) :: Cell
    type(CELL_GEOMETRY), dimension(:), pointer         :: Cell_Collate

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (p_info%IOP) then
       ALLOCATE(Cell_Collate(ncells_tot))
    else
       ALLOCATE(Cell_Collate(0))
    end if
    
    call COLLATE_CELL(Cell_Collate, Cell)

  END FUNCTION CELL_COLLATE

  SUBROUTINE COLLATE_CELL (Collated_Cell, Local_Cell)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed cell into a single large cell on IO PE
    !==================================================================
    use parameter_module,     only: ndim, nfc
    use pgslib_module,        only: PGSLib_COLLATE

    ! Arguments
    type(CELL_GEOMETRY), dimension(:), intent(IN   ) :: Local_Cell
    type(CELL_GEOMETRY), dimension(:), intent(  OUT) :: Collated_Cell
    
    ! Local variables
    integer :: f, n

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
    integer, dimension(:), intent(IN) :: Permuter
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

    ! Arguments
    type(CELL_GEOMETRY), dimension(:), intent(IN   ) :: Orig_Cell
    type(CELL_GEOMETRY), dimension(:), intent(  OUT) :: Permuted_Cell
    integer, dimension(:), intent(IN) :: Permuter

    ! Local variables
    integer :: f, n
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
    integer, dimension(:), intent(IN) :: Permuter

    ! Local variables
    integer :: cell

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do cell = 1, SIZE(Permuter)
       Permuted_Cell(Permuter(cell)) = Orig_Cell(Cell)
    end do

  END SUBROUTINE PERMUTE_CELL_LOCAL

  FUNCTION MESH_DISTRIBUTE (Mesh_Tot)
    !====================================================================
    ! Purpose(s):
    !
    !   Distribute the components of the mesh data structure. If the
    !   mesh data structure changes, this routine must be updated too.
    !
    !====================================================================
    use parameter_module, only: ncells, nfc, nvc
    use pgslib_module,    only: PGSLib_DIST

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(:),     intent(IN) :: Mesh_Tot
    type(MESH_CONNECTIVITY), dimension(ncells)            :: Mesh_Distribute

    ! Local Variables
    integer :: f, v

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

  END FUNCTION MESH_DISTRIBUTE

  FUNCTION VERTEX_DISTRIBUTE (Vertex_Tot)
    !====================================================================
    ! Purpose(s):
    !
    !   Distribute the components of the vertex data structure. If that
    !   data structure changes this routine must be updated too.
    !
    !====================================================================
    use parameter_module,     only: ndim, nnodes
    use pgslib_module,        only: PGSLib_DIST

    ! Arguments
    type(VERTEX_DATA), dimension(:),     intent(IN) :: Vertex_Tot
    type(VERTEX_DATA), dimension(nnodes)            :: Vertex_Distribute

    ! Local Variables
    integer :: n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Distribute the vertex coordinates
    DIST_VERTEX_COORD: do n = 1,ndim

       call PGSLib_DIST (Vertex_Distribute%Coord(n), Vertex_Tot%Coord(n))

    end do DIST_VERTEX_COORD

    call PGSLib_DIST (Vertex_Distribute%Rsum_Rvol, Vertex_Tot%Rsum_Rvol)

  END FUNCTION VERTEX_DISTRIBUTE
  
  FUNCTION CELL_DISTRIBUTE (Cell_Tot)
    !====================================================================
    ! Purpose(s):
    !
    !   Distribute the components of the cell data structure. If that
    !   data structure changes this routine must be updated too.
    !
    !====================================================================
    use parameter_module,     only: ncells
    use pgslib_module,        only: PGSLib_DIST

    ! Arguments
    type(CELL_GEOMETRY), dimension(:), intent(IN) :: Cell_Tot
    type(CELL_GEOMETRY), dimension(ncells)        :: Cell_Distribute

    ! Local Variables
    integer :: n, f

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
