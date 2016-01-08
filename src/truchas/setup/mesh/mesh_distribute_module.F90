!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MESH_DISTRIBUTE_MODULE
  !======================================================================
  ! Purpose(s):
  !
  !   Distribute mesh to all the processors.
  !
  ! Contains: DISTRIBUTE_MESH_VERTEX
  !           DISTRIBUTE_MESH
  !           DISTRIBUTE_VERTEX
  !           COLLATE_MESH
  !
  ! Author(s): Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !
  !======================================================================
  use kinds, only: r8
  use parameter_module, only: max_domains
  use truchas_logging_services
  implicit none
  private

  ! Public Variables
  public :: NCells_List, NNodes_List, number_domains

  ! Public Procedures
  public :: DISTRIBUTE_MESH_VERTEX, COLLATE_MESH, &
            DISTRIBUTE_MESH, DISTRIBUTE_VERTEX

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Global data provided by this module
  integer :: number_domains
  integer, dimension(max_domains) :: NCells_List = -1, &
                                                      NNodes_List = -1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE DISTRIBUTE_MESH_VERTEX (Mesh_Tot, Vertex_Tot, Mesh, Vertex)
    !=======================================================================
    ! Purpose(s):
    !
    !   Distribute a mesh from the IO root pe to all other PEs.
    !   Assume that decomposition is already determined.
    !
    !=======================================================================
    use mesh_module,          only: MESH_CONNECTIVITY, VERTEX_DATA
    use parallel_info_module, only: p_info
    use parameter_module,     only: ncells, ncells_tot, nnodes, nnodes_tot
#ifdef USE_PGSLIB
    use pgslib_module,        only: PGSLIB_COLLATE
#endif

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(ncells_tot), intent(IN)  :: Mesh_Tot
    type(VERTEX_DATA)      , dimension(nnodes_tot), intent(IN)  :: Vertex_Tot
    type(VERTEX_DATA)      , dimension(nnodes),     intent(OUT) :: Vertex
    type(MESH_CONNECTIVITY), dimension(ncells),     intent(OUT) :: Mesh

    ! Local Variables
    logical :: fatal
    character(128), dimension(4) :: error_string
    integer :: ncells_tot_test, nnodes_tot_test
    integer, allocatable, dimension(:) :: Mesh_Extents, Vertex_Extents

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! We need the extents of the arrays to be distributed. Each PE knows its
    ! local extent so we need to collate these to IO PE.
    ALLOCATE (Mesh_Extents(p_info%nPE), Vertex_Extents(p_info%nPE))

#ifdef USE_PGSLIB
    call PGSLIB_COLLATE (Mesh_Extents, ncells)
    call PGSLIB_COLLATE (Vertex_Extents, nnodes)
#else
    Mesh_Extents = ncells
    Vertex_Extents = nnodes
#endif

    ! Check that the sizes add up properly
    fatal = .false.
    CHECK_DIMS: if (p_info%IOP) then

       ncells_tot_test = SUM(Mesh_Extents)
       nnodes_tot_test = SUM(Vertex_Extents)
       fatal = (ncells_tot_test /= ncells_tot) .OR. (nnodes_tot_test /= nnodes_tot)

    end if CHECK_DIMS

    call TLS_fatal_if_any (fatal, 'DISTRIBUTE_MESH_VERTEX: inconsistent array sizes')

    ! Assuming the sizes are correct, now we can distribute the mesh and vertex.
    ! These are compound data structures, so we have to make sure we get each component.
    call DISTRIBUTE_MESH (Mesh, Mesh_Tot, Mesh_Extents)

    call DISTRIBUTE_VERTEX (Vertex, Vertex_Tot, Vertex_Extents)

    deallocate (Vertex_Extents, Mesh_Extents)

  END SUBROUTINE DISTRIBUTE_MESH_VERTEX

  SUBROUTINE DISTRIBUTE_MESH (Mesh, Mesh_Tot, Mesh_Extents)
    !====================================================================
    ! Purpose(s):
    !
    !   Distribute the components of the mesh data structure. If the
    !   mesh data structure changes, this routine must be updated too.
    !
    !====================================================================
    use mesh_module,          only: MESH_CONNECTIVITY
    use parameter_module,     only: ncells, ncells_tot, nfc, nvc
#ifdef USE_PGSLIB
    use pgslib_module,        only: PGSLIB_DIST
#endif

    ! Argument List
    type(MESH_CONNECTIVITY), dimension(ncells_tot), intent(IN) :: Mesh_Tot
    integer, dimension(:), optional, intent(IN) :: Mesh_Extents
    type(MESH_CONNECTIVITY), dimension(ncells), intent(OUT) :: Mesh

    ! Local Variables
    integer :: i
    integer, dimension(ncells)     :: Local_Temp_Int
    integer, dimension(ncells_tot) :: Total_Temp_Int

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Distribute the ngbr_cell structure
    DIST_NGBR_CELL: do i = 1, nfc

       Total_Temp_Int = Mesh_Tot%Ngbr_Cell(i)
#ifdef USE_PGSLIB
       call PGSLIB_DIST (Local_Temp_Int, Total_Temp_Int, Mesh_Extents)
#else
       Local_Temp_Int = Total_Temp_Int
#endif
       Mesh%Ngbr_Cell(i) = Local_Temp_Int

    end do DIST_NGBR_CELL

    ! Distribute the ngbr_face structure
    DIST_NGBR_FACE: do i = 1, nfc

       Total_Temp_Int = Mesh_Tot%Ngbr_Face(i)
#ifdef USE_PGSLIB
       call PGSLIB_DIST (Local_Temp_Int, Total_Temp_Int, Mesh_Extents)
#else
       Local_Temp_Int = Total_Temp_Int
#endif
       Mesh%Ngbr_Face(i) = Local_Temp_Int

    end do DIST_NGBR_FACE

    ! Distribute the ngbr_vrtx structure
    DIST_NGBR_VRTX: do i = 1, nvc

       Total_Temp_Int = Mesh_Tot%Ngbr_Vrtx(i)
#ifdef USE_PGSLIB
       call PGSLIB_DIST (Local_Temp_Int, Total_Temp_Int, Mesh_Extents)
#else
       Local_Temp_Int = Total_Temp_Int
#endif
       Mesh%Ngbr_Vrtx(i) = Local_Temp_Int

    end do DIST_NGBR_VRTX

  END SUBROUTINE DISTRIBUTE_MESH

  SUBROUTINE DISTRIBUTE_VERTEX (Vertex, Vertex_Tot, Vertex_Extents)
    !====================================================================
    ! Purpose(s):
    !
    !   Distribute the components of the vertex data structure. If that
    !   data structure changes this routine must be updated too.
    !
    !====================================================================
    use mesh_module,          only: VERTEX_DATA
    use parameter_module,     only: ndim, nnodes, nnodes_tot 
#ifdef USE_PGSLIB
    use pgslib_module,        only: PGSLIB_DIST
#endif

    ! Arguments
    type(VERTEX_DATA), dimension(nnodes_tot), intent(IN) :: Vertex_Tot
    integer, dimension(:), optional, intent(IN) :: Vertex_Extents
    type(VERTEX_DATA), dimension(nnodes), intent(OUT) :: Vertex

    ! Local Variables
    integer :: n
    real(r8), dimension(nnodes)     :: Local_Temp_Real
    real(r8), dimension(nnodes_tot) :: Total_Temp_Real

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Distribute the vertex coordinates
    DIST_VERTEX_COORD: do n = 1,ndim

       Total_Temp_Real = Vertex_Tot%Coord(n)
#ifdef USE_PGSLIB
       call PGSLIB_DIST (Local_Temp_Real, Total_Temp_Real, Vertex_Extents)
#else
       Local_Temp_Real = Total_Temp_Real
#endif
       Vertex%Coord(n) = Local_Temp_Real

    end do DIST_VERTEX_COORD

    Total_Temp_Real = Vertex_Tot%Rsum_rvol
#ifdef USE_PGSLIB
    call PGSLIB_DIST (Local_Temp_Real, Total_Temp_Real, Vertex_Extents)
#else
    Local_Temp_Real = Total_Temp_Real
#endif
    Vertex%Rsum_rvol = Local_Temp_Real

  END SUBROUTINE DISTRIBUTE_VERTEX
  
  SUBROUTINE COLLATE_MESH (Mesh_Tot, Mesh)
    !====================================================================
    ! Purpose(s):
    !
    !   Collate the components of the mesh data structure. If the
    !   mesh data structure changes, this routine must be updated too.
    !
    !====================================================================
    use mesh_module,          only: MESH_CONNECTIVITY
    use parameter_module,     only: nfc, nvc
#ifdef USE_PGSLIB
    use pgslib_module,        only: PGSLIB_COLLATE
#endif

    ! Argument List
    type(MESH_CONNECTIVITY), dimension(:), intent(IN)  :: Mesh
    type(MESH_CONNECTIVITY), dimension(:), intent(OUT) :: Mesh_Tot

    ! Local Variables
    integer :: i
    integer, dimension(SIZE(Mesh,1))     :: Local_Temp_Int
    integer, dimension(SIZE(Mesh_Tot,1)) :: Total_Temp_Int

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Collate the ngbr_cell structure
    COLLATE_NGBR_CELL: do i = 1, nfc

       Local_Temp_Int = Mesh%Ngbr_cell(i)
#ifdef USE_PGSLIB
       call PGSLIB_COLLATE (Total_Temp_Int, Local_Temp_Int)
#else
       Total_Temp_Int = Local_Temp_Int 
#endif
       Mesh_Tot%Ngbr_cell(i) = Total_Temp_Int

    end do COLLATE_NGBR_CELL

    ! Collate the ngbr_face structure
    COLLATE_NGBR_FACE: do i = 1, nfc

       Local_Temp_Int = Mesh%Ngbr_face(i)
#ifdef USE_PGSLIB
       call PGSLIB_COLLATE (Total_Temp_Int, Local_Temp_Int)
#else
       Total_Temp_Int = Local_Temp_Int 
#endif
       Mesh_Tot%Ngbr_face(i) = Total_Temp_Int

    end do COLLATE_NGBR_FACE

    ! Collate the ngbr_vrtx structure
    COLLATE_NGBR_VRTX: do i = 1, nvc

       Local_Temp_Int = Mesh%Ngbr_vrtx(i)
#ifdef USE_PGSLIB
       call PGSLIB_COLLATE (Total_Temp_Int, Local_Temp_Int)
#else
       Total_Temp_Int = Local_Temp_Int 
#endif
       Mesh_Tot%Ngbr_vrtx(i) = Total_Temp_Int

    end do COLLATE_NGBR_VRTX

  END SUBROUTINE COLLATE_MESH

END MODULE MESH_DISTRIBUTE_MODULE
