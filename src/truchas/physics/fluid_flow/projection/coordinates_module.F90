MODULE COORDINATES_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Gather and store cell-vertex and neighbor cell-centroid coordinates.
  !
  ! Contains: GATHER_COORDINATES
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  ! Public Variables
  public :: Cell_Ngbr_Coord, Vertex_Coord

  ! Public Subroutines
  public :: GATHER_COORDINATES

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Gathered Coordinates
  ! Cell-Vertex Coordinates (Vertex_Coord(ndim,nvc,ncells))
  real(r8), dimension(:,:,:), pointer, save :: Vertex_Coord

  ! Neighbor Cell-Centroid Coordinates (Cell_Ngbr_Coord(ndim,nfc,ncells))
  real(r8), dimension(:,:,:), pointer, save :: Cell_Ngbr_Coord

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE GATHER_COORDINATES ()
    !=======================================================================
    ! Purpose(s):
    !   Gather cell-vertex and neighboring cell-centroid coordinates;
    !   store the data in arrays Cell_Ngbr_Coord(ndim,nfc,ncells) and
    !   Vertex_Coord(ndim,nvc,ncells).
    !=======================================================================
    use gs_module,        only: EE_GATHER, EN_GATHER
    use mesh_module,      only: Cell, Vertex
    use parameter_module, only: ndim

    ! Local Variables
    integer :: n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Coordinate Arrays
    do n = 1, ndim
       ! Cell-Vertex Coordinates
       call EN_GATHER (Vertex_Coord(n,:,:),Vertex%Coord(n))

       ! Neighbor Cell-Centroid Coordinates
       call EE_GATHER (Cell_Ngbr_Coord(n,:,:), Cell%Centroid(n))
    end do

  END SUBROUTINE GATHER_COORDINATES

END MODULE COORDINATES_MODULE
