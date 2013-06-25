! set to 1 to turn on local debugging code
#define DEBUG_FACE_VECTORS 0

MODULE CELL_GEOMETRY_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures necessary to charcterize a cell's geometry.
  !   Specifically, this module contains routines which initialize
  !   the cell's volume, centroid, face areas, face unit normals,
  !   face centroids, face sweep direction, jacobian, and width.
  !
  !   Public Interface:
  !
  !     * call GET_CELL_GEOMETRY (Cell, Mesh, Vertex)
  !
  !         Returns the Cell derived type initialized with
  !         all appropriate cell geometry information.
  !
  !     * call FACE_CENTROID_PHYS ()
  !
  !         stores the face centroid physical coordinates
  !         in cell(c)%face_centroid(ndim,nfc)
  !
  ! Contains: GET_CELL_GEOMETRY
  !           CELL_CENTROID
  !           CELL_TYPE
  !           CELL_VOLUME
  !           CELL_WIDTH
  !           FACE_AREA
  !           FACE_CENTROID_PHYSICAL
  !           FACE_CENTROID_LOGICAL
  !           JACOBIAN
  !
  ! Author(s): Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
  use truchas_logging_services
  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: GET_CELL_GEOMETRY

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE GET_CELL_GEOMETRY ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Get the cell geometry for each cell, which is done by filling up all
    !   components of the Cell derived type. The following routines are called:
    !      CELL_VOLUME   - Sets Cell%Volume  
    !      CELL_CENTROID - Sets Cell%Centroid
    !      FACE_AREA     - Sets Cell%Area, Cell%Face_Normal
    !      JACOBIAN      - Checks for orthogonality
    !      FACE_CENTROID_LOGICAL
    !      FACE_CENTROID_PHYSICAL
    !
    !=======================================================================
    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Vertex
    use parameter_module,     only: ndim
    use pgslib_module,        only: PGSLib_Global_MINVAL, PGSLib_Global_MAXVAL
    use truchas_logging_services
    implicit none

    ! local variables
    integer (int_kind) :: n
    real (real_kind)   :: mincoord
    real (real_kind)   :: maxcoord
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Inform the user of cell geometry calculation.
    call TLS_info ('')
    call TLS_info ('Computing cell geometry ... ', advance=.false.)

    ! Compute and place cell geometry. in the Cell derived type.
    call GAP_ELEMENT_TYPE ()
    call CELL_VOLUME ()
    call CELL_CENTROID ()
    call FACE_AREA ()
    call FACE_CENTROID_LOGICAL ()
    call FACE_CENTROID_PHYSICAL ()
    call CELL_WIDTH ()
    call CELL_TYPE ()

    ! Compute and use cell Jacobians to search for orthogonal cells.
    call JACOBIAN ()

    ! Inform the user of successful geometry computation.
    call TLS_info ('done.')

    ! optionally report some information about the coordinate extents
    if (TLS_VERB_NORMAL <= TLS_verbosity) then
       call TLS_info ('')
       call TLS_info ('                    Min Coord        Max Coord')
       call TLS_info ('                    ---------        ---------')
       
       do n = 1,ndim
          mincoord = PGSLib_Global_MINVAL (Vertex%Coord(n))
          maxcoord = PGSLib_Global_MAXVAL (Vertex%Coord(n))
          write (message, 10) mincoord, maxcoord
10        format (18x,1pe11.4,6x,1pe11.4)
          call TLS_info (message)
       end do
    end if

    return

  END SUBROUTINE GET_CELL_GEOMETRY

  SUBROUTINE FACE_CENTROID_PHYSICAL ()
    !=======================================================================
    ! Purpose:
    !
    !   convert the logical face centroid coordinates (Cell%Face_centroid_L)
    !   to physical coordinates and store in cell(:)%face_centroid(ndim,nfc)
    !
    !=======================================================================
    use gs_module,        only: EN_GATHER
    use kind_module,      only: int_kind, real_kind
    use linear_module,    only: LINEAR_PROP
    use mesh_module,      only: Cell, Vertex, Vrtx_Bdy
    use parameter_module, only: ncells, nfc, ndim, nvc

    implicit none

    ! local variables
    integer (int_kind) :: n
    integer (int_kind) :: f

    real (real_kind), dimension(nvc,ncells) :: Coord

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Compute face centroid physical coordinates for each dimension
    NDIM_LOOP: do n = 1,ndim
       ! Gather the appropriate vertex coordinates
       call EN_GATHER (Coord, Vertex%Coord(n), BOUNDARY=Vrtx_Bdy(n)%Data)

       ! Convert logical coordinates to physical
       FACE_LOOP: do f = 1,nfc

          call LINEAR_PROP (f, Coord, cell(:)%face_centroid(n,f))

       end do FACE_LOOP

    end do NDIM_LOOP

    return

  END SUBROUTINE FACE_CENTROID_PHYSICAL

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE CELL_CENTROID ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the cell centroid, which is part of the Cell structure.
    !     Input  - Vertex%Coord(ndim)
    !     Output - Cell%Centroid(ndim)
    !
    !=======================================================================
    use ArrayAllocate_Module, only: ARRAYCREATE, ARRAYDESTROY
    use constants_module,     only: one_half, one_fourth, one_twelfth,  &
                                    twenty_four, zero, two, three, six, &
                                    four, twelve
    use gs_module,            only: EN_GATHER
    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Cell, Vertex, Vrtx_Bdy
    use parameter_module,     only: ncells, ndim, nvc, nrot

    implicit none

    ! Local Variables
    integer(KIND = int_kind) :: i, i1, i2, v1, v2, v3, v4, v5, v6, v7, v8
    real(KIND = real_kind), pointer, dimension(:,:) :: L, M, N, LxD3, MxD2, NxD1
    real(KIND = real_kind), pointer, dimension(:,:) :: Tmp, D1, D2, D3, Dv, &
                                                       D1xDv, D2xDv, D3xDv
    real(KIND = real_kind), pointer, dimension(:,:,:) :: Xv

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Explicitly allocate temporaries
    call ARRAYCREATE (L, 1, ndim, 1, ncells, 'Array L(ndim,ncells)')
    call ARRAYCREATE (M, 1, ndim, 1, ncells, 'Array M(ndim,ncells)')
    call ARRAYCREATE (N, 1, ndim, 1, ncells, 'Array N(ndim,ncells)')
    call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, ncells, 'Array Xv(ndim,nvc,ncells)')
    call ARRAYCREATE (LxD3, 1, nrot, 1, ncells, 'Array LxD3(nrot,ncells)')
    call ARRAYCREATE (MxD2, 1, nrot, 1, ncells, 'Array MxD2(nrot,ncells)')
    call ARRAYCREATE (NxD1, 1, nrot, 1, ncells, 'Array NxD1(nrot,ncells)')

    ! Only need these for 3-D.
    if (ndim == 3) then
       call ARRAYCREATE (Tmp, 1, ndim, 1, ncells, 'Array Tmp(ndim,ncells)')
       call ARRAYCREATE (D1, 1, ndim, 1, ncells, 'Array D1(ndim,ncells)')
       call ARRAYCREATE (D2, 1, ndim, 1, ncells, 'Array D2(ndim,ncells)')
       call ARRAYCREATE (D3, 1, ndim, 1, ncells, 'Array D3(ndim,ncells)')
       call ARRAYCREATE (Dv, 1, ndim, 1, ncells, 'Array Dv(ndim,ncells)')
       call ARRAYCREATE (D1xDv, 1, nrot, 1, ncells, 'Array D1xDv(nrot,ncells)')
       call ARRAYCREATE (D2xDv, 1, nrot, 1, ncells, 'Array D2xDv(nrot,ncells)')
       call ARRAYCREATE (D3xDv, 1, nrot, 1, ncells, 'Array D3xDv(nrot,ncells)')
    end if

    ! Set vertex numbers
    select case(ndim)
       case (2)
          v1 = 1; v2 = 2; v3 = 3; v4 = 4
       case (3)
          v1 = 1; v2 = 2; v3 = 3; v4 = 4
          v5 = 5; v6 = 6; v7 = 7; v8 = 8
    end select

    ! Gather vertex coordinates
    do i = 1,ndim
       call EN_GATHER (Xv(i,:,:), Vertex%Coord(i), BOUNDARY=Vrtx_Bdy(i)%Data)
    end do

    ! Compute quantities needed for the centroid coordinates
    do i = 1,ndim

       select case (ndim)

          case (2)

             L(i,:) = Xv(i,v1,:) - Xv(i,v4,:)
             M(i,:) = Xv(i,v3,:) - Xv(i,v4,:)
             N(i,:) = - Xv(i,v1,:) + Xv(i,v2,:) - Xv(i,v3,:) + Xv(i,v4,:)

          case (3)

             L(i,:) = Xv(i,v1,:) + Xv(i,v2,:) + Xv(i,v6,:) + Xv(i,v5,:) &
                    - Xv(i,v3,:) - Xv(i,v4,:) - Xv(i,v8,:) - Xv(i,v7,:)
             M(i,:) = Xv(i,v2,:) + Xv(i,v3,:) + Xv(i,v7,:) + Xv(i,v6,:) &
                    - Xv(i,v4,:) - Xv(i,v1,:) - Xv(i,v5,:) - Xv(i,v8,:)
             N(i,:) = Xv(i,v8,:) + Xv(i,v5,:) + Xv(i,v6,:) + Xv(i,v7,:) &
                    - Xv(i,v3,:) - Xv(i,v2,:) - Xv(i,v1,:) - Xv(i,v4,:)
             D1(i,:) = Xv(i,v2,:) + Xv(i,v6,:) + Xv(i,v4,:) + Xv(i,v8,:) &
                     - Xv(i,v1,:) - Xv(i,v5,:) - Xv(i,v3,:) - Xv(i,v7,:)
             D2(i,:) = Xv(i,v3,:) + Xv(i,v4,:) + Xv(i,v5,:) + Xv(i,v6,:) &
                     - Xv(i,v1,:) - Xv(i,v2,:) - Xv(i,v7,:) - Xv(i,v8,:)
             D3(i,:) = Xv(i,v1,:) + Xv(i,v4,:) + Xv(i,v6,:) + Xv(i,v7,:) &
                     - Xv(i,v2,:) - Xv(i,v3,:) - Xv(i,v5,:) - Xv(i,v8,:)
             Dv(i,:) = Xv(i,v1,:) + Xv(i,v3,:) + Xv(i,v6,:) + Xv(i,v8,:) &
                     - Xv(i,v2,:) - Xv(i,v4,:) - Xv(i,v5,:) - Xv(i,v7,:)

       end select

    end do

    select case (ndim)

       case (2)

          i1 = 1; i2 = 2 
          do i = 1,nrot
             LxD3(i,:) = L(i1,:)*M(i2,:) - L(i2,:)*M(i1,:)
             MxD2(i,:) = L(i1,:)*N(i2,:) - L(i2,:)*N(i1,:)
             NxD1(i,:) = N(i1,:)*M(i2,:) - N(i2,:)*M(i1,:)
          end do

       case (3)

          L = one_fourth*L
          M = one_fourth*M
          N = one_fourth*N
          Tmp = zero

          do i = 1,nrot

             select case (i)
             case (1)
                i1 = 2; i2 = 3
             case (2)
                i1 = 3; i2 = 1
             case (3)
                i1 = 1; i2 = 2
             end select

             LxD3(i,:) = L(i1,:)*D3(i2,:) - L(i2,:)*D3(i1,:)
             MxD2(i,:) = M(i1,:)*D2(i2,:) - M(i2,:)*D2(i1,:)
             NxD1(i,:) = N(i1,:)*D1(i2,:) - N(i2,:)*D1(i1,:)

             D1xDv(i,:) = D1(i1,:)*Dv(i2,:) - D1(i2,:)*Dv(i1,:)
             D2xDv(i,:) = D2(i1,:)*Dv(i2,:) - D2(i2,:)*Dv(i1,:)
             D3xDv(i,:) = D3(i1,:)*Dv(i2,:) - D3(i2,:)*Dv(i1,:)

          end do

          do i = 1,ndim

             Tmp(v1,:) = Tmp(v1,:) + L(i,:)*(MxD2(i,:) - NxD1(i,:)) &
                        + one_twelfth*(N(i,:)*D2xDv(i,:) - M(i,:)*D1xDv(i,:))
             Tmp(v2,:) = Tmp(v2,:) + M(i,:)*(NxD1(i,:) - LxD3(i,:)) &
                        + one_twelfth*(L(i,:)*D1xDv(i,:) - N(i,:)*D3xDv(i,:))
             Tmp(v3,:) = Tmp(v3,:) + N(i,:)*(LxD3(i,:) - MxD2(i,:)) &
                        + one_twelfth*(M(i,:)*D3xDv(i,:) - L(i,:)*D2xDv(i,:))

          end do

          do i = 1,ndim
             Tmp(i,:) = one_half + Tmp(i,:)/(twenty_four*Cell%Volume)
          end do

    end select

    ! Compute the centroid
    do i = 1,ndim

       select case (ndim)

          case (2)

             Cell%Centroid(i) = LxD3(1,:)*(twelve*Xv(i,v4,:) + six*(L(i,:) + M(i,:)) + three*N(i,:)) + &
                                MxD2(1,:)*(six*Xv(i,v4,:) + four*L(i,:) + three*M(i,:) + two*N(i,:)) + &
                                NxD1(1,:)*(six*Xv(i,v4,:) + three*L(i,:) + four*M(i,:) + two*N(i,:))
             Cell%Centroid(i) = Cell%Centroid(i)/(twelve*Cell%Volume)

          case (3)

             Cell%Centroid(i) = Xv(i,v4,:) &
                              + Tmp(v1,:)*(Xv(i,v1,:) - Xv(i,v4,:)) &
                              + Tmp(v2,:)*(Xv(i,v3,:) - Xv(i,v4,:)) &
                              + Tmp(v1,:)*Tmp(v2,:)*(Xv(i,v2,:) + Xv(i,v4,:) - Xv(i,v1,:) - Xv(i,v3,:)) &
                              + Tmp(v3,:)*(Xv(i,v8,:) - Xv(i,v4,:) &
                              + Tmp(v1,:)*(Xv(i,v4,:) + Xv(i,v5,:) - Xv(i,v1,:) - Xv(i,v8,:)) &
                              + Tmp(v2,:)*(Xv(i,v4,:) + Xv(i,v7,:) - Xv(i,v3,:) - Xv(i,v8,:)) &
                              + Tmp(v1,:)*Tmp(v2,:)*Dv(i,:))

       end select

    end do

    ! Explicitly deallocate temporaries
    call ARRAYDESTROY (L, 'Array L(ndim,ncells)')
    call ARRAYDESTROY (M, 'Array M(ndim,ncells)')
    call ARRAYDESTROY (N, 'Array N(ndim,ncells)')
    call ARRAYDESTROY (Xv, 'Array Xv(ndim,nvc,ncells)')
    call ARRAYDESTROY (LxD3, 'Array LxD3(nrot,ncells)')
    call ARRAYDESTROY (MxD2, 'Array MxD2(nrot,ncells)')
    call ARRAYDESTROY (NxD1, 'Array NxD1(nrot,ncells)')
    if (ndim == 3) then
       call ARRAYDESTROY (Tmp, 'Array Tmp(nvc,ncells)')
       call ARRAYDESTROY (D1, 'Array D1(ndim,ncells)')
       call ARRAYDESTROY (D2, 'Array D2(ndim,ncells)')
       call ARRAYDESTROY (D3, 'Array D3(ndim,ncells)')
       call ARRAYDESTROY (Dv, 'Array Dv(ndim,ncells)')
       call ARRAYDESTROY (D1xDv, 'Array D1xDv(nrot,ncells)')
       call ARRAYDESTROY (D2xDv, 'Array D2xDv(nrot,ncells)')
       call ARRAYDESTROY (D3xDv, 'Array D3xDv(nrot,ncells)')
    end if

    return

  END SUBROUTINE CELL_CENTROID

  SUBROUTINE CELL_VOLUME ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the cell volume, which is part of the Cell structure.
    !   The volume is computed according to the prescription given
    !   by J. Dukowicz, JCP 74: 493-496 (1988).
    !     Input  - Vertex%Coord(ndim)
    !     Output - Cell%Volume
    !
    !=======================================================================
    use constants_module,     only: one, one_half, zero
    use cutoffs_module,       only: alittle
    use gs_module,            only: EN_GATHER, EN_SUM_SCATTER
    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Cell, Vertex, Vrtx_Bdy, Mesh,       &
                                    volume_min, volume_max, GAP_ELEMENT_1
    use parameter_module,     only: ncells, ndim, nfc, nnodes, nvc, nec
    use pgslib_module,        only: PGSLib_GLOBAL_MINLOC, &
                                    PGSLib_GLOBAL_MAXLOC, &
                                    PGSLib_Global_SUM,    &
                                    PGSLib_GLOBAL_MINVAL, &
                                    PGSLib_GLOBAL_MAXVAL, &
                                    PGSLib_GLOBAL_ANY
    use truchas_logging_services

    implicit none

    ! Local Variables
    integer(KIND = int_kind) :: f, n, v1, v2, v3, v4, v5, v6, icell
    integer(KIND = int_kind), dimension(1)               :: MinLoc_L, MaxLoc_L
    real(KIND = real_kind),   dimension(nnodes)          :: Coord
    real(KIND = real_kind),   dimension(ndim,nvc,ncells) :: Xn
    real(KIND = real_kind),   dimension(ndim,ncells)     :: X1, X2, X3
    character(128) :: message
    real(KIND = real_kind) :: total_volume

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Gather vertex coordinates.
    do n = 1,ndim
       call EN_GATHER (Xn(n,:,:), Vertex%Coord(n), BOUNDARY=Vrtx_Bdy(n)%Data)
    end do

    ! Loop over faces, accumulating the cell volume.
    FACE_LOOP: do f = 1,nfc

       select case (f)

       case (1)

          ! Side 1
          v1 = 8; v2 = 4; v3 = 7; v4 = 8; v5 = 3; v6 = 4

       case (2)

          ! Side 2
          v1 = 6; v2 = 2; v3 = 5; v4 = 6; v5 = 1; v6 = 2

       case (3)

          ! Side 3
          v1 = 5; v2 = 1; v3 = 8; v4 = 5; v5 = 4; v6 = 1

       case (4)

          ! Side 4
          v1 = 7; v2 = 3; v3 = 6; v4 = 7; v5 = 2; v6 = 3

       case (5)

          ! Side 5
          v1 = 3; v2 = 4; v3 = 2; v4 = 3; v5 = 1; v6 = 4

       case (6)

          ! Side 6
          v1 = 6; v2 = 5; v3 = 7; v4 = 6; v5 = 8; v6 = 5

       end select

       select case (ndim)

       case (2)
             ! 2-D
          X1 = one_half
          v1 = v2
          v3 = v5
          v4 = v6
          do n = 1,ndim
             X2(n,:) = Xn(n,v3,:) + Xn(n,v5,:)
             X3(n,:) = Xn(n,v2,:) + Xn(n,v4,:)
          end do
          v1 = 1; v2 = 1; v3 = 2
          Cell%Volume = Cell%Volume + X1(v1,:)*(X2(v2,:)*X3(v3,:) - X3(v2,:)*X2(v3,:))

       case (3)
          ! 3-D
          do n = 1,ndim
             X1(n,:) = Xn(n,v1,:) + Xn(n,v2,:)
             X2(n,:) = Xn(n,v3,:) + Xn(n,v4,:)
             X3(n,:) = Xn(n,v5,:) + Xn(n,v6,:)
          end do
          do n = 1,ndim
             select case (n)
             case (1)
                v1 = 1; v2 = 2; v3 = 3
             case (2)
                v1 = 2; v2 = 3; v3 = 1
             case (3)
                v1 = 3; v2 = 1; v3 = 2
             end select
             Cell%Volume = Cell%Volume + X1(v1,:)*(X2(v2,:)*X3(v3,:) - X3(v2,:)*X2(v3,:))
          end do

       end select

    end do FACE_LOOP

    Cell%Volume = Cell%Volume/REAL(nec)

    ! Look for gap elements.  Maybe make this conditional based on some input flag?
    GAP_ELE_CHECK: do icell = 1,ncells
       if (abs(Cell(icell)%Volume) < (10.0*alittle)) then
          if (Mesh(icell)%Cell_Shape >= GAP_ELEMENT_1) then
             Cell(icell)%Volume = 2.0*alittle
          else
             call TLS_panic ('CELL_VOLUME: mesh contains cells with very small volume that do not appear to be gap elements')
          end if
       end if
    end do GAP_ELE_CHECK

    ! Make sure volumes are OK.
    if (PGSLib_Global_ANY(Cell%Volume < zero)) then
       call TLS_panic ('CELL_VOLUME: mesh contains cells with negative volumes')
    end if

    ! Write out cell volume extrema and total.
    MinLoc_L     = PGSLib_Global_MINLOC(Cell%Volume)
    MaxLoc_L     = PGSLib_Global_MAXLOC(Cell%Volume)
    total_volume = PGSLib_Global_SUM(Cell%Volume)
    volume_min    = PGSLib_Global_MINVAL(Cell%Volume)
    volume_max    = PGSLib_Global_MAXVAL(Cell%Volume)
    
    call TLS_info ('')
    call TLS_info ('                    Min        Cell       Max        Cell       Total')
    call TLS_info ('                    ---        ----       ---        ----       -----')
    write (message,5) 'Volumes', volume_min, MinLoc_L, volume_max, MaxLoc_L, total_volume  
5   format (7x,a,2x,1pe11.4,2x,i8,2x,1pe11.4,2x,i8,3x,1pe11.4)
    call TLS_info (message)

    ! Scatter and store the reciprocal sum of the inverse volumes.
    X1(1,:) = one/Cell%Volume
    
    call EN_SUM_SCATTER (Coord, X1(1,:))
    Vertex%Rsum_rvol = one/Coord

    return

  END SUBROUTINE CELL_VOLUME

  SUBROUTINE CELL_WIDTH ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the distance from each cell face centroid to its
    !   own cell centroid.
    !
    !=======================================================================
    use kind_module,          only: int_kind
    use mesh_module,          only: Cell
    use parameter_module,     only: ndim, nfc

    implicit none

    ! Local Variables
    integer(KIND = int_kind) :: f, n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Loop over faces, calculating the distance from the center of the cell to
    ! the center of each face.
    do f = 1,nfc
       ! Loop over dimensions, calculating the sum of the squares of the coordinates.
       do n = 1,ndim
          Cell(:)%Halfwidth(f) = Cell(:)%Halfwidth(f) + (Cell(:)%Centroid(n) - cell(:)%face_centroid(n,f))**2
       end do
       ! Take the square root of the sum of squares to make a distance to each face.
       Cell%Halfwidth(f) = SQRT(Cell%Halfwidth(f))
    end do

    return

  END SUBROUTINE CELL_WIDTH

  SUBROUTINE CELL_TYPE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Evaluate Mesh%Cell_Shape based on the number of degenerate faces.
    !   Prisms & pyramids differentiated by the number
    !
    !=======================================================================
    use kind_module,          only: int_kind
    use mesh_module,          only: Mesh, DEGENERATE_FACE, &
                                    CELL_HEX, CELL_PYRAMID, CELL_PRISM, &
                                    CELL_TET, GAP_ELEMENT_1, GAP_ELEMENT_3, GAP_ELEMENT_5
    use parameter_module,     only: ncells

    implicit none

    ! Local Variables
    integer(KIND = int_kind) :: i

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do i = 1,ncells
       ! gap elements are already identified
       if (Mesh(i)%Cell_Shape == GAP_ELEMENT_1) cycle
       if (Mesh(i)%Cell_Shape == GAP_ELEMENT_3) cycle
       if (Mesh(i)%Cell_Shape == GAP_ELEMENT_5) cycle
       ! hex cell if face 6 isn't degenerate
       if (Mesh(i)%Ngbr_Face(6) /= DEGENERATE_FACE) then
          Mesh(i)%Cell_Shape = CELL_HEX

       else
          ! tet cell if face 2 is degenerate
          if (Mesh(i)%Ngbr_Face(2) == DEGENERATE_FACE) then
             Mesh(i)%Cell_Shape = CELL_TET

          else
             ! prism cell if vertex 5 isn't the same as vertex 6;
             ! ortherwise a prism
             if (Mesh(i)%Ngbr_Vrtx(5) /= Mesh(i)%Ngbr_Vrtx(6)) then
                Mesh(i)%Cell_Shape = CELL_PRISM

             else
                Mesh(i)%Cell_Shape = CELL_PYRAMID

             end if
          end if
       end if
    end do

    return

  END SUBROUTINE CELL_TYPE

  SUBROUTINE FACE_AREA ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face areas (Face_Area) and unit normals (Face_Normal),
    !   which are part of the Cell derived type.  Face unit normal vectors
    !   are assigned for each face of the cell, and are assumed to be 
    !   pointing out of each face into the neighboring cell.
    !
    !=======================================================================
    use constants_module, only: one_half, zero
    use cutoffs_module,   only: alittle
    use gs_module,        only: EN_GATHER
    use kind_module,      only: int_kind, real_kind
    use mesh_module,      only: Cell, Vertex, Vrtx_Bdy
    use parameter_module, only: ncells, ndim, nfc, nvc

    implicit none

    ! Local Variables
    integer(KIND = int_kind) :: f, i, v1, v2, v3, v4
    real(KIND = real_kind), dimension(ndim,nvc,ncells) :: Xn
    real(KIND = real_kind), dimension(ndim,ncells)     :: X1, X2

#if DEBUG_FACE_VECTORS
    integer(KIND = int_kind) :: ff
    real(KIND = real_kind), dimension(nfc,ncells) :: Face_coord
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Gather vertex coordinates.
    do i = 1,ndim
       call EN_GATHER (Xn(i,:,:), Vertex%Coord(i), BOUNDARY=Vrtx_Bdy(i)%Data)
    end do

    ! Loop over faces, computing face unit normals and areas.
    FACE_AREA_LOOP: do f = 1,nfc

       select case (f)

          case (1)

          ! Side 1 (vertices 4-8-7-3)
          v1 = 8; v2 = 3; v3 = 7; v4 = 4

          case (2)

          ! Side 2 (vertices 1-2-6-5)
          v1 = 6; v2 = 1; v3 = 5; v4 = 2

          case (3)

          ! Side 3 (vertices 4-1-5-8)
          v1 = 5; v2 = 4; v3 = 8; v4 = 1

          case (4)

          ! Side 4 (vertices 3-7-6-2)
          v1 = 7; v2 = 2; v3 = 6; v4 = 3

          case (5)

          ! Side 5 (vertices 4-3-2-1)
          v1 = 3; v2 = 1; v3 = 2; v4 = 4

          case (6)

          ! Side 6 (vertices 8-7-6-5)
          v1 = 6; v2 = 8; v3 = 7; v4 = 5

       end select

       if (ndim == 2) then
          v1 = v4
          v3 = v2 
       end if

       do i = 1,ndim
          X1(i,:) = Xn(i,v1,:) - Xn(i,v2,:)
          X2(i,:) = Xn(i,v3,:) - Xn(i,v4,:)
       end do

       ! These are face area vectors for now.
       do i = 1,ndim
          select case (i)
             case (1)
                v1 = 2; v2 = 3
                if (ndim == 2) Cell%Face_Normal(i,f) = X1(v1,:)
             case (2)
                v1 = 3; v2 = 1
                if (ndim == 2) Cell%Face_Normal(i,f) = X2(v2,:)
             case (3)
                v1 = 1; v2 = 2
          end select
          if (ndim == 3) Cell%Face_Normal(i,f) = one_half*(X1(v1,:)*X2(v2,:) - X2(v1,:)*X1(v2,:))
       end do

       ! Set components to zero if they're small; accumulate areas.
       do i = 1,ndim
          Cell%Face_Normal(i,f) = MERGE(zero, Cell%Face_Normal(i,f), &
                                        ABS(Cell%Face_Normal(i,f)) < alittle)
          Cell%Face_Area(f) = Cell%Face_Area(f) + Cell%Face_Normal(i,f)**2
       end do

       ! Compute areas.
       Cell%Face_Area(f) = SQRT(Cell%Face_Area(f))

       ! Convert the face vectors to unit normals.
       do i = 1,ndim
          where (Cell%Face_Area(f) >= alittle)
             Cell%Face_Normal(i,f) = Cell%Face_Normal(i,f) / Cell%Face_Area(f)
          elsewhere
             Cell%Face_Normal(i,f) = zero
             Cell%Face_Area(f) = zero
          end where
       end do

    end do FACE_AREA_LOOP

#if DEBUG_FACE_VECTORS
    ! Check to see if adjacent faces have equal and opposite area vectors.
    CHECK_FACE_VECTORS: do f = 1,nfc

       ! Gather face area vector for the neighbor of this face f
       X1 = zero; X2 = zero
       do i = 1,ndim
          call EE_GATHER (DEST = Face_coord, SRC = Cell%Face_Normal(i,f))
          do ff = 1,nfc
             where (Mesh%Ngbr_face(ff) == f) 
                X1(i,:) = Face_coord(ff,:)
                X2(i,:) = Cell%Face_Normal(i,ff)
             end where
          end do
       end do

       ! Take the cross product, which should be zero
       do i = 1,ndim
          select case (i)
             case (1)
                v1 = 2; v2 = 3
                if (ndim == 2) v2 = 1
             case (2)
                v1 = 3; v2 = 1
                if (ndim == 2) v1 = 2
             case (3)
                v1 = 1; v2 = 2
          end select
          Xn(i,1,:) = X1(v1,:)*X2(v2,:) - X2(v1,:)*X1(v2,:)
          Xn(i,1,:) = MERGE(zero, Xn(i,1,:), ABS(Xn(i,1,:)) < ten_tominus14)
       end do

       if (PGSLib_Global_ANY(Xn(:,1,:) /= zero)) then
          call TLS_panic ('FACE_AREA: Some adjacent face area vectors are not equal and opposite')
       end if

    end do CHECK_FACE_VECTORS
#endif

    return

  END SUBROUTINE FACE_AREA

  SUBROUTINE FACE_CENTROID_LOGICAL ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face centroid (Face_centroid), which is part
    !   of the Cell structure. A face centroid is assigned for
    !   each face of the cell, and is a LOGICAL coordinate.
    !     Input  - Vertex%Coord(ndim)
    !     Output - Cell%Face_centroid_L(ndim,nfc)
    !
    !=======================================================================
    use ArrayAllocate_Module, only: ARRAYCREATE, ARRAYDESTROY
    use constants_module,     only: one, one_half, one_twelfth, zero
    use cutoffs_module,       only: alittle
    use gs_module,            only: EN_GATHER
    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Cell, Vertex, Vrtx_Bdy
    use parameter_module,     only: ncells, ndim, nfc, nvc

    implicit none

    ! Local Variables
    integer(KIND = int_kind) :: f, i, i1, i2, n, v11, v12, v13, v14,    &
                                v21, v22, v23, v24, v31, v32, v33, v34, &
                                nn
    integer(KIND = int_kind) :: n1 = 1, n2 = 2, n3 = 3
    real(KIND = real_kind), dimension(ndim)           :: Coef
    real(KIND = real_kind), dimension(ncells)         :: Face_Area
    real(KIND = real_kind), pointer, dimension(:,:,:) :: Xv

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Only do this if it's 3-D.
    if (ndim == 3) then
       ! Allocate temporaries
       call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, ncells, 'Array Xv(ndim,nvc,ncells)')

       ! Gather vertex coordinates
       do i = 1,ndim
          call EN_GATHER (Xv(i,:,:), Vertex%Coord(i), BOUNDARY=Vrtx_Bdy(i)%Data)
       end do
    end if

    ! Loop over faces
    FACE_LOOP: do f = 1,nfc

       ! Preinitialize the face centroid
       do n = 1,ndim
          Cell%Face_Centroid_L(n,f) = one_half
          Coef(n)                   = one
       end do

       ! We're done if this is 2-D.
       if (ndim == 2) cycle FACE_LOOP

       Face_Area = one_twelfth / (Cell%Face_Area(f) + alittle)

       select case (f)

          case (1)

          ! Side 1 (vertices 4-8-7-3)
          v11 = 1; v12 = 1; v13 = 1; v14 = 1
          v21 = 7; v22 = 8; v23 = 3; v24 = 4
          v31 = 8; v32 = 4; v33 = 7; v34 = 3
          nn  = 1

          Cell%Face_centroid_L(nn,f) = zero; Coef(nn) = zero

          case (2)

          ! Side 2 (vertices 1-2-6-5)
          v11 = 1; v12 = 1; v13 = 1; v14 = 1
          v21 = 2; v22 = 1; v23 = 6; v24 = 5
          v31 = 6; v32 = 2; v33 = 5; v34 = 1
          nn  = 1

          Cell%Face_centroid_L(nn,f) = one; Coef(nn) = zero

          case (3)

          ! Side 3 (vertices 4-1-5-8)
          v11 = 1; v12 = 4; v13 = 5; v14 = 8
          v21 = 1; v22 = 1; v23 = 1; v24 = 1
          v31 = 5; v32 = 1; v33 = 8; v34 = 4
          nn  = 2

          Cell%Face_centroid_L(nn,f) = zero; Coef(nn) = zero

          case (4)

          ! Side 4 (vertices 3-7-6-2)
          v11 = 6; v12 = 7; v13 = 2; v14 = 3
          v21 = 1; v22 = 1; v23 = 1; v24 = 1
          v31 = 7; v32 = 3; v33 = 6; v34 = 2
          nn  = 2

          Cell%Face_centroid_L(nn,f) = one; Coef(nn) = zero

          case (5)

          ! Side 5 (vertices 4-3-2-1)
          v11 = 2; v12 = 3; v13 = 1; v14 = 4
          v21 = 3; v22 = 4; v23 = 2; v24 = 1
          v31 = 1; v32 = 1; v33 = 1; v34 = 1
          nn  = 3

          Cell%Face_centroid_L(nn,f) = zero; Coef(nn) = zero

          case (6)

          ! Side 6 (vertices 8-5-6-7)
          v11 = 5; v12 = 8; v13 = 6; v14 = 7
          v21 = 6; v22 = 5; v23 = 7; v24 = 8
          v31 = 1; v32 = 1; v33 = 1; v34 = 1
          nn  = 3

          Cell%Face_centroid_L(nn,f) = one; Coef(nn) = zero

       end select

       NDIM_LOOP: do i = 1,ndim

          select case (i)
             case (1)
                i1 = 2; i2 = 3
             case (2)
                i1 = 3; i2 = 1
             case (3)
                i1 = 1; i2 = 2
          end select

          Cell%Face_centroid_L(n1,f) = Cell%Face_centroid_L(n1,f) + &
                Coef(n1)*Face_Area*Cell%Face_Normal(i,f)* &
                ((Xv(i1,v11,:) - Xv(i1,v12,:))*(Xv(i2,v13,:) - Xv(i2,v14,:)) &
               - (Xv(i2,v11,:) - Xv(i2,v12,:))*(Xv(i1,v13,:) - Xv(i1,v14,:)))

          Cell%Face_centroid_L(n2,f) = Cell%Face_centroid_L(n2,f) + &
               Coef(n2)*Face_Area*Cell%Face_Normal(i,f)* &
               ((Xv(i1,v21,:) - Xv(i1,v22,:))*(Xv(i2,v23,:) - Xv(i2,v24,:)) &
              - (Xv(i2,v21,:) - Xv(i2,v22,:))*(Xv(i1,v23,:) - Xv(i1,v24,:)))

          Cell%Face_centroid_L(n3,f) = Cell%Face_centroid_L(n3,f) + &
               Coef(n3)*Face_Area*Cell%Face_Normal(i,f)* &
               ((Xv(i1,v31,:) - Xv(i1,v32,:))*(Xv(i2,v33,:) - Xv(i2,v34,:)) &
              - (Xv(i2,v31,:) - Xv(i2,v32,:))*(Xv(i1,v33,:) - Xv(i1,v34,:)))

       end do NDIM_LOOP

    end do FACE_LOOP

    ! Deallocate temporaries
    if (ASSOCIATED(Xv)) call ARRAYDESTROY (Xv, 'Array Xv(ndim,nvc,ncells)')

    return

  END SUBROUTINE FACE_CENTROID_LOGICAL

  SUBROUTINE JACOBIAN ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute and analyze the inverse transpose Jacobian for each cell.
    !   The Jacobian is computed as follows:
    !             -                        -
    !             |  dx/di   dy/di   dz/di |
    !       J =   |  dx/dj   dy/dj   dz/dj |
    !             |  dx/dk   dy/dk   dz/dk |
    !             -                        -
    !   where di = face 2 - face 1, dj = face 4 - face 3, and
    !   dk = face 6 - face 5.
    !
    !=======================================================================
    use ArrayAllocate_Module, only: ARRAYCREATE, ARRAYDESTROY
    use constants_module,     only: zero
    use cutoffs_module,       only: alittle
    use discrete_ops_data,    only: use_ortho_face_gradient, discrete_ops_type
    use kind_module,          only: int_kind, log_kind, real_kind
    use mesh_module,          only: Cell, orthogonal_mesh
    use parameter_module,     only: ncells, ncells_tot, ndim
    use pgslib_module,        only: PGSLib_GLOBAL_COUNT
    use truchas_logging_services
    implicit none

    ! Local Variables
    integer (int_kind) :: l
    integer (int_kind) :: m
    integer (int_kind) :: n
    integer (int_kind) :: orthogonal_cells
    integer (int_kind) :: r
    logical (log_kind), pointer, dimension(:)   :: Mask
    real (real_kind), pointer, dimension(:,:)   :: Xn
    real (real_kind), pointer, dimension(:,:,:) :: jacob
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! allocate temporaries
    call ARRAYCREATE (mask,  1, ncells, 'jacobian: mask(ncells)')
    call ARRAYCREATE (Xn,    1, ndim, 1, ncells, 'jacobian: Xn(ndim,ncells)')
    call ARRAYCREATE (Jacob, 1, ndim, 1, ndim, 1, ncells, 'jacobian: Jacob(ndim,ndim,ncells)')

    ! Loop over the number of dimensions, and compute an approximate
    ! Jacobian by taking differences of face coordinates
    do n = 1,ndim
       do m = 1,ndim
          r = m*2
          l = r-1
          Jacob(m,n,:) = cell(:)%face_centroid(n,l) - cell(:)%face_centroid(n,r)
          Jacob(m,n,:) = MERGE(zero, Jacob(m,n,:), ABS(Jacob(m,n,:)) <= alittle)
       end do
    end do

    ! Compute the dot product of each Jacobian vector
    ! with each of other Jacobian vectors: 1*2 in 2D; 1*2, 2*3, 1*3 in 3D
    Xn = zero
    select case (ndim)
       case (1)
          ! a 1D mesh is inherently orthogonal
       case (2)
          Xn(1,:) = Jacob(1,1,:)*Jacob(2,1,:) + Jacob(1,2,:)*Jacob(2,2,:)
       case (3)
          do l = 1,3
             Xn(1,:) = Xn(1,:) + Jacob(1,l,:)*Jacob(2,l,:)
             Xn(2,:) = Xn(2,:) + Jacob(2,l,:)*Jacob(3,l,:)
             Xn(3,:) = Xn(3,:) + Jacob(1,l,:)*Jacob(3,l,:)
          end do
    end select
    do n = 1, ndim
       Xn(n,:) = MERGE(zero, Xn(n,:), ABS(Xn(n,:)) <= alittle)
    end do

    ! Now count the number of cells that have
    ! all 3 Jacobian vectors mutually orthogonal.
    Mask = .true.
    do n = 1,ndim
       Mask = Mask .and. Xn(n,:) == zero
    end do
    orthogonal_cells = PGSLib_Global_COUNT(Mask)

    ! Print out diagnostic mesh information and 
    ! set orthogonal mesh mask if conditions are right.
    if (orthogonal_cells == ncells_tot) then
       orthogonal_mesh = .true.
       if (discrete_ops_type == 'default') use_ortho_face_gradient = .true.
       call TLS_info ('')
       call TLS_info ('    Entire mesh is orthogonal')
       call TLS_info ('')
    else
       call TLS_info ('')
       orthogonal_mesh = .false.
       if (discrete_ops_type == 'default') use_ortho_face_gradient = .false.
       write (message,20) orthogonal_cells, ncells_tot
20     format (9x,i8,' out of ',i8,' total cells are orthogonal')
       call TLS_info ('')
       call TLS_info (message)
       call TLS_info ('')
    end if

    write (message,'(4x,a,l1)') 'use_ortho_face_gradient = ',use_ortho_face_gradient
    call TLS_info (message)
    call TLS_info ('')
    write (message,'(4x,a)') 'Using full pivoting for LSLR_ operators'
    call TLS_info (message)
    call TLS_info ('')

    ! deallocate temporaries
    call ARRAYDESTROY (jacob, 'jacobian: jacob(ndim,ndim,ncells)')
    call ARRAYDESTROY (xn,    'jacobian: xn(nvc,ncells)')
    call ARRAYDESTROY (mask,  'jacobian: mask(ncells)')

    return

  END SUBROUTINE JACOBIAN

  SUBROUTINE GAP_ELEMENT_TYPE ()
    ! Purpose:
    !  This routine used the body input and mesh data to assign gap element 
    !  types to the cell_shape attribute.  It is called before the rest of the
    !  cell geometry calls so that cell_volume can use this info to ignore zero
    !  volume cells that are gap elements.
    !
    !  Author(s):
    !  Dave Korzekwa (dak@lanl.gov)
    !___________________________________________________________________________


    use mesh_module,       only: Mesh, GAP_ELEMENT_1, GAP_ELEMENT_3, GAP_ELEMENT_5
    use parameter_module,  only: ncells, ndim, nfc
    use kind_module,       only: int_kind, log_kind
    use gs_module,         only: EE_GATHER
    use mesh_input_module, only: gap_element_blocks

    ! Local variables
    integer                                          :: status
    integer(kind = int_kind)                         :: icell, idim, gface, nbid, gap_type
    integer(kind = int_kind),pointer, dimension(:,:) :: NbrBlkID
    logical(log_kind)                                :: found_gap, nbgap
    

    if (ANY(gap_element_blocks > 0)) then
       allocate(NbrBlkID(nfc,ncells), STAT=status)
       CALL TLS_fatal_if_any ((status /= 0), 'GAP_ELEMENT_TYPE: error allocating array NbrBlkID(nfc,ncells)')
       call EE_GATHER (NbrBlkID, Mesh%CBlockID)

       do icell = 1,ncells
          ! Is this cell in this gap element body?
          if (ANY(gap_element_blocks == Mesh(icell)%CBlockID)) then
             found_gap = .false.
             do idim = 1,ndim
                ! Figure out which flavor of gap element we have
                gface = 2 * idim - 1
                nbid = NbrBlkID(gface,icell)
                nbgap = .false.
                if ((.not. ANY(gap_element_blocks == nbid)) .and. (nbid /= 0)) then
                   if (.not. found_gap) then
                      gap_type = idim
                      found_gap = .true.
                   else
                      CALL TLS_panic ('GAP_ELEMENT_TYPE: Error finding gap element type')
                   end if
                end if
             end do
             if (found_gap) then
                select case(gap_type)
                case(1)
                   Mesh(icell)%Cell_Shape = GAP_ELEMENT_1
                case(2)
                   Mesh(icell)%Cell_Shape = GAP_ELEMENT_3
                case(3)
                   Mesh(icell)%Cell_Shape = GAP_ELEMENT_5
                end select
             else
                CALL TLS_panic ('GAP_ELEMENT_TYPE: Gap element type not found')
             end if
          end if
       end do

       deallocate(NbrBlkID)
    end if
          
  END SUBROUTINE GAP_ELEMENT_TYPE

END MODULE CELL_GEOMETRY_MODULE
