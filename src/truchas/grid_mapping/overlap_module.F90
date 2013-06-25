!     Author:  S.J. Mosso

!     Last change:  SJM  March 27 2005
!     Corrected problem where grazing plane-polyhedra intersections
!     in PLANE_POLY_INT3D resulted in '... facial edge not found...'

MODULE overlap_module

  IMPLICIT NONE

  SAVE

   ! number of spatial dimensions
   INTEGER, PARAMETER :: ndim = 3

   ! number of vector dimensions
   INTEGER, PARAMETER :: mdim = 3

!   INTEGER, PARAMETER :: RK = 8
   INTEGER, PARAMETER :: RK = KIND(1.0d0)

   REAL(RK), PARAMETER :: S0   =   0.0_RK
   REAL(RK), PARAMETER :: S1   =   1.0_RK
   REAL(RK), PARAMETER :: S2   =   2.0_RK
   REAL(RK), PARAMETER :: S6   =   6.0_RK
   REAL(RK), PARAMETER :: S100 = 100.0_Rk
   REAL(RK), PARAMETER :: F1o2 = S1 / S2
   REAL(RK), PARAMETER :: F1o6 = S1 / S6
   
   ! Define a 'special' number of undefined
   REAL(RK), PARAMETER :: undef = 123451234._RK

   ! Define a small number
   REAL(RK), PARAMETER :: small = 1.e-10_RK

   INTEGER, PARAMETER :: maxFaces = 18
   INTEGER, PARAMETER :: maxCrnrs = 50
   INTEGER, PARAMETER :: maxCperF =  8
   INTEGER, PARAMETER :: maxEperF = 10
   INTEGER, PARAMETER :: maxEperV =  8
   INTEGER, PARAMETER :: maxEdges = 50

  TYPE vector_type
     REAL(RK), DIMENSION(mdim) :: x
  END TYPE vector_type

  TYPE(vector_type), PARAMETER :: vector_undef = &
                     vector_type( (/ UNDEF,UNDEF,UNDEF /) )
  TYPE(vector_type), PARAMETER :: Origin = &
                     vector_type ( (/ 0.0_RK, 0.0_RK, 0.0_RK /) )

  TYPE PLANE_TYPE
     TYPE(VECTOR_TYPE) :: n
     REAL(RK)          :: p
  END TYPE PLANE_TYPE

  TYPE (PLANE_TYPE), PARAMETER :: Plane_Undef = &
          PLANE_TYPE ( VECTOR_UNDEF, HUGE( S1 ) )

  TYPE poly3D_type
     INTEGER                                         :: numCrnr
     TYPE(VECTOR_TYPE), DIMENSION(maxCrnrs)          :: Crnr
     INTEGER                                         :: numFaces
     INTEGER                                         :: numEdges
     INTEGER,           DIMENSION(         maxFaces) :: nEdgeofFace
     INTEGER,           DIMENSION(maxEperF,maxFaces) :: EdgeofFace
     INTEGER,           DIMENSION(   2    ,maxEdges) :: CrnrofEdge
     INTEGER,           DIMENSION(         maxCrnrs) :: nEdgeofCrnr
     INTEGER,           DIMENSION(maxEperV,maxCrnrs) :: EdgeofCrnr
  END TYPE poly3D_type

  INTEGER :: ic_undef

  TYPE(POLY3D_TYPE), PARAMETER :: POLY3D_Init = POLY3D_TYPE( 0,  &
                       (/ (VECTOR_UNDEF,ic_undef=1,maxCrnrs) /), &
                                                              0, &
                                                              0, &
                       (/ (0,           ic_undef=1,maxFaces) /), &
              RESHAPE( (/ (0,  ic_undef=1,maxEperF*maxFaces) /), &
                       (/ maxEperF,maxFaces /)                ), &
              RESHAPE( (/ (0,  ic_undef=1,    2   *maxEdges) /), &
                       (/    2    ,maxEdges /)                ), &
                       (/ (0,           ic_undef=1,maxCrnrs) /), &
              RESHAPE( (/ (0,  ic_undef=1,maxEperV*maxCrnrs) /), &
                       (/ maxEperV,maxCrnrs /)                )  )

  INTERFACE Plane_Poly_Int
     MODULE PROCEDURE Plane_Poly_Int3D
  END INTERFACE

  INTERFACE Poly_Poly_Int
     MODULE PROCEDURE Poly_Poly_Int3D
  END INTERFACE

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE Poly3D_2_Poly3D
  END INTERFACE

  INTERFACE Volume_Poly
     MODULE PROCEDURE volm_Poly3D
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE negate_vector
     MODULE PROCEDURE negate_vector_array
     MODULE PROCEDURE subtract_vector
  END INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_vector
  END INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE scalar_mult_vector
     MODULE PROCEDURE vector_mult_scalar
     MODULE PROCEDURE vector_dot_vector
  END INTERFACE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE vector_div_scalar
  END INTERFACE

  INTERFACE ABS
     MODULE PROCEDURE abs_of_vector
  END INTERFACE

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE Vector_2_Vector
  END INTERFACE

CONTAINS

   SUBROUTINE Vol_Poly_Poly_Int (Tet, Hex, Vol)

   IMPLICIT NONE

   TYPE(POLY3D_TYPE), INTENT(IN)  :: Tet
   TYPE(POLY3D_TYPE), INTENT(IN)  :: Hex
   REAL(RK),          INTENT(OUT) :: Vol

   !------------------------------------------------------------------

   TYPE(POLY3D_TYPE) :: Int

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

   call Poly_Poly_Int3D (Tet, Hex, Int)
   Vol = Volm_Poly3D (Int)

   END SUBROUTINE Vol_Poly_Poly_Int

   SUBROUTINE Poly_Poly_Int3D (Poly1, Poly2, Int)

   IMPLICIT NONE

   TYPE(POLY3D_TYPE), INTENT(IN)  :: Poly1
   TYPE(POLY3D_TYPE), INTENT(IN)  :: Poly2
   TYPE(POLY3D_TYPE), INTENT(OUT) :: Int

   !------------------------------------------------------------------

   INTEGER :: f
   TYPE(plane_type) :: facet1
   INTEGER :: ie1, ie2, iP1, iP2, iP3
   TYPE(vector_type) :: P1, P2, P3
   TYPE(POLY3D_TYPE) :: Tmp

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

   Tmp = Poly2
   do f = 1, Poly1%numFaces

     ie1 = Poly1%EdgeofFace(1,f)
     ie2 = Poly1%EdgeofFace(2,f)
     call Get_Edge_Crnrs (Poly1, ie1, iP1, iP2)
     call Get_Edge_Crnrs (Poly1, ie2, iP2, iP3)

     P1 = Poly1%Crnr(iP1)
     P2 = Poly1%Crnr(iP2)
     P3 = Poly1%Crnr(iP3)

     facet1 = Plane_from_Points (P1, P2, P3)

     call Plane_Poly_Int3D (facet1, Tmp, Int)

     if (Int%numCrnr < 4) then
       EXIT
     else
       Tmp = Int
     end if
   end do

   END SUBROUTINE poly_poly_int3D


    SUBROUTINE PLANE_POLY_INT3D (plane, Poly, Int)

    !==================================================================
    !  This subroutine calculates the vertices of a polygon, Poly, that
    !  define the portion of the polygon that lies behind plane.
    !==================================================================

    IMPLICIT NONE

    TYPE(plane_type),  INTENT(IN)    :: plane
    TYPE(poly3D_type), INTENT(IN)    :: Poly
    TYPE(poly3D_type), INTENT(INOUT) :: Int

    !------------------------------------------------------------------

    REAL(RK), DIMENSION(maxCrnrs) :: CrnrDist
    REAL(RK) :: maxDist, minDist
    INTEGER :: c, f, e

    INTEGER :: n, nCrnr

    INTEGER                      :: iP1, iP2, ie
    TYPE(vector_type)            :: P1, P2, Vint
    REAL(RK)                     :: dP1, dP2

    INTEGER, DIMENSION(maxCrnrs) :: newCrnr
    INTEGER, DIMENSION(maxEdges) :: newEdge, facet_newEdge
    INTEGER :: newC, newE, last_newE, nOrigC
    INTEGER :: nEdgeBehind

    INTEGER :: PolyiP1, PolyiP2
    INTEGER :: IntiP1, IntiP2, last_IntiP1, last_IntiP2
    INTEGER :: oldE, Polyf, Intf
    INTEGER :: ie1, et, nleft, efrst

    INTEGER, DIMENSION(maxEperV,MaxCrnrs) :: EdgeofCrnr
    INTEGER, DIMENSION(maxEperF,maxFaces) :: EdgeofFace
    LOGICAL :: lError

    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Int = POLY3D_Init

    newCrnr = 0
    newEdge = 0

    ! Compute each of Poly's corners'
    ! signed distance from the plane.
    call Poly_Crnr_Plane_Dist (Poly, plane, CrnrDist)

    nCrnr   = Poly%numCrnr
    maxDist = MAXVAL( CrnrDist(1:nCrnr) )
    minDist = MINVAL( CrnrDist(1:nCrnr) )

    if (maxDist <= S0) then

      ! If all corners of Poly are behind the plane,
      ! return Poly as the intersection.
      Int = Poly
      RETURN

    ELSEIF (minDist >= S0) then

      ! If all corners of Poly are in front of the plane,
      ! return an empty Polyhedra as the intersection.
      RETURN

    end if

    ! Some of Poly is behind Plane and some of Poly
    ! is in-front of Plane. Intersect Poly with Plane
    ! and store the portion behind or on Plane in
    ! the intersection polyhedron, Int.

    ! Loop over Poly's corners. If they reside
    ! behind or on the plane, assign them
    ! corner storage in the intersection polygon.
    newC = 0
    do c = 1, Poly%numCrnr
      if (CrnrDist(c) <= S0) then
        call Store_Poly_Crnr (Int, Poly%Crnr(c), newC)
        newCrnr(c)     = newC
      end if
    end do
    Int%numCrnr = newC

    ! The corners that existed in Poly, also
    ! exist in Int from Corner(1:nOrigC).
    ! Intersection corners will be
    ! stored following these.
    nOrigC = newC

    ! Loop over Poly's edges. If they reside
    ! behind or on the plane, assign them
    ! edge storage in the intersection polygon.
    newE            = 0
    Int%nEdgeofCrnr = 0
    do e = 1, Poly%numEdges

      call Get_Edge_Crnrs (Poly, e, PolyiP1, PolyiP2)

      if (CrnrDist(PolyiP1) <= S0) then
        if (CrnrDist(PolyiP2) <= S0) then

          ! The entire edge is behind the plane.
          ! Store the entire edge in the intersection Poly.
          IntiP1     = newCrnr(PolyiP1)
          IntiP2     = newCrnr(PolyiP2)
          call Store_Poly_Edge (Int, IntiP1, IntiP2, iE)
          newE       = iE
          newEdge(e) = iE

        else if (CrnrDist(PolyiP1) < S0) then

          ! The edge intersects the plane. Compute the
          ! intersection point, add it to Int's list
          ! of corners and store this modified edge in Int.
          newE = newE + 1

          P1   = Poly%Crnr(PolyiP1)
          P2   = Poly%Crnr(PolyiP2)
          dP1  = CrnrDist(PolyiP1)
          dP2  = CrnrDist(PolyiP2)
          call Edge_IntwPlane (P1, dP1, P2, dP2, Vint)

          ! add the intersection point
          ! to Int's list of corners
          call Store_Poly_Crnr (Int, Vint, IntiP2)
          IntiP1     = newCrnr(PolyiP1)
          newC       = IntiP2

          ! Store the edge in Int's edge list.
          call Store_Poly_Edge (Int, IntiP1, IntiP2, iE)

          newE       = iE
          newEdge(e) = iE

        end if

      elseif (CrnrDist(PolyiP2) < S0) then

        ! The edge intersects the plane.
        ! Compute the intersection point,
        P1   = Poly%Crnr(PolyiP1)
        P2   = Poly%Crnr(PolyiP2)
        dP1  = CrnrDist(PolyiP1)
        dP2  = CrnrDist(PolyiP2)
        call Edge_IntwPlane (P1, dP1, P2, dP2, Vint)

        ! add the intersection point
        ! to Int's list of corners
        call Store_Poly_Crnr (Int, Vint, IntiP1)
        IntiP2     = newCrnr(PolyiP2)
        newC       = IntiP1

        ! Store the edge in Int's edge list.
        call Store_Poly_Edge (Int, IntiP1, IntiP2, iE)

        newE       = iE
        newEdge(e) = iE

      end if
    end do

    ! Examine each face of Poly. If all or part of
    ! the face lies behind the plane, store the
    ! proper portion of the face in Int.
    do Polyf = 1,Poly%numFaces

      ! If more than a single edge of face Polyf of Poly
      ! lies behind the plane, store the face in Int.
      ! Form a list of the edges for this facet
      ! that are on or behind the plane.
      nEdgeBehind   = 0
      facet_newEdge = 0
      do e = 1, Poly%nEdgeofFace(Polyf)
        oldE = Poly%EdgeofFace(e,Polyf)
        if (oldE > 0) then
          ie = newEdge(oldE)
        else
          ie = -newEdge(-oldE)
        end if
        if (ie /= 0) then
          nEdgeBehind                = nEdgeBehind + 1
          facet_newEdge(nEdgeBehind) = ie
        end if
      end do

      if (nEdgeBehind > 1) then
        ! This face comprises a part of Poly
        ! behind the plane. Associate the
        ! edge with this face in the
        ! intersection polyhedra, Int.
        ! If the plane cuts through the face,
        ! a new edge will be created. This
        ! new edge will extend from iPI1 to
        ! iPI2.

        ! Allocate the face in Int.
        Intf         = Int%numFaces + 1
        Int%numFaces = Intf
        if (Intf > maxFaces) then
          write (*,*) 'ERROR PLANE_POLY_INT3D: Too many faces'
          STOP
        end if

        ! We have formed a list of the edges for this face
        ! that are already stored in the polyhedra Int.
        ! These edges now will be associated with this
        ! facet of Int. There may be an intersection edge
        ! that cuts across the facet in the input polyhedra.
        ! This cut will become a new edge. To check for
        ! this cut edge, check that the starting and
        ! ending vertices are the same. If they are
        ! different, a cut edge should be constructed
        ! that connects these two vertices.
        last_newE = facet_newEdge(nEdgeBehind)
        call Get_Edge_Crnrs (Int,  last_newE, last_IntiP1, last_IntiP2)

        do e = 1, nEdgeBehind

          newE = facet_newEdge(e)    ! the edge in Poly

          if (newE == 0) then
            ! The edge is not part of the intersection.
            ! This is an error because the edges that
            ! are part of the intersection have been
            ! compressed in a loop above.
            write (*,*) 'ERROR PLANE_POLY_INT3D: Compressed list '// &
                        'facet_oldedge contains a zero'
            STOP
          end if

          call Get_Edge_Crnrs (Int,  newE, IntiP1,  IntiP2)

          ! Is there an intersection edge cutting across the
          ! Poly face?
          if (IntiP1 /= last_IntiP2) then
              call Store_Poly_Edge (Int, last_IntiP2, IntiP1, iE)
              call Store_Poly_Edge_of_Face (Int, iE, Intf)
          end if

          ! Store this edge
          call Store_Poly_Edge_of_Face (Int, newE, Intf)

          last_IntiP1 = IntiP1
          last_IntiP2 = IntiP2

        end do

      end if
    enddo

    ! If the plane intersected the polygon at other than a single point,
    ! there is an intersection face that needs to be constructed in Int.
    ! This face is composed of the unpaired edges in Int. For instance,
    ! most edges +iE and -iE will be a part of adjacent faces of Int.
    ! We will compose the intersection face by eliminating all paired
    ! edges. Only the intersection face will remain.
    EdgeofCrnr (1:maxEperV,1:maxCrnrs) = Int%EdgeofCrnr (1:maxEperV,1:maxCrnrs)
    EdgeofFace (1:maxEperF,1:maxFaces) = Int%EdgeofFace (1:maxEperF,1:maxFaces)
    do f = 1, Int%numFaces
      do e = 1, Int%nEdgeofFace(f)

        ! Grab an edge of Int that is part of face f
        ! and its endpoints.
        ie1 = EdgeofFace(e,f)
        call Get_Edge_Crnrs (Int, ie1, iP1, iP2)

        ! Delete the edge from the list of edges
        ! emanating from iP1
        lError = .TRUE.
        do et = 1, Int%nEdgeofCrnr(iP1)
          if (EdgeofCrnr(et,iP1) == ie1) then
            EdgeofCrnr(et,iP1) = 0
            lError             = .FALSE.
            EXIT
          end if
        end do

        if (lError) then
          write (*,*) 'ERROR PLANE_POLY_INT3D: facial edge not found at start corner'
          write (*,*) 'Input Poly is: '
          call Print_Poly (Poly)
          write (*,*) ' '
          write (*,*) 'Current Intersection Poly is:'
          call Print_Poly (Int)
          STOP
        end if

      end do
    end do

    do
      ! If any edges remain in the list, arrange their
      ! opposites into the intersection face.
      nleft = 0
      efrst = 0
      do c = 1, Int%numCrnr
        do e = 1, Int%nEdgeofCrnr(c)
          if (EdgeofCrnr(e,c) /= 0) then
            nleft = nleft + 1
            efrst = EdgeofCrnr(e,c)
          end if
        end do
      end do

      if (nleft == 0) EXIT

      f            = Int%numFaces + 1
      Int%numFaces = f
      if (f > maxFaces) then
        write (*,*) 'ERROR PLANE_POLY_INT3D: maxFaces too small'
        STOP
      end if
      ie1 = efrst
      call Store_Poly_Edge_of_Face (Int, ie1, f)
      do c = 1, nleft-1
        call Get_Edge_Crnrs (Int, ie1, iP1, iP2)
        do e = 1, INT%nEdgeofCrnr(iP2)
          if (EdgeofCrnr(e,iP2) /= 0) then
            ie1 = EdgeofCrnr(e,iP2)
            CYCLE
          end if
        end do
        if (ie1 /= efrst) then
          call Store_Poly_Edge_of_Face (Int, ie1, f)
        else
          EXIT
        endif
      enddo

      ! We've closed the face. Zero the edges we've used
      ! in constructing this face.
      do iE = 1, Int%nEdgeofFace(f)
        ie1 = Int%EdgeofFace(iE,f)
        call Get_Edge_Crnrs (Int, ie1, iP1, iP2)
        lError = .TRUE.
        do n = 1, Int%nEdgeofCrnr(iP1)
          if (EdgeofCrnr(n,iP1) == ie1) then
            EdgeofCrnr(n,iP1) = 0
            lError = .FALSE.
            EXIT
          endif
        enddo
        if (lError) then
          write (*,*) 'ERROR PLANE_POLY_INT3D: Lost the edge trail!'
          STOP
        endif
      enddo
    enddo

    END SUBROUTINE PLANE_POLY_INT3D

   SUBROUTINE Poly_Crnr_Plane_Dist (Poly, plane, CrnrDist)

  !==================================================================
  ! This routine is passed a three-dimensional Polyhedra, Poly, and
  ! a plane. It calculates the signed distance of each corner in
  ! Poly from the plane.
  !==================================================================

   IMPLICIT NONE

   TYPE(Poly3D_Type),               INTENT(IN)    :: Poly
   TYPE(Plane_Type),                INTENT(IN)    :: plane
   REAL(RK),          DIMENSION(:), INTENT(INOUT) :: CrnrDist

   !------------------------------------------------------------------

   INTEGER :: Size_of_Dist
   INTEGER :: c, nCrnr

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

   Size_of_Dist = SIZE( CrnrDist )
   if (Size_of_Dist < Poly%numCrnr) then
     write (*,*) 'ERROR Poly_Plane_Dist: not enough entries in CrnrDist'
     STOP
   end if

   CrnrDist(1:Size_of_Dist) = S0
   nCrnr                    = Poly%numCrnr
   do c = 1, nCrnr
     CrnrDist(c) = Point_Dist_Plane (Poly%Crnr(c), plane)
   end do

   END SUBROUTINE Poly_Crnr_Plane_Dist

   SUBROUTINE Edge_IntwPlane (P1, dP1, P2, dP2, Vint)

   IMPLICIT NONE

   TYPE(vector_type), INTENT(IN)  :: P1
   REAL(RK),          INTENT(IN)  :: dP1
   TYPE(vector_type), INTENT(IN)  :: P2
   REAL(RK),          INTENT(IN)  :: dP2
   TYPE(vector_type), INTENT(OUT) :: Vint

   !------------------------------------------------------------------

   REAL(RK) :: small_dist
   REAL(RK) :: dsign, denom, frac

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!   small_dist = S100 * SPACING( S1 )
   small_dist = 0.0_RK

   ! Avoid roundoff by picking the point closest to the plane.
   if (ABS( dP2 ) < ABS( dP1 )) then
     denom = dP1 - dP2
     dsign = SIGN( S1, denom )
     denom = MAX( ABS(denom), small_dist )
     denom = dsign * denom
     frac  = -dP2 / denom
     Vint  = P2 + frac * (P1 - P2)
   else
     denom = dP2 - dP1
     dsign = SIGN( S1, denom )
     denom = MAX( ABS(denom), small_dist )
     denom = dsign * denom
     frac  = -dP1 / denom
     Vint  = P1 + frac * (P2 - P1)
   end if

   END SUBROUTINE Edge_IntwPlane

  TYPE(Poly3D_Type) FUNCTION Load_Tet (Th1, Th2, Th3, Th4)

  TYPE(VECTOR_TYPE), INTENT(IN) :: Th1, Th2, Th3, Th4

  Load_Tet = Poly3D_Init

  Load_Tet%numCrnr = 4
  Load_Tet%Crnr(1) = Th1
  Load_Tet%Crnr(2) = Th2
  Load_Tet%Crnr(3) = Th3
  Load_Tet%Crnr(4) = Th4

  Load_Tet%numFaces = 4
  Load_Tet%numEdges = 6

  Load_Tet%nEdgeofFace(1:4) = 3

  Load_Tet%EdgeofFace(1:3,1) = (/  1,  5, -4 /)
  Load_Tet%EdgeofFace(1:3,2) = (/  2,  6, -5 /)
  Load_Tet%EdgeofFace(1:3,3) = (/  3,  4, -6 /)
  Load_Tet%EdgeofFace(1:3,4) = (/ -1, -3, -2 /)

  Load_Tet%CrnrofEdge(1:2,1) = (/ 1, 2 /)
  Load_Tet%CrnrofEdge(1:2,2) = (/ 2, 3 /)
  Load_Tet%CrnrofEdge(1:2,3) = (/ 3, 1 /)
  Load_Tet%CrnrofEdge(1:2,4) = (/ 1, 4 /)
  Load_Tet%CrnrofEdge(1:2,5) = (/ 2, 4 /)
  Load_Tet%CrnrofEdge(1:2,6) = (/ 3, 4 /)

  Load_Tet%nEdgeofCrnr(1:4)  = 3

  Load_Tet%EdgeofCrnr(1:3,1) = (/  1, -3,  4 /)
  Load_Tet%EdgeofCrnr(1:3,2) = (/  2, -1,  5 /)
  Load_Tet%EdgeofCrnr(1:3,3) = (/  3, -2,  6 /)
  Load_Tet%EdgeofCrnr(1:3,4) = (/ -4, -5, -6 /)

  END FUNCTION Load_Tet

  TYPE(Poly3D_Type) FUNCTION Load_Hex (Hx1, Hx2, Hx3, Hx4, &
                                       Hx5, Hx6, Hx7, Hx8)

  TYPE(VECTOR_TYPE), INTENT(IN) :: Hx1, Hx2, Hx3, Hx4
  TYPE(VECTOR_TYPE), INTENT(IN) :: Hx5, Hx6, Hx7, Hx8

  Load_Hex = Poly3D_Init

  Load_Hex%numCrnr = 8

  Load_Hex%Crnr(1) = Hx1
  Load_Hex%Crnr(2) = Hx2
  Load_Hex%Crnr(3) = Hx3
  Load_Hex%Crnr(4) = Hx4
  Load_Hex%Crnr(5) = Hx5
  Load_Hex%Crnr(6) = Hx6
  Load_Hex%Crnr(7) = Hx7
  Load_Hex%Crnr(8) = Hx8

  Load_Hex%numFaces =  6
  Load_Hex%numEdges = 12

  Load_Hex%nEdgeofFace(1:6) = 4

  Load_Hex%EdgeofFace(1:4,1) = (/   3,   2,  -4,  -1 /)
  Load_Hex%EdgeofFace(1:4,2) = (/   5,   8,  -6,  -7 /)
  Load_Hex%EdgeofFace(1:4,3) = (/   1,  11,  -5,  -9 /)
  Load_Hex%EdgeofFace(1:4,4) = (/  10,   6, -12,  -2 /)
  Load_Hex%EdgeofFace(1:4,5) = (/   9,   7, -10,  -3 /)
  Load_Hex%EdgeofFace(1:4,6) = (/   4,  12,  -8, -11 /)

  Load_Hex%CrnrofEdge(1:2, 1) = (/ 1, 2 /)
  Load_Hex%CrnrofEdge(1:2, 2) = (/ 4, 3 /)
  Load_Hex%CrnrofEdge(1:2, 3) = (/ 1, 4 /)
  Load_Hex%CrnrofEdge(1:2, 4) = (/ 2, 3 /)
  Load_Hex%CrnrofEdge(1:2, 5) = (/ 5, 6 /)
  Load_Hex%CrnrofEdge(1:2, 6) = (/ 8, 7 /)
  Load_Hex%CrnrofEdge(1:2, 7) = (/ 5, 8 /)
  Load_Hex%CrnrofEdge(1:2, 8) = (/ 6, 7 /)
  Load_Hex%CrnrofEdge(1:2, 9) = (/ 1, 5 /)
  Load_Hex%CrnrofEdge(1:2,10) = (/ 4, 8 /)
  Load_Hex%CrnrofEdge(1:2,11) = (/ 2, 6 /)
  Load_Hex%CrnrofEdge(1:2,12) = (/ 3, 7 /)

  Load_Hex%nEdgeofCrnr( 1:8) = 3

  Load_Hex%EdgeofCrnr(1:3,1) = (/   1,   3,   9 /)
  Load_Hex%EdgeofCrnr(1:3,2) = (/  -1,   4,  11 /)
  Load_Hex%EdgeofCrnr(1:3,3) = (/  -4,  -2,  12 /)
  Load_Hex%EdgeofCrnr(1:3,4) = (/   2,  -3,  10 /)
  Load_Hex%EdgeofCrnr(1:3,5) = (/   5,   7,  -9 /)
  Load_Hex%EdgeofCrnr(1:3,6) = (/  -5,   8, -11 /)
  Load_Hex%EdgeofCrnr(1:3,7) = (/  -6,  -8, -12 /)
  Load_Hex%EdgeofCrnr(1:3,8) = (/  -7,   6, -10 /)

  END FUNCTION Load_Hex

  SUBROUTINE Get_Edge_Crnrs (Poly, iE, iP1, iP2)

  IMPLICIT NONE

  !==================================================================
  ! A Three-Dimensional Polyhedra, Poly, is passed to this routine.
  ! iE is an edge index in Poly. This routine returns the starting
  ! and ending corner indices of the edge.
  !==================================================================

  TYPE(Poly3D_Type), INTENT(IN)  :: Poly
  INTEGER,           INTENT(IN)  :: iE
  INTEGER,           INTENT(OUT) :: iP1, iP2

  !------------------------------------------------------------------
  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  if (ABS( iE ) > Poly%numEdges .OR. iE == 0) then
    write (*,*) 'ERROR GET_EDGE_CRNRS: illegal edge number passed in'
    STOP
  end if

  if (iE > 0) then
    iP1 = Poly%CrnrofEdge(1,iE)
    iP2 = Poly%CrnrofEdge(2,iE)
  else
    iP1 = Poly%CrnrofEdge(2,-iE)
    iP2 = Poly%CrnrofEdge(1,-iE)
  end if

  END SUBROUTINE Get_Edge_Crnrs

  SUBROUTINE Store_Poly_Crnr (Poly, C, iP)

  IMPLICIT NONE

  !==================================================================
  ! A Three-Dimensional Polyhedra, Poly, is passed to this routine.
  ! The position vector of acorner of the polyhedra is passed in
  ! as C. This routine stores the new Corner in Poly and returns
  ! its corner index in Poly as iP.
  !
  ! It is assumed that the corner DOES NOT currently exist in Poly.
  ! If it is possible that the corner already exists in Poly, a more
  ! extensive routine should be used that will check for the it's
  ! existence before storing it.
  !
  !==================================================================

  TYPE(Poly3D_Type), INTENT(INOUT) :: Poly
  TYPE(Vector_Type), INTENT(IN)    :: C
  INTEGER,           INTENT(OUT)   :: iP

  !------------------------------------------------------------------

  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  iP = Poly%numCrnr + 1
  if (iP > maxCrnrs) then
    write (*,*) 'ERROR STORE_POLY_CRNR: Parameter maxCrnrs too small'
    STOP
  end if

  Poly%numCrnr = iP

  Poly%Crnr(iP) = C

  END SUBROUTINE Store_Poly_Crnr

  SUBROUTINE Store_Poly_Edge (Poly, iP1, iP2, iE)

  IMPLICIT NONE

  !==================================================================
  ! A Three-Dimensional Polyhedra, Poly, is passed to this routine.
  ! Two edge corner indices are also passed. The edge starting
  ! index is iP1 and the ending index is iP2. This routine stores
  ! the edge in Poly and returns its edge index, iE.
  !
  !==================================================================

  TYPE(Poly3D_Type), INTENT(INOUT) :: Poly
  INTEGER,           INTENT(IN)    :: iP1
  INTEGER,           INTENT(IN)    :: iP2
  INTEGER,           INTENT(OUT)   :: iE

  !------------------------------------------------------------------

  INTEGER :: e, numE, tiP1, tiP2

  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  if (iP1 > Poly%numCrnr .OR. iP1 < 0) then
    write (*,*) 'ERROR STORE_POLY_EDGE: illegal starting corner index'
    STOP
  end if

  if (iP2 > Poly%numCrnr .OR. iP2 < 0) then
    write (*,*) 'ERROR STORE_POLY_EDGE: illegal ending corner index'
    STOP
  end if

  ! The edge may already exist in Poly. An economical test
  ! is to look at the edges associated with corner iP1.
  ! If any of these edges stop at iP2, return this edge
  ! as iE and don't re-store the edge in Poly.
  do e = 1, Poly%nEdgeofCrnr(iP1)

    iE = Poly%EdgeofCrnr(e,iP1)

    call Get_Edge_Crnrs (Poly, iE, tiP1, tiP2)

    if (tiP2 == iP2) then
      RETURN
    endif
  end do

  ! The edge doens't currently exist in Poly,
  ! create it.
  iE = Poly%numEdges + 1
  if (iE > maxEdges) then
    write (*,*) 'ERROR STORE_POLY_EDGE: parameter maxEdges too small'
    STOP
  end if

  Poly%numEdges = iE

  ! Store the edge in Poly's edge list.
  Poly%CrnrofEdge(1,iE) = iP1
  Poly%CrnrofEdge(2,iE) = iP2

  ! Associate the edge with the
  ! starting corner of the edge.
  numE = Poly%nEdgeofCrnr(iP1) + 1
  if (numE > maxEperV) then
    write (*,*) 'ERROR STORE_POLY_EDGE: parameter maxEperV too small'
    STOP
  end if
  Poly%nEdgeofCrnr(iP1)     = numE
  Poly%EdgeofCrnr(numE,iP1) = iE

  ! Associate the edge with the
  ! ending corner of the edge.
  numE = Poly%nEdgeofCrnr(iP2) + 1
  if (numE > maxEperV) then
    write (*,*) 'ERROR STORE_POLY_EDGE: parameter maxEperV too small'
    STOP
  end if
  Poly%nEdgeofCrnr(iP2)     = numE
  Poly%EdgeofCrnr(numE,iP2) = -iE

  END SUBROUTINE Store_Poly_Edge

  SUBROUTINE Store_Poly_Edge_of_Face (Poly, e, f)

  IMPLICIT NONE

  !==================================================================
  ! A Three-Dimensional Polyhedra, Poly, is passed to this routine.
  ! It is also passed an edge index, e, and a face index, f, for
  ! for entries in Poly. The edge e is stored as the 'next' edge
  ! in face, f.
  !
  !==================================================================

  TYPE(Poly3D_Type), INTENT(INOUT) :: Poly
  INTEGER,           INTENT(IN)    :: e
  INTEGER,           INTENT(IN)    :: f

  !------------------------------------------------------------------

  INTEGER :: n

  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  if (Poly%numEdges < ABS( e )) then
    write (*,*) 'ERROR STORE_POLY_EDGE_OF_FACE: index of edge ',  &
                'passed is > # of edges in Poly'
    STOP
  end if
  if (Poly%numFaces < f) then
    write (*,*) 'ERROR STORE_POLY_EDGE_OF_FACE: index of face ',  &
                'passed is > # of faces in Poly'
    STOP
  end if
  if (f < 0) then
    write (*,*) 'ERROR STORE_POLY_EDGE_OF_FACE: index of face ',  &
                'passed is < 0'
    STOP
  end if

  ! Determine if the edge is already stored as part of this face.
  n = Poly%nEdgeofFace(f)
  if (n > 0) then
    if (Poly%EdgeofFace(1,f) == e .OR. &
        Poly%EdgeofFace(n,f) == e      ) then
      RETURN
    endif
  end if

  n = Poly%nEdgeofFace(f) + 1
  if (n > maxEperF) then
    write (*,*) 'ERROR STORE_POLY_EDGE_OF_FACE: too many edges ',  &
                'stored on face'
    STOP
  endif

  Poly%nEdgeofFace(  f) = n
  Poly%EdgeofFace (n,f) = e

  END SUBROUTINE Store_Poly_Edge_of_Face

  REAL(RK) FUNCTION volm_poly3D (Poly)

  !==================================================================
  ! Compute the volume of the polygon, Poly
  !==================================================================

  IMPLICIT NONE

  TYPE(POLY3D_TYPE), INTENT(IN)  :: Poly

  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  volm_poly3D = volm_poly_3d_cart (Poly)

  END FUNCTION volm_poly3D

  REAL(RK) FUNCTION volm_poly_3d_cart (Poly)

  !==================================================================
  ! Compute the volume of the polygon for 3D in Cartesian coordinates
  !==================================================================

  IMPLICIT NONE

  TYPE(POLY3D_TYPE), INTENT(IN)  :: Poly

  !------------------------------------------------------------------

  INTEGER :: f, e, ie, iP1, iP2
  TYPE(VECTOR_TYPE) :: Pm, Fm, Fa, P1, P2
  REAL(RK) :: Vol

  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  volm_poly_3d_cart = 0.0_RK

  if (Poly%numCrnr > 3) then

    Pm = Poly%Crnr(1)

    Vol = 0.0_RK
    do f = 1, Poly%numFaces

      if (Poly%nEdgeofFace(f) > 2) then
        ie = Poly%EdgeofFace(1,f)
        call Get_Edge_Crnrs (Poly, ie, iP1, iP2)
        Fm = Poly%Crnr(iP1)

        Fa = ORIGIN
        do e = 2, Poly%nEdgeofFace(f)

          ie = Poly%EdgeofFace(e,f)
          call Get_Edge_Crnrs (Poly, ie, iP1, iP2)

          P1 = Poly%Crnr(iP1)
          P2 = Poly%Crnr(iP2)

          Fa = Fa + Vector_Cross_Prod (P1-Fm, P2-Fm)
        end do
        Vol = Vol + Fa*(Fm - Pm)

      end if
    end do

    volm_poly_3d_cart = Vol / 6.0_RK

  endif

  END FUNCTION volm_poly_3d_cart

  SUBROUTINE Poly3D_2_Poly3D (Trg, Src)

  !==================================================================
  !  This subroutine places the contents of the Polygon Src into
  !  the Polygon Trg.
  !==================================================================

  IMPLICIT NONE

  TYPE(POLY3D_TYPE), INTENT(OUT) :: Trg
  TYPE(POLY3D_TYPE), INTENT(IN)  :: Src

  INTEGER :: c, f, e

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  Trg%numCrnr = Src%numCrnr
  do c = 1,maxCrnrs
    Trg%Crnr(c) = Src%Crnr(c)
  enddo
  Trg%numFaces = Src%numFaces
  Trg%numEdges = Src%numEdges
  do f = 1, maxFaces
    Trg%nEdgeofFace(f) = Src%nEdgeofFace(f)
    do e = 1, maxEperF
      Trg%EdgeofFace(e,f) = Src%EdgeofFace(e,f)
    end do
  enddo
  do e = 1, maxEdges
    Trg%CrnrofEdge(1,e) = Src%CrnrofEdge(1,e)
    Trg%CrnrofEdge(2,e) = Src%CrnrofEdge(2,e)
  end do
  do c = 1, maxCrnrs
    Trg%nEdgeofCrnr(c) = Src%nEdgeofCrnr(c)
    do e = 1, maxEperV
      Trg%EdgeofCrnr(e,c) = Src%EdgeofCrnr(e,c)
    end do
  end do

  END SUBROUTINE Poly3D_2_Poly3D

  REAL(RK) FUNCTION Point_Dist_Plane (point, plane)
  
  !==================================================================
  ! This logical function determines if the location in point is
  ! physically behind the plane.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: point
  TYPE(PLANE_TYPE),  INTENT(IN) :: plane

  !------------------------------------------------------------------

  INTEGER :: small_dist

  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  Point_Dist_Plane = Point * Plane%n - Plane%p

!  small_dist = S100 * SPACING( S1 )
  small_dist = 0.0_RK

  if (ABS( Point_Dist_Plane ) < small_dist) then
    Point_Dist_Plane = S0
  end if

  END FUNCTION Point_Dist_Plane

  TYPE(PLANE_TYPE) FUNCTION Plane_from_Points (Pt1, Pt2, Pt3)

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: Pt1
  TYPE(VECTOR_TYPE), INTENT(IN) :: Pt2
  TYPE(VECTOR_TYPE), INTENT(IN) :: Pt3

  !------------------------------------------------------------------

  TYPE(VECTOR_TYPE) :: nrml
  REAL(RK) :: p

  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  nrml              = Vector_Cross_Prod (Pt3 - Pt2, Pt1 - Pt2)
  nrml              = nrmlz_vector (nrml)
  p                 = nrml * Pt1
  Plane_from_Points = PLANE_TYPE ( nrml, p )

  END FUNCTION Plane_from_Points

  TYPE(VECTOR_TYPE) FUNCTION NEGATE_VECTOR (A)

  !==================================================================
  !  This function returns a directed line segment by negating
  !  vector A.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A

  !------------------------------------------------------------------

  INTEGER :: d

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  do d = 1,mdim
    negate_vector%x(d) = -A%x(d)
  enddo

  END FUNCTION NEGATE_VECTOR

  FUNCTION NEGATE_VECTOR_ARRAY (A) RESULT(B)

  !==================================================================
  !  This function returns an array of vectors by negating an array
  !  of vectors A.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), DIMENSION(:),      INTENT(IN) :: A
  TYPE(VECTOR_TYPE), DIMENSION(SIZE(A))            :: B

  !------------------------------------------------------------------

  INTEGER :: l

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  do l = 1,SIZE(A)
    B(l) = -A(l)
  enddo

  END FUNCTION NEGATE_VECTOR_ARRAY

  TYPE(VECTOR_TYPE) FUNCTION SUBTRACT_VECTOR (A, B)

  !==================================================================
  !  This function returns a directed line segment by subtracting
  !  vector B from vector A.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A, B

  !------------------------------------------------------------------

  INTEGER :: d

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  do d = 1,mdim
     subtract_vector%x(d) = A%x(d) - B%x(d)
  enddo

  END FUNCTION SUBTRACT_VECTOR

  TYPE(VECTOR_TYPE) FUNCTION ADD_VECTOR (A, B)

  !==================================================================
  !  This function returns a directed line segment by adding
  !  point B to point A.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A, B

  !------------------------------------------------------------------

  INTEGER :: d

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  do d = 1,mdim
     add_vector%x(d) = A%x(d) + B%x(d)
  enddo

  END FUNCTION ADD_VECTOR

  TYPE(VECTOR_TYPE) FUNCTION scalar_mult_vector (A, B)

  !==================================================================
  !  This function returns a directed line segment by multiplying
  !  the scalar A with the vector B.
  !==================================================================

  IMPLICIT NONE

  REAL(RK),          INTENT(IN) :: A
  TYPE(VECTOR_TYPE), INTENT(IN) :: B

  !------------------------------------------------------------------

  INTEGER :: d

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  do d = 1,mdim
     scalar_mult_vector%x(d) = A * B%x(d)
  enddo

  END FUNCTION scalar_mult_vector

  TYPE(VECTOR_TYPE) FUNCTION vector_mult_scalar (A, B)

  !==================================================================
  !  This function returns a directed line segment by multiplying
  !  the scalar B with the vector A.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A
  REAL(RK),          INTENT(IN) :: B

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  vector_mult_scalar = scalar_mult_vector (B, A)

  END FUNCTION vector_mult_scalar

  REAL(RK) FUNCTION vector_dot_vector (A, B)

  !==================================================================
  !  This function returns the dot product of vector A with 
  !  vector B.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A
  TYPE(VECTOR_TYPE), INTENT(IN) :: B

  !------------------------------------------------------------------

  REAL(RK) :: dot
  INTEGER :: n

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  dot = 0.0_RK
  do n = 1,mdim
     dot = dot + A%x(n)*B%x(n)
  enddo
  vector_dot_vector = dot

  END FUNCTION vector_dot_vector

  TYPE(VECTOR_TYPE) FUNCTION vector_div_scalar (A, B)

  !==================================================================
  !  This function returns a directed line segment by dividing
  !  the vector A by the scalar B.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A
  REAL(RK),          INTENT(IN) :: B

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  vector_div_scalar = &
               scalar_mult_vector (S1/SIGN(MAX(ABS(B),small),B), A)

  END FUNCTION vector_div_scalar

  TYPE(VECTOR_TYPE) FUNCTION abs_of_vector (A)

  !==================================================================
  !  This function returns the absolute value of each
  !  component of the vector A.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A

  !------------------------------------------------------------------

  INTEGER :: d

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  do d = 1,mdim
     abs_of_vector%x(d) = ABS (A%x(d))
  enddo

  END FUNCTION abs_of_vector

  REAL(RK) FUNCTION vector_lngth (A)

  !==================================================================
  !  This function returns the length of the vector A.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A

  !------------------------------------------------------------------

  INTEGER :: d
  REAL(RK) :: vector_sum

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  vector_sum = 0.0_RK
  do d = 1,mdim
     vector_sum = vector_sum + A%x(d)**2
  enddo
  vector_lngth = SQRT (vector_sum)

  END FUNCTION vector_lngth

  TYPE(VECTOR_TYPE) FUNCTION nrmlz_vector (A)

  !==================================================================
  !  This function is passed a vector A and this routine normalizes
  !  the components.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A

  !------------------------------------------------------------------

  INTEGER :: d
  REAL(RK) :: Length

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  Length = 1.0_RK / MAX (vector_lngth(A), small)

  do d = 1,mdim
     nrmlz_vector%x(d) = A%x(d) * Length
  enddo

  END FUNCTION nrmlz_vector

  FUNCTION vector_cross_prod (A, B) RESULT (C)

  !==================================================================
  !  This function returns the vector cross product of A x B.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(IN) :: A, B
  TYPE(VECTOR_TYPE)             :: C

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  C%x(1) =  A%x(2) * B%x(3) - A%x(3) * B%x(2)
  C%x(2) =  A%x(3) * B%x(1) - A%x(1) * B%x(3)
  C%x(3) =  A%x(1) * B%x(2) - A%x(2) * B%x(1)

  END FUNCTION vector_cross_prod

  SUBROUTINE Vector_2_Vector (Trg, Src)

  !==================================================================
  !  This subroutine places the contents of the vector Src into
  !  the vector Trg.
  !==================================================================

  IMPLICIT NONE

  TYPE(VECTOR_TYPE), INTENT(OUT) :: Trg
  TYPE(VECTOR_TYPE), INTENT(IN)  :: Src

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  Trg%x(1:mdim) = Src%x(1:mdim)

  END SUBROUTINE Vector_2_Vector

  SUBROUTINE Print_Poly_Mathematica (Poly)

   IMPLICIT NONE

   !===================================================================
   ! This routine prints the poly_type description of a poly.
   !===================================================================

   TYPE(poly3D_type), INTENT(IN) :: Poly

   !-------------------------------------------------------------------

   INTEGER :: face, edge, e, iv1, iv2

   ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

   write (*,*) ' '
   do face = 1, Poly%numFaces
     write (*,*) 'Polygon[{'
     do edge = 1, Poly%nEdgeofFace(face)
       e = Poly%EdgeofFace(edge,face)
       call Get_Edge_Crnrs (Poly, e, iv1, iv2)
       write (*,*) '{',Poly%Crnr(iv1)%x(1),',', Poly%Crnr(iv1)%x(2),',', Poly%Crnr(iv1)%x(3),'},'
     enddo
     write (*,*) '}],'
   enddo

  END SUBROUTINE Print_Poly_Mathematica

  SUBROUTINE Check_Poly_Edges (Poly, Plane)

  IMPLICIT NONE

  !===================================================================
  ! This routine prints the poly_type description of a poly.
  !===================================================================

  TYPE(poly3D_type), INTENT(IN) :: Poly
  TYPE(plane_type),  INTENT(IN) :: Plane

  !-------------------------------------------------------------------

  INTEGER :: face, edge, e
  INTEGER, DIMENSION(maxEperF,maxFaces) :: tmp_EdgeofFace
  LOGICAL, DIMENSION(-Poly%numEdges:Poly%numEdges) :: Possible_Edges
  LOGICAL :: Error_Found

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! ALLOCATE (Possible_Edges(-Poly%numEdgesPoly%:numEdges))

  ! Copy the poly's edges to a temporary array
  tmp_EdgeofFace = 0
  do face = 1, Poly%numFaces
    do edge = 1, Poly%nEdgeofFace(face)
      tmp_EdgeofFace(edge,face) = Poly%EdgeofFace(edge,face)
    enddo
  enddo

  ! Initialize the list of possible edges
  do edge = 1, Poly%numEdges
    Possible_Edges(-edge) = .TRUE.
    Possible_Edges( edge) = .TRUE.
  enddo

  ! Look for each edge between -numEdges:-1 and 1:numEdges
  ! Zero each edge once.
  Error_Found = .FALSE.
  do face = 1, Poly%numFaces
    do edge = 1, Poly%nEdgeofFace(face)
      e = Poly%EdgeofFace(edge,face)
      if (e < -Poly%numEdges .OR. &
          e == 0             .OR. &
          e >  Poly%numEdges      ) then
        write (*,*) 'ERROR Check_Poly_Edges: illegal edge found'
        Error_Found = .TRUE.
      endif
      if (Possible_Edges(e)) then
        Possible_Edges(e) = .FALSE.
      else
        write (*,*) 'ERROR Check_Poly_Edges: Edge ',e,' used multiple times'
        Error_Found = .TRUE.
      endif
    enddo
  enddo

  do edge = -Poly%numEdges,-1
    if (Possible_Edges(edge)) then
      Error_Found = .TRUE.
      write (*,*) 'ERROR Check_Poly_Edges: edge ',edge,' not used'
    endif
  enddo

  do edge = 1, Poly%numEdges
    if (Possible_Edges(edge)) then
      Error_Found = .TRUE.
      write (*,*) 'ERROR Check_Poly_Edges: edge ',edge,' not used'
    endif
  enddo

  if (Error_Found) then
    call Print_Poly (Poly)
    STOP
  endif
   
  END SUBROUTINE Check_Poly_Edges

  SUBROUTINE Print_Poly (Poly)

   IMPLICIT NONE

   !===================================================================
   ! This routine prints the poly_type description of a poly.
   !===================================================================

   TYPE(poly3D_type), INTENT(IN) :: Poly

   !-------------------------------------------------------------------

   INTEGER  :: crnr, face, edge, numEdges
   INTEGER  :: e1, e2, iv1, iv2, iv3
   REAL(RK) :: FaceArea, SubArea
   TYPE(Vector_type) :: Pvt, V1, V2, V3, A123, Ann

   TYPE(Vector_type) :: Pt9, Pt13, Pt16
   TYPE(Vector_type) :: Pt4, Pt14, Pt15
   TYPE(Plane_type)  :: Face11, Face12

   ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

   write (*,*)
   write (*,*) 'Poly has ',Poly%numCrnr,' corners'
   do crnr = 1, Poly%numCrnr
     write (*,5) ' Crnr(',crnr,') = {',Poly%Crnr(crnr)%x(1), &
                               ', ',   Poly%Crnr(crnr)%x(2),'}'
5    format (a,i2,a,1pe12.5,a,1pe12.5,a)
   enddo

   do edge = 1, Poly%numEdges
     write (*,10) ' Edge(',edge,') = {',Poly%CrnrofEdge(1,edge), &
                                ', ',   Poly%CrnrofEdge(2,edge),'}'
10   format (a,i3,a,i4,a,i4,a)
   enddo

   write (*,*) 'Poly Edges:'
   do face = 1, Poly%numFaces
     write (*,*) 'Polygon face: ', face
     numEdges = Poly%nEdgeofFace(face)
     write (*,15) '        edges: ',Poly%EdgeofFace(1:numEdges,face)
15   format (a, 10i4)
     FaceArea = 0.0_RK
     e1 = Poly%EdgeofFace(1,face)
     e2 = Poly%EdgeofFace(2,face)
     call Get_Edge_Crnrs (Poly, e1, iV1, iV2)
     call Get_Edge_Crnrs (Poly, e2, iV2, iV3)
     V1 = Poly%Crnr(iV1)
     V2 = Poly%Crnr(iV2)
     V3 = Poly%Crnr(iV3)
     FaceArea =  Vector_Lngth( Vector_Cross_Prod (V2-V1, V3-V1) )
     if (Poly%nEdgeofFace(face) > 3) then
       Pvt  = V1
       A123 = Vector_Cross_Prod (V2-V1, V3-V1)
       A123 = nrmlz_vector (A123)
       do edge = 3, Poly%nEdgeofFace(face)-1
         e1       = Poly%EdgeofFace(edge,face)
         call Get_Edge_Crnrs (Poly, e1, iV1, iV2)
         V1       = Poly%Crnr(iV1)
         V2       = Poly%Crnr(iV2)
         Ann      = Vector_Cross_Prod (V1-Pvt, V2-Pvt)
         SubArea  = Vector_Lngth( Ann )
         FaceArea = FaceArea + SubArea
         Ann      = nrmlz_vector (Ann)
         write (*,*) 'Face Sub Dot = ',A123*Ann,' SubArea = ',SubArea   
       enddo
     endif
     write (*,*) 'Face area = ',FaceArea
   enddo

   Face11 = Plane_from_Points (Pt9, Pt13, Pt16)
   Face12 = Plane_from_Points (Pt15, Pt14, Pt4)

   write (*,*) 'Pt9  dist = ',Point_Dist_Plane (Pt9,  Face12)
   write (*,*) 'Pt13 dist = ',Point_Dist_Plane (Pt13, Face12)
   write (*,*) 'Pt16 dist = ',Point_Dist_Plane (Pt16, Face12)

   write (*,*) 'Pt4  dist = ',Point_Dist_Plane (Pt4,  Face11)
   write (*,*) 'Pt14 dist = ',Point_Dist_Plane (Pt14, Face11)
   write (*,*) 'Pt15 dist = ',Point_Dist_Plane (Pt15, Face11)

  END SUBROUTINE Print_Poly

END MODULE OVERLAP_MODULE
