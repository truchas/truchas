MODULE MESH_QUALITY_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Ensure hex mesh is legitimate using a criterion based on the
  !   Jacobian determinant criterion.
  !
  !   Public Interface:
  !
  !     * call TWISTED_CELL_TEST ()
  !
  !         Checks positivity of Jacobian determinant, whose scaled value
  !         may be used as a measure for hex mesh quality
  !
  ! Contains: TWISTED_CELL_TEST
  !
  ! Author(s): Kin Lam, LANL ESA-EA (klam@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: TWISTED_CELL_TEST

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE TWISTED_CELL_TEST ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate Jacobian volume of 24 tetrahedra associated with 12 edges
    !   (right-handed and left-handed orientation) of hexahedron cells
    !
    !   The volume of a tetrahedron is 1/6 that of the Jacobian determinant
    !
    !      det(J) = ( e1 cross e2 ) dot e3
    !
    !   where e1, e2, and e3 are three base/edge vectors defined by the four
    !   vertices of the tetrahedron.
    !
    !   For an untangled mesh, det(J) must be positive.
    !
    !   A common criterion used for hex mesh quality is
    !
    !      0.5 < scaled det(J) = det(J) / |e1|*|e2|*|e3| < 1.0
    !
    !   with scaled det(J) = 1.0 being the ideal case where the three
    !   edge vectors are orthogonal to each other.
    !
    !   References
    !   ----------
    !
    !   Handbook of Grid Generation, ed. Joe F. Thompson, Bharat K. Suni,
    !   Nigel P. Weatherill, CRC Press, 1999, Chapter 8, p. 28.
    !
    !   Patrick M. Knupp, Int. J. Numer. Meth. Engng. 2000; 48:1165-1185
    !
    !=======================================================================
    use mesh_module,          only: Vertex
    use gs_module,            only: EN_GATHER
    use parameter_module,     only: ncells, ndim, nvc
    use PGSLIB_Module,        only: pgslib_global_any

    ! Local variables
    integer :: c, edge, i
    integer :: vr1, vr2, vr3, vr4     ! vertices of right-handed tet
    integer :: vl1, vl2, vl3, vl4     ! vertices of left-handed tet
    integer :: v1, v2, v3, v4, v5, v6, v7, v8

    real(r8)   :: det_r, det_l
    real(r8)   :: e1e2e3_r, e1e2e3_l
    real(r8) :: jac_scaled_r
    real(r8) :: jac_scaled_l

    real(r8), dimension(3)   :: er1, er2, er3    ! edge/basis vectors for right-handed tet
    real(r8), dimension(3)   :: el1, el2, el3    ! edge/basis vectors for left-handed tet
    real(r8), dimension(3,3) :: jac_r, jac_l

    real(r8), pointer, dimension(:,:,:) :: Coord
    logical, pointer, dimension(:,:) :: twisted
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Do this only for hex (3D) mesh
    if (ndim /= 3) return

    ! Set vertex numbers
    v1 = 1; v2 = 2; v3 = 3; v4 = 4; v5 = 5; v6 = 6; v7 = 7; v8 = 8

    ! allocate temporary arrays
    allocate(coord(ndim,nvc,ncells), stat=status)
    if (status /= 0) call TLS_panic ('TWISTED_CELL_TEST: allocation error: coord')
    allocate(twisted(ncells,12), stat=status)
    if (status /= 0) call TLS_panic ('TWISTED_CELL_TEST: allocation error: twisted')
    twisted = .false.

    ! Gather vertex coordinates
    do i = 1, ndim
       call EN_GATHER (Coord(i,:,:), Vertex%Coord(i))
    end do

    ! Loop through edges of all cells
    do c = 1, ncells
       do edge = 1, 12
          select case (edge)
             case(1)                   ! edge 12
                vr1 = v4
                vr2 = v1
                vr3 = v2
                vr4 = v6
                vl1 = v5
                vl2 = v1
                vl3 = v2
                vl4 = v3
             case(2)                   ! edge 23
                vr1 = v1
                vr2 = v2
                vr3 = v3
                vr4 = v7
                vl1 = v6
                vl2 = v2
                vl3 = v3
                vl4 = v4
             case(3)                   ! edge 34
                vr1 = v2
                vr2 = v3
                vr3 = v4
                vr4 = v8
                vl1 = v7
                vl2 = v3
                vl3 = v4
                vl4 = v1
             case(4)                   ! edge 41
                vr1 = v3
                vr2 = v4
                vr3 = v1
                vr4 = v5
                vl1 = v8
                vl2 = v4
                vl3 = v1
                vl4 = v2
             case(5)                   ! edge 56
                vr1 = v1
                vr2 = v5
                vr3 = v6
                vr4 = v7
                vl1 = v8
                vl2 = v5
                vl3 = v6
                vl4 = v2
             case(6)                   ! edge 67
                vr1 = v2
                vr2 = v6
                vr3 = v7
                vr4 = v8
                vl1 = v5
                vl2 = v6
                vl3 = v7
                vl4 = v3
             case(7)                   ! edge 78
                vr1 = v3
                vr2 = v7
                vr3 = v8
                vr4 = v5
                vl1 = v6
                vl2 = v7
                vl3 = v8
                vl4 = v4
             case(8)                   ! edge 85
                vr1 = v4
                vr2 = v8
                vr3 = v5
                vr4 = v6
                vl1 = v7
                vl2 = v8
                vl3 = v5
                vl4 = v1
             case(9)                   ! edge 15
                vr1 = v2
                vr2 = v1
                vr3 = v5
                vr4 = v8
                vl1 = v4
                vl2 = v1
                vl3 = v5
                vl4 = v6
             case(10)                  ! edge 26
                vr1 = v3
                vr2 = v2
                vr3 = v6
                vr4 = v5
                vl1 = v1
                vl2 = v2
                vl3 = v6
                vl4 = v7
             case(11)                  ! edge 37
                vr1 = v4
                vr2 = v3
                vr3 = v7
                vr4 = v6
                vl1 = v2
                vl2 = v3
                vl3 = v7
                vl4 = v8
             case(12)                  ! edge 48
                vr1 = v1
                vr2 = v4
                vr3 = v8
                vr4 = v7
                vl1 = v3
                vl2 = v4
                vl3 = v8
                vl4 = v5
          end select
   
          do i = 1, 3
             er1(i) = coord(i,vr2,c) - coord(i,vr1,c)
             er2(i) = coord(i,vr3,c) - coord(i,vr2,c)
             er3(i) = coord(i,vr4,c) - coord(i,vr3,c)
             jac_r(1,i) = er1(i)
             jac_r(2,i) = er2(i)
             jac_r(3,i) = er3(i)
             el1(i) = coord(i,vl2,c) - coord(i,vl1,c)
             el2(i) = coord(i,vl3,c) - coord(i,vl2,c)
             el3(i) = coord(i,vl4,c) - coord(i,vl3,c)
             jac_l(1,i) = el1(i)
             jac_l(2,i) = el2(i)
             jac_l(3,i) = el3(i)
          end do

          det_r =   jac_r(1,1)*jac_r(2,2)*jac_r(3,3) &
                  + jac_r(3,1)*jac_r(1,2)*jac_r(2,3) &
                  + jac_r(1,3)*jac_r(3,2)*jac_r(2,1) &
                  - jac_r(1,3)*jac_r(2,2)*jac_r(3,1) &
                  - jac_r(1,1)*jac_r(2,3)*jac_r(3,2) &
                  - jac_r(3,3)*jac_r(1,2)*jac_r(2,1)
   
          det_l =   jac_l(1,1)*jac_l(2,2)*jac_l(3,3) &
                  + jac_l(3,1)*jac_l(1,2)*jac_l(2,3) &
                  + jac_l(1,3)*jac_l(3,2)*jac_l(2,1) &
                  - jac_l(1,3)*jac_l(2,2)*jac_l(3,1) &
                  - jac_l(1,1)*jac_l(2,3)*jac_l(3,2) &
                  - jac_l(3,3)*jac_l(1,2)*jac_l(2,1)
   
          det_l = -det_l     ! Jacobian volume is minus (e1 x e2) . e3
                             ! for base/edge vectors forming left-handed coordinate system

          e1e2e3_r =  SQRT(er1(1)**2+er1(2)**2+er1(3)**2) &
                    * SQRT(er2(1)**2+er2(2)**2+er2(3)**2) &
                    * SQRT(er3(1)**2+er3(2)**2+er3(3)**2)
   
          e1e2e3_l =  SQRT(el1(1)**2+el1(2)**2+el1(3)**2) &
                    * SQRT(el2(1)**2+el2(2)**2+el2(3)**2) &
                    * SQRT(el3(1)**2+el3(2)**2+el3(3)**2)
   
          jac_scaled_r = det_r / (e1e2e3_r + 1.0d-14)
          jac_scaled_l = det_l / (e1e2e3_l + 1.0d-14)
          if (jac_scaled_l <= -1.0d-6 .or. jac_scaled_r <= -1.0d-6) then
             twisted(c,edge) = .true.
          end if
        end do
    end do

    ! see if we have any problems
    if (PGSLIB_GLOBAL_ANY(twisted)) then
       ! we should now map the bad cells back to their original numbers and report

       ! Could print out minimum value of scaled Jacobian
       ! and issue warning if it's smaller than, say, 0.2.
       !
       ! However, for tet meshes, some edges and vertices
       ! are degenerate, in Truchas representation,
       ! hence leading to zero values for jacobian determinant.

       ! instead, we complain and punt with minimal diagnostics
       call TLS_fatal ('TWISTED_CELL_TEST: at least one cell in mesh is tangled/folded')
    end if

    ! Explicitly deallocate temporary array
    deallocate(coord)
    deallocate(twisted)

  end subroutine twisted_cell_test

END MODULE MESH_QUALITY_MODULE
