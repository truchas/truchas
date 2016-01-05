!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE TRIANGLE_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables and procedures necessary to compute interface
  !   triangulated polygons from the interface dump data.
  !
  ! Public Interface(s):
  !
  !   * call INTERFACE_TRIANGLES (nicells, Xv, Rho, Normal)
  !
  ! Contains: INTERFACE_TRIANGLES
  !           POLYGON_VERT
  !           Angl_Permute
  !
  ! Author(s): Douglas B. Kothe (LANL Group MST-8, dbk@lanl.gov)
  !            Matthew Williams (LANL Group MST-8, mww@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  ! Public Procedures
  public :: INTERFACE_TRIANGLES


  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE INTERFACE_TRIANGLES (nicells, Xv, Ro, Normal, N_Order, Nvrt, Pv,  &
                                  Perm, Mdpnt, ortho_mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute triangle vertices associated with each interface plane
    !
    !=======================================================================
    use constants_module,  only: pi
    use parameter_module,  only: ndim, nec, nvc

    ! Arguments
    integer,                              intent(IN)  :: nicells
    real(r8),   dimension(nicells),          intent(IN)  :: Ro
    integer, dimension(nicells),          intent(OUT) :: Nvrt
    real(r8),   dimension(ndim,nicells),     intent(IN)  :: Normal
    real(r8),   dimension(nvc,ndim,nicells), intent(IN)  :: Xv
    integer, dimension(nec,nicells),      intent(OUT) :: N_Order
    integer, dimension(nicells,nec),      intent(OUT) :: Perm
    real(r8),   dimension(ndim,nec,nicells), intent(OUT) :: Pv
    real(r8),   dimension(nicells,ndim),     intent(OUT) :: Mdpnt
    logical,                              intent(IN)  :: ortho_mesh

    ! Local Variables
    integer                              :: edge, i, n
    logical  ,dimension(nec,nicells)     :: Edge_Hit
    integer, dimension(nicells)          :: Cntr
    real(r8),   dimension(nicells,nec,ndim) :: V_vec
    real(r8),   dimension(nicells,nec)      :: Angl, Dist
    real(r8),   dimension(nicells,ndim)     :: Y_vec

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize arrays.
    Edge_Hit = .false.
    Pv       = 0.0_r8
    Nvrt     = 0
    N_Order  = 0

   ! Compute the polygon vertices (arrays Pv and Edge_Hit)
    call POLYGON_VERT (nicells, Xv, Normal, Ro, Edge_Hit, Pv, ortho_mesh)

    ! Compute the number of vertices for each interface.
    do edge = 1, nec
       do i = 1, nicells
          if (Edge_Hit(edge,i)) Nvrt(i) = Nvrt(i) + 1
       end do
    end do

    ! Find which edges are intersected
    ! can be made into subroutine order_vertices

     Cntr = 1
      do edge = 1,nec
      do n=1,nicells
        if(Edge_Hit(edge,n))then
         N_Order(Cntr(n),n) = edge
         Cntr(n) = Cntr(n) + 1
        endif
      enddo
      enddo

    ! Find the midpoint between the first two vertices listed (in Pv)
      do n=1,nicells

      do i=1,ndim
   Mdpnt(n,i) =   Pv(i,N_Order(2,n),n) +          &
                       (Pv(i,N_Order(1,n),n) - Pv(i,N_Order(2,n),n))/2
      enddo
      do i = 1, Nvrt(n)
   V_vec(n,N_Order(i,n),:) = Pv(:,N_Order(i,n),n) - Mdpnt(n,:)
      enddo
      enddo

    ! Normalize V_vec array
   Do n=1,nicells
    Do i = 1,nec
     If(Edge_Hit(i,n))then
       Dist(n,i) = DSqrt(V_vec(n,i,1)**2+V_vec(n,i,2)**2+V_vec(n,i,3)**2)
       V_vec(n,i,1) = V_vec(n,i,1)/Dist(n,i)
       V_vec(n,i,2) = V_vec(n,i,2)/Dist(n,i)
       V_vec(n,i,3) = V_vec(n,i,3)/Dist(n,i)
     end if
    End Do
   End Do

    ! Determine a unit vector (Y_vec) orthogonal to both the normal and V_vec(n,1,:)

   Do n = 1,nicells
    If(V_vec(n,N_Order(1,n),2)*Normal(3,n) - V_vec(n,N_Order(1,n),3)*Normal(2,n) &
           .ne. 0.0_r8)then
      Y_vec(n,1) = 1.0_r8
     If(Normal(3,n)*V_vec(n,N_Order(1,n),2) .ne. 0.0_r8)then
      Y_vec(n,3) = (-Normal(1,n)+(Normal(2,n)*V_vec(n,N_Order(1,n),1))/    &
                     V_vec(n,N_Order(1,n),2))/  &
                   (Normal(3,n)-(Normal(2,n)*V_vec(n,N_Order(1,n),3))/     &
                    V_vec(n,N_Order(1,n),2))
      Y_vec(n,2) = -(V_vec(n,N_Order(1,n),1)+V_vec(n,N_Order(1,n),3)*Y_vec(n,3))/  &
                     V_vec(n,N_Order(1,n),2)
     Else
      Y_vec(n,2) = (-Normal(1,n)+(Normal(3,n)*V_vec(n,N_Order(1,n),1))/    &
                    V_vec(n,N_Order(1,n),3))/  &
                   (Normal(2,n)-(Normal(3,n)*V_vec(n,N_Order(1,n),2))/     &
                    V_vec(n,N_Order(1,n),3))
      Y_vec(n,3) = -(V_vec(n,N_Order(1,n),1)+V_vec(n,N_Order(1,n),2)*Y_vec(n,2))/  &
                     V_vec(n,N_Order(1,n),3)
     End If

    Else
    If(V_vec(n,N_Order(1,n),1)*Normal(3,n) - V_vec(n,N_Order(1,n),3)*Normal(1,n) &
                    .ne. 0.0_r8)then
      Y_vec(n,2) = 1.0_r8
     If(Normal(3,n)*V_vec(n,N_Order(1,n),1) .ne. 0.0_r8)then
      Y_vec(n,3) = (-Normal(2,n)+(Normal(1,n)*V_vec(n,N_Order(1,n),2))/    &
                    V_vec(n,N_Order(1,n),1))/  &
                   (Normal(3,n)-(Normal(1,n)*V_vec(n,N_Order(1,n),3))/     &
                    V_vec(n,N_Order(1,n),1))
      Y_vec(n,1) = -(V_vec(n,N_Order(1,n),2)+V_vec(n,N_Order(1,n),3)*Y_vec(n,3))/   &
                     V_vec(n,N_Order(1,n),1)
     Else
      Y_vec(n,1) = (-Normal(2,n)+(Normal(3,n)*V_vec(n,N_Order(1,n),2))/    &
                    V_vec(n,N_Order(1,n),3))/  &
                   (Normal(1,n)-(Normal(3,n)*V_vec(n,N_Order(1,n),1))/     &
                    V_vec(n,N_Order(1,n),3))
      Y_vec(n,3) = -(V_vec(n,N_Order(1,n),2)+V_vec(n,N_Order(1,n),1)*Y_vec(n,1))/   &
                     V_vec(n,N_Order(1,n),3)
     EndIf

    Else
    If(V_vec(n,N_Order(1,n),1)*Normal(2,n) - V_vec(n,N_Order(1,n),2)*Normal(1,n) &
                   .ne. 0.0_r8)then
      Y_vec(n,3) = 1.0_r8
     If(Normal(2,n)*V_vec(n,N_Order(1,n),1) .ne. 0.0_r8)then
      Y_vec(n,2) = (-Normal(3,n)+(Normal(1,n)*V_vec(n,N_Order(1,n),3))/    &
                    V_vec(n,N_Order(1,n),1))/  &
                   (Normal(2,n)-(Normal(1,n)*V_vec(n,N_Order(1,n),2))/     &
                    V_vec(n,N_Order(1,n),1))
      Y_vec(n,1) = -(V_vec(n,N_Order(1,n),3)+V_vec(n,N_Order(1,n),2)*Y_vec(n,2))/   &
                     V_vec(n,N_Order(1,n),1)
     Else
      Y_vec(n,1) = (-Normal(3,n)+(Normal(2,n)*V_vec(n,N_Order(1,n),3))/    &
                    V_vec(n,N_Order(1,n),2))/  &
                   (Normal(1,n)-(Normal(2,n)*V_vec(n,N_Order(1,n),1))/     &
                    V_vec(n,N_Order(1,n),2))
      Y_vec(n,2) = -(V_vec(n,N_Order(1,n),3)+V_vec(n,N_Order(1,n),1)*Y_vec(n,1))/  &
                     V_vec(n,N_Order(1,n),2)
     Endif

    Endif
    Endif
    End If
   End Do

    ! Normalize Y_vec array
   Do n=1,nicells
       Dist(n,1) = DSqrt(Y_vec(n,1)**2+Y_vec(n,2)**2+Y_vec(n,3)**2)
       Y_vec(n,1) = Y_vec(n,1)/Dist(n,1)
       Y_vec(n,2) = Y_vec(n,2)/Dist(n,1)
       Y_vec(n,3) = Y_vec(n,3)/Dist(n,1)
   End Do


    ! Choose initial vector and find the angle between this vector and others

      Angl = 100.0

   Do n=1,nicells
    Do i = 1, nec

   If(Edge_Hit(i,n))then

  Angl(n,i) =  DOT_PRODUCT(V_vec(n,N_Order(1,n),:),V_vec(n,i,:)) 
  Angl(n,i) =  Angl(n,i)  /  ( &
       DSqrt(V_vec(n,N_Order(1,n),1)**2  + V_vec(n,N_Order(1,n),2)**2 +  &
             V_vec(n,N_Order(1,n),3)**2) *                               &
       DSqrt(V_vec(n,i,1)**2 + V_vec(n,i,2)**2 +   &
              V_vec(n,i,3)**2)) 
   if(ABS(Angl(n,i)+1.0) .lt. 1.0e-08)then
  Angl(n,i) = pi
   else
   if(ABS(Angl(n,i)-1.0) .lt. 1.0e-08)then
  Angl(n,i) = 0.0_r8
   else
  Angl(n,i) = ACOS(Angl(n,i))
   endif
   endif
     If(Y_vec(n,1)*V_vec(n,i,1) + Y_vec(n,2)*V_vec(n,i,2) +   &
             Y_vec(n,3)*V_vec(n,i,3) .lt. 0.0_r8)then
   Angl(n,i) = 2*pi - Angl(n,i)
     End If

   Endif
    Enddo
   Enddo

    ! Find a permutation vector which reorders N_Order
    ! Equivalent to the old Re_Order_Vertices

   Perm = Angl_Permute(Angl,Nvrt,N_Order,nicells,nec)

  END SUBROUTINE INTERFACE_TRIANGLES

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE POLYGON_VERT (nicells, Xv, Normal, Rho, Edge_Hit, Pv, ortho_mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate all the intersection
    !   points of an arbitrary plane and the 12 edges
    !   of an arbitrary hexahedron (potentially nonorthogonal).
    !   The equation describing the plane is assumed to be:
    !             Normal * X = Rho,
    !   where Normal is the normal vector to the plane
    !   and Rho is a constant.  The equation for each edge
    !   of the hexahedron is written in parametric form as:
    !               x = x_0 + a*t
    !               y = y_0 + b*t
    !               z = z_0 + c*t
    !   where (x_0,y_0,z_0) is a point on the line,
    !   (a,b,c) are the components of a tangent to the line,
    !   and t is the parametric variable that is allowed to
    !   vary.
    !=======================================================================
    use parameter_module,  only: ndim, nec, nvc

    ! Arguments
    integer,                              intent(IN)  :: nicells
    real(r8), dimension(nicells),          intent(IN)  :: Rho
    real(r8), dimension(ndim,nicells),     intent(IN)  :: Normal
    real(r8), dimension(nvc,ndim,nicells), intent(IN)  :: Xv
    logical, dimension(nec,nicells),      intent(OUT) :: Edge_Hit
    real(r8),   dimension(ndim,nec,nicells), intent(OUT) :: Pv
    logical,                              intent(IN)  :: ortho_mesh

    ! Local Variables
    integer :: i, edge, n, v1, v2
    real(r8)                          :: epsilon = 1.0e-08
    integer, dimension(2,nec)        :: Edge_Vertex
    real(r8), dimension(ndim,nicells) :: X_0, X_1, Delta
    real(r8), dimension(nec,nicells)  :: T
    real(r8), dimension(nicells)      :: N_Dot_Delta

    ! Set Edge-Vertex connectivity
    data Edge_Vertex /1, 2, 2, 3, 3, 4, 4, 1, 2, 6, 3, 7, &
                      4, 8, 1, 5, 5, 6, 6, 7, 7, 8, 8, 5  /

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Vertex and edge numbering convention:

    !        7*------10-------*6
    !        /|               /| 
    !      11 |              9 |
    !      /  6             /  5
    !    8*------12--------*5  |
    !     |   |            |   |
    !     |  3*-------2----|---*2
    !     7  /             8  /
    !     | 3              | 1
    !     |/               |/
    !     *-------4--------*
    !     4                1

    ! Loop over edges, searching for plane/edge intersections.
    EDGE_LOOP: do edge = 1, nec

       v1 = Edge_Vertex(1,edge)
       v2 = Edge_Vertex(2,edge)

       do n = 1,ndim
          X_0(n,:)   = Xv(v1,n,:)
          X_1(n,:)   = Xv(v2,n,:)
          Delta(n,:) = X_1(n,:) - X_0(n,:)
       end do

       N_Dot_Delta = 0.0_r8
       do n = 1,ndim
          N_Dot_Delta = N_Dot_Delta + Delta(n,:)*Normal(n,:)
       end do

       T(edge,:) = 0.0_r8
       do n = 1,ndim
          T(edge,:) = T(edge,:) + Normal(n,:)*X_0(n,:)
       end do

     where (ABS(N_Dot_Delta) > epsilon)  T(edge,:) = (Rho - T(edge,:))/N_Dot_Delta
     where (ABS(T(edge,:)) .lt. epsilon) T(edge,:) = ABS(T(edge,:))
     where (ABS(1.0_r8 - T(edge,:)) .lt. epsilon) T(edge,:) = 1.0_r8 

     where (ABS(N_Dot_Delta) > epsilon)
       Pv(1,edge,:) = X_0(1,:) + Delta(1,:)*T(edge,:)
       Pv(2,edge,:) = X_0(2,:) + Delta(2,:)*T(edge,:)
       Pv(3,edge,:) = X_0(3,:) + Delta(3,:)*T(edge,:)
       Edge_Hit(edge,:) = T(edge,:) >= 0.0_r8 .and. T(edge,:) <= 1.0_r8
     end where

     if(.not. ortho_mesh)   &
    where(Sqrt(Delta(1,:)**2+Delta(2,:)**2+Delta(3,:)**2) <= epsilon) &
                        Edge_Hit(edge,:) = .false.

    end do EDGE_LOOP

DEGENERATE_CHECK:      if(.not. ortho_mesh)  then

  where( Pv(1,1,:)  .eq. Pv(1,3,:) .and.  &
         Pv(2,1,:)  .eq. Pv(2,3,:) .and.  &
         Pv(3,1,:)  .eq. Pv(3,3,:))  Edge_Hit(3,:) = .false.

  where( Pv(1,1,:)  .eq. Pv(1,9,:) .and.  &
         Pv(2,1,:)  .eq. Pv(2,9,:) .and.  &
         Pv(3,1,:)  .eq. Pv(3,9,:))  Edge_Hit(9,:) = .false.

  where( Pv(1,2,:)  .eq. Pv(1,4,:) .and.  &
         Pv(2,2,:)  .eq. Pv(2,4,:) .and.  &
         Pv(3,2,:)  .eq. Pv(3,4,:))  Edge_Hit(4,:) = .false.

  where( Pv(1,2,:)  .eq. Pv(1,10,:) .and.  &
         Pv(2,2,:)  .eq. Pv(2,10,:) .and.  &
         Pv(3,2,:)  .eq. Pv(3,10,:))  Edge_Hit(10,:) = .false.

  where( Pv(1,3,:)  .eq. Pv(1,11,:) .and.  &
         Pv(2,3,:)  .eq. Pv(2,11,:) .and.  &
         Pv(3,3,:)  .eq. Pv(3,11,:))  Edge_Hit(11,:) = .false.

  where( Pv(1,4,:)  .eq. Pv(1,12,:) .and.  &
         Pv(2,4,:)  .eq. Pv(2,12,:) .and.  &
         Pv(3,4,:)  .eq. Pv(3,12,:))  Edge_Hit(12,:) = .false.

  where( Pv(1,5,:)  .eq. Pv(1,6,:) .and.  &
         Pv(2,5,:)  .eq. Pv(2,6,:) .and.  &
         Pv(3,5,:)  .eq. Pv(3,6,:))  Edge_Hit(6,:) = .false.

  where( Pv(1,5,:)  .eq. Pv(1,8,:) .and.  &
         Pv(2,5,:)  .eq. Pv(2,8,:) .and.  &
         Pv(3,5,:)  .eq. Pv(3,8,:))  Edge_Hit(8,:) = .false.

  where( Pv(1,6,:)  .eq. Pv(1,7,:) .and.  &
         Pv(2,6,:)  .eq. Pv(2,7,:) .and.  &
         Pv(3,6,:)  .eq. Pv(3,7,:))  Edge_Hit(7,:) = .false.

  where( Pv(1,7,:)  .eq. Pv(1,8,:) .and.  &
         Pv(2,7,:)  .eq. Pv(2,8,:) .and.  &
         Pv(3,7,:)  .eq. Pv(3,8,:))  Edge_Hit(8,:) = .false.

  where( Pv(1,9,:)  .eq. Pv(1,11,:) .and.  &
         Pv(2,9,:)  .eq. Pv(2,11,:) .and.  &
         Pv(3,9,:)  .eq. Pv(3,11,:))  Edge_Hit(11,:) = .false.

  where( Pv(1,10,:)  .eq. Pv(1,12,:) .and.  &
         Pv(2,10,:)  .eq. Pv(2,12,:) .and.  &
         Pv(3,10,:)  .eq. Pv(3,12,:))  Edge_Hit(12,:) = .false.

                       end if DEGENERATE_CHECK

  END SUBROUTINE POLYGON_VERT

  Function Angl_Permute (vector,Nvrt,N_Order,nicells,nec)
    !---------------------------------------------------------------------------
    ! purpose:
    !    return a permutation vector that will permute an unsorted
    !    integer vector into a vector sorted into ascending order
    !
    !    vector(permute(i)) is the i'th element of the sorted list
    !
    !    warning: this is a simple n^2 sort, and should not be used
    !    if n is greater than about 10
    !---------------------------------------------------------------------------
    ! arguments
    integer, intent(IN)  :: nicells
    integer, intent(IN)  :: nec
    real(r8),   dimension(nicells, nec),  intent(IN) :: vector
    integer, dimension(nicells),       intent(IN) :: Nvrt
    integer,   dimension(nec,nicells), intent(IN) :: N_Order
    


    ! function return
    integer , Dimension(nicells, nec) :: Angl_Permute

    ! local variables
    integer :: i, j, n, tmp

    !---------------------------------------------------------------------------

    ! initialize permutation vector to do nothing
    Do n = 1,nicells
     Do i = 1,nec
       Angl_Permute(n,i) = i
     End Do
    End Do

    ! bubble sort on the vector, reordering the permutation indices
  Do n = 1,nicells

    Do j = nec-1, 1, -1

       Do i = 1, j
          If (vector(n,Angl_Permute(n,i)) >      &
                vector(n,Angl_Permute(n,i+1))) Then
             ! exchange the permutation indices
             tmp          = Angl_Permute(n,i+1)
             Angl_Permute(n,i+1) = Angl_Permute(n,i)
             Angl_Permute(n,i)   = tmp
          End If
       End Do

    End Do

  End Do

  End Function Angl_Permute


END MODULE TRIANGLE_MODULE
