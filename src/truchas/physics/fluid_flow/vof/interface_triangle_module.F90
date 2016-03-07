!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE INTERFACE_TRIANGLE_MODULE
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
  !           Angle_Permute
  !           INTERFACE_AREA
  !
  ! Author(s): Douglas B. Kothe (LANL Group MST-8, dbk@lanl.gov)
  !            Matthew Williams (LANL Group MST-8, mww@lanl.gov)  
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  public :: INTERFACE_TRIANGLES

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE INTERFACE_TRIANGLES (Xv)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute triangle vertices associated with each interface plane
    !
    !=======================================================================
    use constants_module,     only: pi
    use interface_module,     only: Int_Geom
    use parameter_module,     only: nicells
    use legacy_mesh_api,      only: ndim, nec, nvc

    ! Arguments
    real(r8), dimension(ndim,nvc,nicells), intent(IN) :: Xv

    ! Local Variables
    integer :: edge, i, n
    integer, allocatable, dimension(:,:) :: N_Order, Perm
    integer, allocatable, dimension(:) :: Nvrt, Cntr
    logical, allocatable, dimension(:,:) :: Edge_Hit
    real(r8), allocatable, dimension(:,:,:) :: Pv, V_vec
    real(r8), allocatable, dimension(:,:) :: Mdpnt, Angl, Y_vec, Dist
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Deallocate/Allocate polygon arrays,
    ! Note some arrays can be defined as ragged arrays (nicells,Nvrt)
    ALLOCATE (N_Order(nec,nicells),    &
              Nvrt(nicells),           &
              Pv(ndim,nec,nicells),    &
              Edge_Hit(nec,nicells),   &
              Mdpnt(nicells,ndim),     &
              Cntr(nicells),           &
              V_vec(nicells,nec,ndim), &
              Angl(nicells,nec),       &
              Perm(nicells,nec),       &
              Y_vec(nicells,ndim),     &
              Dist(nicells,nec))

    ! Initialize arrays.
    Edge_Hit = .false.
    Pv = 0.0_r8
    Nvrt = 0
    N_Order = 0
    V_vec = 0.0_r8 

    ! Compute the polygon vertices (arrays Pv and Edge_Hit)
    call POLYGON_VERT (Xv, Pv, Edge_Hit)

    ! Compute the number of vertices for each interface.
    do edge = 1,nec
       do i = 1,nicells
          if (Edge_Hit(edge,i)) Nvrt(i) = Nvrt(i) + 1
       end do
    end do

    ! Find which edges are intersected
    Cntr = 1
    do edge = 1,nec
       do n = 1,nicells
          if (Edge_Hit(edge,n)) then
             N_Order(Cntr(n),n) = edge
             Cntr(n) = Cntr(n) + 1
          end if
       end do
    end do

    ! Find the midpoint between the first two vertices listed (in Pv)
    do n = 1,nicells
       do i = 1,ndim
          Mdpnt(n,i) = Pv(i,N_Order(2,n),n) + (Pv(i,N_Order(1,n),n) - Pv(i,N_Order(2,n),n))/2
       end do
       do i = 1,Nvrt(n)
          V_vec(n,N_Order(i,n),:) = Pv(:,N_Order(i,n),n) - Mdpnt(n,:)
       end do
    end do

    ! Determine a unit vector (Y_vec) orthogonal to both the normal and V_vec(n,1,:)
    do n = 1,nicells

       if (  V_vec(n,N_Order(1,n),2)*Int_Geom(n)%Normal(3) &
           - V_vec(n,N_Order(1,n),3)*Int_Geom(n)%Normal(2) .ne. 0.0_r8) then
          Y_vec(n,1) = 1.0_r8
          if (Int_Geom(n)%Normal(3)*V_vec(n,N_Order(1,n),2) .ne. 0.0_r8) then
             Y_vec(n,3) = (-Int_Geom(n)%Normal(1) + Int_Geom(n)%Normal(2)*V_vec(n,N_Order(1,n),1)/V_vec(n,N_Order(1,n),2)) &
                        / ( Int_Geom(n)%Normal(3) - Int_Geom(n)%Normal(2)*V_vec(n,N_Order(1,n),3)/V_vec(n,N_Order(1,n),2))
             Y_vec(n,2) = -(V_vec(n,N_Order(1,n),1) + V_vec(n,N_Order(1,n),3)*Y_vec(n,3)) / V_vec(n,N_Order(1,n),2)
          else
             Y_vec(n,2) = (-Int_Geom(n)%Normal(1) + Int_Geom(n)%Normal(3)*V_vec(n,N_Order(1,n),1)/V_vec(n,N_Order(1,n),3)) &
                        / ( Int_Geom(n)%Normal(2) - Int_Geom(n)%Normal(3)*V_vec(n,N_Order(1,n),2)/V_vec(n,N_Order(1,n),3))
             Y_vec(n,3) = -(V_vec(n,N_Order(1,n),1) + V_vec(n,N_Order(1,n),2)*Y_vec(n,2)) / V_vec(n,N_Order(1,n),3)
          end if

       else if ( V_vec(n,N_Order(1,n),1)*Int_Geom(n)%Normal(3) &
               - V_vec(n,N_Order(1,n),3)*Int_Geom(n)%Normal(1) .ne. 0.0_r8) then
          Y_vec(n,2) = 1.0_r8
          if (Int_Geom(n)%Normal(3)*V_vec(n,N_Order(1,n),1) .ne. 0.0_r8) then
             Y_vec(n,3) = (-Int_Geom(n)%Normal(2) + Int_Geom(n)%Normal(1)*V_vec(n,N_Order(1,n),2)/V_vec(n,N_Order(1,n),1)) &
                        / ( Int_Geom(n)%Normal(3) - Int_Geom(n)%Normal(1)*V_vec(n,N_Order(1,n),3)/V_vec(n,N_Order(1,n),1))
             Y_vec(n,1) = -(V_vec(n,N_Order(1,n),2) + V_vec(n,N_Order(1,n),3)*Y_vec(n,3)) / V_vec(n,N_Order(1,n),1)
          else
             Y_vec(n,1) = (-Int_Geom(n)%Normal(2) + Int_Geom(n)%Normal(3)*V_vec(n,N_Order(1,n),2)/V_vec(n,N_Order(1,n),3)) &
                        / ( Int_Geom(n)%Normal(1) - Int_Geom(n)%Normal(3)*V_vec(n,N_Order(1,n),1)/V_vec(n,N_Order(1,n),3))
             Y_vec(n,3) = -(V_vec(n,N_Order(1,n),2) + V_vec(n,N_Order(1,n),1)*Y_vec(n,1)) / V_vec(n,N_Order(1,n),3)
          end if

       else if ( V_vec(n,N_Order(1,n),1)*Int_Geom(n)%Normal(2) &
               - V_vec(n,N_Order(1,n),2)*Int_Geom(n)%Normal(1) .ne. 0.0_r8) then
          Y_vec(n,3) = 1.0_r8
          if (Int_Geom(n)%Normal(2)*V_vec(n,N_Order(1,n),1) .ne. 0.0_r8) then
             Y_vec(n,2) = (-Int_Geom(n)%Normal(3) + Int_Geom(n)%Normal(1)*V_vec(n,N_Order(1,n),3)/V_vec(n,N_Order(1,n),1)) &
                        / ( Int_Geom(n)%Normal(2) - Int_Geom(n)%Normal(1)*V_vec(n,N_Order(1,n),2)/V_vec(n,N_Order(1,n),1))
             Y_vec(n,1) = -(V_vec(n,N_Order(1,n),3) + V_vec(n,N_Order(1,n),2)*Y_vec(n,2)) / V_vec(n,N_Order(1,n),1)
          else
             Y_vec(n,1) = (-Int_Geom(n)%Normal(3) + Int_Geom(n)%Normal(2)*V_vec(n,N_Order(1,n),3)/V_vec(n,N_Order(1,n),2)) &
                        / ( Int_Geom(n)%Normal(1) - Int_Geom(n)%Normal(2)*V_vec(n,N_Order(1,n),1)/V_vec(n,N_Order(1,n),2))
             Y_vec(n,2) = -(V_vec(n,N_Order(1,n),3) + V_vec(n,N_Order(1,n),1)*Y_vec(n,1)) / V_vec(n,N_Order(1,n),2)
          end if

       end if
    end do

    ! Normalize Y_vec array
    do n = 1,nicells
       Dist(n,1) = DSqrt(Y_vec(n,1)**2+Y_vec(n,2)**2+Y_vec(n,3)**2)
       Y_vec(n,1) = Y_vec(n,1)/Dist(n,1)
       Y_vec(n,2) = Y_vec(n,2)/Dist(n,1)
       Y_vec(n,3) = Y_vec(n,3)/Dist(n,1)
    end do

    ! Choose initial vector and find the angle between this vector and others
    Angl = 100.d0

    do n = 1,nicells
       do i = 1,nec
          if (Edge_Hit(i,n)) then
             Angl(n,i) =  DOT_PRODUCT(V_vec(n,N_Order(1,n),:),V_vec(n,i,:)) 
             Angl(n,i) =  Angl(n,i) / ( DSqrt(V_vec(n,N_Order(1,n),1)**2 + V_vec(n,N_Order(1,n),2)**2 + V_vec(n,N_Order(1,n),3)**2) &
                                      * DSqrt(V_vec(n,i,1)**2 + V_vec(n,i,2)**2 + V_vec(n,i,3)**2)) 
             if (ABS(Angl(n,i)+1.0) .lt. 1.0e-08) then
                Angl(n,i) = pi
             else if (ABS(Angl(n,i)-1.0) .lt. 1.0e-08) then
                Angl(n,i) = 0.0_r8
             else
                Angl(n,i) = ACOS(Angl(n,i))
             end if
          end if
          if (Y_vec(n,1)*V_vec(n,i,1) + Y_vec(n,2)*V_vec(n,i,2) + Y_vec(n,3)*V_vec(n,i,3) .lt. 0.0_r8) then
             Angl(n,i) = 2*pi - Angl(n,i)
          end if
       end do
    end do

    ! Find a permutation vector which reorders N_Order
    Perm = Angle_Permute (Angl, Nvrt, N_Order, nicells, nec)

    call INTERFACE_AREA (V_Vec, Nvrt, Perm)
     
    DEALLOCATE (N_Order,  &
                Pv,       &
                Nvrt,     &
                Edge_Hit, &
                Mdpnt,    &
                Cntr,     &
                V_vec,    &
                Angl,     &
                Perm,     &
                Y_vec,    &
                Dist)
    
  END SUBROUTINE INTERFACE_TRIANGLES

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE POLYGON_VERT (Xv, Pv, Edge_Hit)
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
    use interface_module,  only: Int_Geom
    use parameter_module,  only: nicells
    use legacy_mesh_api,   only: orthogonal_mesh, ndim, nec, nvc

    ! Arguments
    real(r8), dimension(ndim,nvc,nicells), intent(IN) :: Xv
    real(r8), dimension(ndim, nec, nicells), intent(INOUT) :: Pv
    logical, dimension(nec,nicells), intent(INOUT) :: Edge_Hit

    ! Local Variables
    integer :: edge, n, v1, v2
    integer, dimension(2,nec) :: Edge_Vertex
    real(r8) :: epsilon = 1.0e-07
    real(r8), dimension(ndim,nicells) :: X_0, X_1, Delta
    real(r8), dimension(nec,nicells) :: T
    real(r8), dimension(nicells) :: N_Dot_Delta

    ! Set Edge-Vertex connectivity
    data Edge_Vertex /1, 2, 2, 3, 3, 4, 4, 1, 2, 6, 3, 7, &
                      4, 8, 1, 5, 5, 6, 6, 7, 7, 8, 8, 5  /

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Vertex and edge numbering convention:

    !        7*------10-------*6
    !        /|               /| 
    !      11 |              9 |     Reference for orthogonal coords:
    !      /  6             /  5
    !    8*------12--------*5  |
    !     |   |            |   |        zeta ^ eta
    !     |  3*-------2----|---*2            | /
    !     7  /             8  /              |/
    !     | 3              | 1               *--> xi
    !     |/               |/
    !     *-------4--------*
    !     4                1

    ! Loop over edges, searching for plane/edge intersections.
    EDGE_LOOP: do edge = 1,nec

       v1 = Edge_Vertex(1,edge)
       v2 = Edge_Vertex(2,edge)

       do n = 1,ndim
          X_0(n,:) = Xv(n,v1,:)
          X_1(n,:) = Xv(n,v2,:)
          Delta(n,:) = X_1(n,:) - X_0(n,:)
       end do

       N_Dot_Delta = 0.0_r8
       do n = 1,ndim
          N_Dot_Delta = N_Dot_Delta + Delta(n,:)*Int_Geom(:)%Normal(n)
       end do

       T(edge,:) = 0.0_r8
       do n = 1,ndim
          T(edge,:) = T(edge,:) + Int_Geom(:)%Normal(n)*X_0(n,:)
       end do

       where (ABS(N_Dot_Delta) > epsilon)  T(edge,:) = (Int_Geom(:)%Rho - T(edge,:))/N_Dot_Delta
       where (ABS(T(edge,:)) < epsilon) T(edge,:) = ABS(T(edge,:))
       where (ABS(1.0_r8 - T(edge,:)) < epsilon) T(edge,:) = 1.0_r8 

       where (ABS(N_Dot_Delta) > epsilon)
          Pv(1,edge,:) = X_0(1,:) + Delta(1,:)*T(edge,:)
          Pv(2,edge,:) = X_0(2,:) + Delta(2,:)*T(edge,:)
          Pv(3,edge,:) = X_0(3,:) + Delta(3,:)*T(edge,:)
          Edge_Hit(edge,:) = T(edge,:) >= 0.0_r8 .and. T(edge,:) <= 1.0_r8
       end where

       if (.not. orthogonal_mesh) &
          where (Sqrt(Delta(1,:)**2 + Delta(2,:)**2 + Delta(3,:)**2) <= epsilon) &
             Edge_Hit(edge,:) = .false.

    end do EDGE_LOOP

    DEGENERATE_CHECK: if (.not. orthogonal_mesh) then

       where (Pv(1,1,:) .eq. Pv(1,3,:) .and. &
              Pv(2,1,:) .eq. Pv(2,3,:) .and. &
              Pv(3,1,:) .eq. Pv(3,3,:))      &
          Edge_Hit(3,:) = .false.

       where (Pv(1,1,:) .eq. Pv(1,9,:) .and. &
              Pv(2,1,:) .eq. Pv(2,9,:) .and. &
              Pv(3,1,:) .eq. Pv(3,9,:))      &
          Edge_Hit(9,:) = .false.

       where (Pv(1,2,:) .eq. Pv(1,4,:) .and. &
              Pv(2,2,:) .eq. Pv(2,4,:) .and. &
              Pv(3,2,:) .eq. Pv(3,4,:))      &
          Edge_Hit(4,:) = .false.

       where (Pv(1,2,:) .eq. Pv(1,10,:) .and. &
              Pv(2,2,:) .eq. Pv(2,10,:) .and. &
              Pv(3,2,:) .eq. Pv(3,10,:))      &
          Edge_Hit(10,:) = .false.
     
       where (Pv(1,3,:) .eq. Pv(1,11,:) .and. &
              Pv(2,3,:) .eq. Pv(2,11,:) .and. &
              Pv(3,3,:) .eq. Pv(3,11,:))      &
          Edge_Hit(11,:) = .false.
     
       where (Pv(1,4,:) .eq. Pv(1,12,:) .and. &
              Pv(2,4,:) .eq. Pv(2,12,:) .and. &
              Pv(3,4,:) .eq. Pv(3,12,:))      &
          Edge_Hit(12,:) = .false.
     
       where (Pv(1,5,:) .eq. Pv(1,6,:) .and. &
              Pv(2,5,:) .eq. Pv(2,6,:) .and. &
              Pv(3,5,:) .eq. Pv(3,6,:))      &
          Edge_Hit(6,:) = .false.
     
       where (Pv(1,5,:) .eq. Pv(1,8,:) .and. &
              Pv(2,5,:) .eq. Pv(2,8,:) .and. &
              Pv(3,5,:) .eq. Pv(3,8,:))      &
          Edge_Hit(8,:) = .false.
     
       where (Pv(1,6,:) .eq. Pv(1,7,:) .and. &
              Pv(2,6,:) .eq. Pv(2,7,:) .and. &
              Pv(3,6,:) .eq. Pv(3,7,:))      &
          Edge_Hit(7,:) = .false.
     
       where (Pv(1,7,:) .eq. Pv(1,8,:) .and. &
              Pv(2,7,:) .eq. Pv(2,8,:) .and. &
              Pv(3,7,:) .eq. Pv(3,8,:))      &
          Edge_Hit(8,:) = .false.
     
       where (Pv(1,9,:) .eq. Pv(1,11,:) .and. &
              Pv(2,9,:) .eq. Pv(2,11,:) .and. &
              Pv(3,9,:) .eq. Pv(3,11,:))      &
          Edge_Hit(11,:) = .false.
     
       where (Pv(1,10,:) .eq. Pv(1,12,:) .and. &
              Pv(2,10,:) .eq. Pv(2,12,:) .and. &
              Pv(3,10,:) .eq. Pv(3,12,:))      &
          Edge_Hit(12,:) = .false.

    end if DEGENERATE_CHECK

  END SUBROUTINE POLYGON_VERT

  Function Angle_Permute (vector, Nvrt, N_Order, nicells, nec)
    !=======================================================================
    ! purpose:
    !    return a permutation vector that will permute an unsorted
    !    integer vector into a vector sorted into ascending order,
    !    yet leave the first cell untouched (reorder for two down).
    !
    !    vector(permute(i)) is the i'th element of the sorted list
    !
    !    warning: this is a simple n^2 sort, and should not be used
    !    if n is greater than about 10
    !=======================================================================

    ! arguments
    integer, intent(IN) :: nicells
    integer, intent(IN) :: nec
    real(r8), dimension(nicells,nec), intent(IN) :: vector
    integer, dimension(nicells), intent(IN) :: Nvrt
    integer, dimension(nec,nicells), intent(IN) :: N_Order
    
    ! function return
    integer, dimension(nicells,nec) :: Angle_Permute

    ! local variables
    integer :: i, j, n, tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! initialize permutation vector to do nothing
    do n = 1,nicells
       do i = 1,nec
          Angle_Permute(n,i) = i
       end do
    end do

    ! bubble sort on the vector, reordering the permutation indices
    do n = 1,nicells
       do j = nec-1,1,-1
          do i = 1,j
             if (vector(n,Angle_Permute(n,i)) > vector(n,Angle_Permute(n,i+1))) then
             ! exchange the permutation indices
                tmp = Angle_Permute(n,i+1)
                Angle_Permute(n,i+1) = Angle_Permute(n,i)
                Angle_Permute(n,i)   = tmp
             end if
          end do
       end do
    end do

  End Function Angle_Permute

  SUBROUTINE INTERFACE_AREA (V_vec, Nvrt, Perm)
    !=======================================================================
    ! Purpose(s):
    !
    ! Calculate area of interface in each cell containing an interface.
    !
    !=======================================================================
    use interface_module,  only: Int_Geom
    use parameter_module,  only: nicells
    use legacy_mesh_api,   only: ndim, nec

    ! Arguments
    real(r8), dimension(nicells,nec,ndim), intent(IN) :: V_vec
    integer, dimension(nicells), intent(IN) :: Nvrt
    integer, dimension(nicells,nec), intent(IN) :: Perm

    ! Local Variables
    integer :: i, n, incr
    real(r8) :: x_1, x_2, y_1, y_2, z_1, z_2
    real(r8), dimension(nicells,nec) :: Pv_1, Pv_2, Pv_3

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! account for possible 2 dimension problem
    select case (ndim)
       case (2)
          Pv_1 = V_vec(:,:,1)
          Pv_2 = V_vec(:,:,2)
       case (3)
          Pv_1 = V_vec(:,:,1)
          Pv_2 = V_vec(:,:,2)
          Pv_3 = V_vec(:,:,3)
    end select

    Int_Geom%Area = 0.0_r8    

    do n = 1,nicells
       incr = 0

       do i = 2,Nvrt(n)
          incr = incr + 1
          select case (ndim)
             case (2)
                x_1 = Pv_1(n,Perm(n,incr))
                x_2 = Pv_1(n,Perm(n,incr+1))
                y_1 = Pv_2(n,Perm(n,incr))
                y_2 = Pv_2(n,Perm(n,incr+1))
                Int_Geom(n)%Area = Int_Geom(n)%Area &
                                 + sqrt((x_1*y_2-y_1*x_2)**2)
             case (3)
                x_1 = Pv_1(n,Perm(n,incr))
                x_2 = Pv_1(n,Perm(n,incr+1))
                y_1 = Pv_2(n,Perm(n,incr))
                y_2 = Pv_2(n,Perm(n,incr+1))
                z_1 = Pv_3(n,Perm(n,incr))
                z_2 = Pv_3(n,Perm(n,incr+1))
                Int_Geom(n)%Area = Int_Geom(n)%Area &
                                 + sqrt((y_1*z_2-z_1*y_2)**2 &
                                      + (x_1*z_2-z_1*x_2)**2 &
                                      + (x_1*y_2-y_1*x_2)**2)
          end select
       end do

       x_1 = Pv_1(n,Perm(n,Nvrt(n)))
       x_2 = Pv_1(n,Perm(n,1))
       y_1 = Pv_2(n,Perm(n,Nvrt(n)))
       y_2 = Pv_2(n,Perm(n,1))
       z_1 = Pv_3(n,Perm(n,Nvrt(n)))
       z_2 = Pv_3(n,Perm(n,1))
       Int_Geom(n)%Area = Int_Geom(n)%Area + &
          sqrt((y_1*z_2-z_1*y_2)**2 + (x_1*z_2-z_1*x_2)**2 + (x_1*y_2-y_1*x_2)**2)
    end do 

    Int_Geom%Area = 0.5_r8 * Int_Geom%Area

  END SUBROUTINE INTERFACE_AREA

END MODULE INTERFACE_TRIANGLE_MODULE
