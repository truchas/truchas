!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DISCRETE_OP_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define various procedures to perform discrete operations
  !   such as gradient, divergence, etc.
  !
  !   Public Interface(s):
  !     * FACE_CENTER (Cell, Mesh, Vertex, Q_Cell, Q_Face)
  !
  !         Evaluate a cell-face centroid quantity given a cell-center quantity.
  !
  !     * GRADIENT (dPhi_dx, dPhi_dy, dPhi_dz, Phi, Phi_Bit, Phi_BC, method)
  !
  !         Returns the cell-centered gradient (dPhi_dx,dPhi_dy,dPhi_dz)
  !         of a cell-centered quantity Phi.
  !
  !     * VERTEX_AVG (X_cell, X_vertex)
  ! 
  !         Returns the vertex avg (X_vertex) of a cell-centered quantity
  !         (X_cell).
  !
  ! Contains: DETERMINANT_VOL_AVG
  !           FACE_CENTER
  !           GRADIENT
  !           LSLR
  !           VERTEX_AVG
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !
  !=======================================================================
! Using gs_module here, instead of in every routine, reduces the
! virtual size required to compile (discrete_operators) file under linux/fujitsu from
! 560+ mb to ~70 mb, and reduces the time required greatly.  Seems we
! must do this to get around a hard kernel limit on rockhopper.

  use kinds, only: r8
  use gs_module, only: EN_GATHER, EE_GATHER, EN_SUM_SCATTER
  implicit none
  private

  ! Public Subroutines
  public :: DETERMINANT_VOL_AVG, &
            FACE_CENTER,         &
            GRADIENT,            &
            GRADIENT_CELL,       &
            VERTEX_AVG,          &
            DYNAMIC_PRESSURE_FACE_GRADIENT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE DETERMINANT_VOL_AVG (Xv, Yv, Zv, Avg)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the volume-averaged value of a Jacobian determinant |J|:
    !   Avg = Integral(|J|dV), where J is given by column vectors X, Y,
    !   and Z:             -                       -
    !                      |  X_xi   Y_xi   Z_xi   |
    !                J =   |  X_eta  Y_eta  Z_eta  |
    !                      |  X_zeta Y_zeta Z_zeta |
    !                      -                       -
    !   where X_xi == dX/dxi, X_eta == dX/deta, X_zeta == dX/dzeta, etc.
    !   The volume integral is converted to a surface integral and
    !   evaluated following J. Dukowicz, JCP 74: 493-496 (1988).
    !
    !=======================================================================
    use parameter_module, only: ncells, nfc, nvc

    ! Arguments
    real(r8), dimension(nvc,ncells), intent(IN)  :: Xv, Yv, Zv
    real(r8), dimension(ncells),     intent(OUT) :: Avg

    ! Local Variables
    integer :: f
    integer :: v1, v2, v3, v4, v5, v6

    real(r8), dimension(ncells) :: X1, Y1, Z1
    real(r8), dimension(ncells) :: X2, Y2, Z2
    real(r8), dimension(ncells) :: X3, Y3, Z3

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant quantities
    Avg = 0.0_r8

    ! Loop over faces, accumulating the average
    do f = 1,nfc
       select case (f)
          case (1)
          ! Side 1 (vertices 4-8-7-3)
          v1 = 8; v2 = 4; v3 = 7; v4 = 8; v5 = 3; v6 = 4
          case (2)
          ! Side 2 (vertices 1-2-6-5)
          v1 = 6; v2 = 2; v3 = 5; v4 = 6; v5 = 1; v6 = 2
          case (3)
          ! Side 3 (vertices 4-1-5-8)
          v1 = 5; v2 = 1; v3 = 8; v4 = 5; v5 = 4; v6 = 1
          case (4)
          ! Side 4 (vertices 3-7-6-2)
          v1 = 7; v2 = 3; v3 = 6; v4 = 7; v5 = 2; v6 = 3
          case (5)
          ! Side 5 (vertices 4-3-2-1)
          v1 = 3; v2 = 4; v3 = 2; v4 = 3; v5 = 1; v6 = 4
          case (6)
          ! Side 6 (vertices 8-5-6-7)
          v1 = 6; v2 = 5; v3 = 7; v4 = 6; v5 = 8; v6 = 5
       end select

       X1 = Xv(v1,:) + Xv(v2,:)
       Y1 = Yv(v1,:) + Yv(v2,:)
       Z1 = Zv(v1,:) + Zv(v2,:)

       X2 = Xv(v3,:) + Xv(v4,:)
       Y2 = Yv(v3,:) + Yv(v4,:)
       Z2 = Zv(v3,:) + Zv(v4,:)

       X3 = Xv(v5,:) + Xv(v6,:)
       Y3 = Yv(v5,:) + Yv(v6,:)
       Z3 = Zv(v5,:) + Zv(v6,:)

       Avg = Avg + X1*(Y2*Z3 - Y3*Z2) + Y1*(X3*Z2 - X2*Z3) + Z1*(X2*Y3 - X3*Y2)
    end do

    Avg = Avg / 12.0_r8

  END SUBROUTINE DETERMINANT_VOL_AVG

  SUBROUTINE FACE_CENTER (Q_Cell, Q_Face)
    !=======================================================================
    ! Purpose(s):
    !   Evaluate a cell-face centroid quantity given a cell-center quantity.
    !=======================================================================
    use linear_module,      only: LINEAR_PROP
    use parameter_module,   only: ncells, nfc, nnodes, nvc

    ! Argument List

    real(r8), dimension(ncells),     intent(IN)  :: Q_Cell
    real(r8), dimension(nfc,ncells), intent(OUT) :: Q_Face

    ! Local Variables
    integer :: f
    real(r8), dimension(nnodes)     :: Q_Vrtx
    real(r8), dimension(nvc,ncells) :: Q_Vrtx_Cell

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Cell-Vertex Quantity ...................................
    ! Vertex Average
    call VERTEX_AVG (Q_Cell, Q_Vrtx)

    ! Gather Quantity
    call EN_GATHER  (Q_Vrtx_Cell, Q_Vrtx)

    ! Cell-Face Quantity .....................................
    ! Initialize Quantity
    Q_Face = 0.0_r8

    ! Evaluate Quantity
    do f = 1, nfc
       call LINEAR_PROP (f, Q_Vrtx_Cell, Q_Face(f,:))
    end do

  END SUBROUTINE FACE_CENTER

  SUBROUTINE GRADIENT (dPhi_dx, dPhi_dy, dPhi_dz, Phi, &
                       Weight_Data, Phi_Bit, Phi_BC, method)
    !=======================================================================
    ! Purpose(s):
    ! 
    !   Compute the cell-centered gradient (dPhi_dx,dPhi_dy,dPhi_dz) of a
    !   cell-centered scalar quantity Phi with one of the following methods:
    !   'volume-average':
    !       By averaging the gradient, defined by the solution
    !       of a 3x3 system of equations, over the volume.  This
    !       averaging is accomplished with an integral over the
    !       cell volume.  The gradient at any logical point within
    !       the cell is a ratio of Jacobian-like determinants.
    !   'green-gauss':
    !       With a discrete approximation to Gauss's theorem, whereby
    !       the integral of (dPhi_dx,dPhi_dy,dPhi_dz) over the cell
    !       volume is converted to an integral of Phi over the cell
    !       surface area vector.  The area integral is approximated
    !       as a discrete sum over cell faces.  This control volume
    !       formulation is discretely conservative, i.e., adjacent
    !       face contributions will telescope, leaving only boundary
    !       contributions.
    !   'lslr':
    !       Least-squares linear reconstruction.  This technique has
    !       been popularized by Tim Barth and coworkers at NASA/Ames
    !       (e.g., see AIAA-93-0688).  A Taylor series expansion is
    !       constructed from each cell to each of its neighbors and
    !       the difference between the expansion and the actual
    !       neighboring value is minimized in the least-squares sense.
    !
    !=======================================================================
    use ArrayAllocate_Module, only: ArrayCreate, ArrayDestroy
    use bc_module,            only: Vel
    use cutoffs_module,       only: alittle
    use linear_module,        only: LINEAR_PROP
    use mesh_module,          only: Cell, Vertex, Vrtx_Bdy
    use parameter_module,     only: ncells, nfc, nnodes, nvc

    ! Arguments
    character(LEN = *), intent(IN), optional :: method
    real(r8), dimension(ncells), intent(OUT) :: dPhi_dx
    real(r8), dimension(ncells), intent(OUT) :: dPhi_dy
    real(r8), dimension(ncells), intent(OUT) :: dPhi_dz
    real(r8), dimension(ncells), intent(IN)  :: Phi
    real(r8), dimension(ncells), intent(IN), optional :: Weight_Data
    real(r8), dimension(nfc,ncells), intent(IN), optional :: Phi_BC
    integer,  dimension(nfc), intent(IN), optional :: Phi_Bit

    ! Local Variables
    character(80) :: gradient_method = 'volume-average'
    logical :: specified_method
    logical :: dirichlet_bc, velocity_bc
    integer :: f

    real(r8), dimension(ncells) :: Phi_f
    real(r8), dimension(:,:), pointer :: Phi_e
    real(r8), dimension(nvc,ncells) :: Phi_v
    real(r8), dimension(nnodes) :: Phi_vtx

    real(r8), dimension(:,:), pointer :: X_v, Y_v, Z_v
    real(r8), dimension(:,:), pointer :: X_e, Y_e, Z_e

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! If the method has been specified, take it
    specified_method = PRESENT(method)
    if (specified_method) gradient_method = method

    ! See if Dirichlet BCs must be taken into account
    dirichlet_bc = PRESENT(Phi_Bit) .and. PRESENT(Phi_BC)

    ! If Dirichlet BCs, see if they're for velocity (inflow)
    ! (just check for equivalency in first face bit location)
    if (dirichlet_bc) velocity_bc = Vel%Face_Bit(1) - Phi_Bit(1) == 0

    ! Compute Phi_vtx, a vertex value of Phi
    call VERTEX_AVG (Phi, Phi_vtx)

    ! Gather vertex Phi_vtx into cell vector Phi_v
    call EN_GATHER (Phi_v, Phi_vtx)

!!$    ! Correct the vertex averages on faces having a Dirichlet BC
!!$    if (dirichlet_bc) then
!!$       do f = 1,nfc
!!$          if (velocity_bc) then
!!$             Mask = In_Flow (BC%Flag, Phi_Bit(f))
!!$          else
!!$             Mask = Dirichlet (BC%Flag, Phi_Bit(f))
!!$          end if
!!$          do v = 1,nvf
!!$             where (Mask) Phi_v(Face_Vrtx(f,v),:) = Phi_BC(f,:)
!!$          end do
!!$       end do
!!$    end if

    ! Allocate arrays for those cases that need them
    select case (gradient_method)
       case default

       case ('volume-average', 'lslr')
          ! Allocate coordinate arrays
          call ARRAYCREATE (X_v, 1, nvc, 1, ncells, 'Array X_v(nvc,ncells)')
          call ARRAYCREATE (Y_v, 1, nvc, 1, ncells, 'Array Y_v(nvc,ncells)')
          call ARRAYCREATE (Z_v, 1, nvc, 1, ncells, 'Array Z_v(nvc,ncells)')

          ! Gather the vertex coordinates into (X_v,Y_v,Z_v)
          call EN_GATHER (X_v, Vertex%Coord(1), BOUNDARY=Vrtx_Bdy(1)%Data)
          call EN_GATHER (Y_v, Vertex%Coord(2), BOUNDARY=Vrtx_Bdy(2)%Data)
          call EN_GATHER (Z_v, Vertex%Coord(3), BOUNDARY=Vrtx_Bdy(3)%Data)
    end select

    ! Choose method by which gradient will be computed
    select case (gradient_method)
       case ('volume-average')
          ! Compute the volume-averaged gradient
          call DETERMINANT_VOL_AVG (Phi_v, Y_v, Z_v, dPhi_dx)
          call DETERMINANT_VOL_AVG (X_v, Phi_v, Z_v, dPhi_dy)
          call DETERMINANT_VOL_AVG (X_v, Y_v, Phi_v, dPhi_dz)

          ! Normalize by cell volume
          Phi_f = 1.0_r8/Cell%Volume
          dPhi_dx = dPhi_dx*Phi_f
          dPhi_dy = dPhi_dy*Phi_f
          dPhi_dz = dPhi_dz*Phi_f

       case ('Green-Gauss')
          ! Loop over faces, accumulating the product
          ! Phi_f*Face_Normal for each area vector component
          dPhi_dx = 0.0_r8; dPhi_dy = 0.0_r8; dPhi_dz = 0.0_r8
          do f = 1,nfc
             ! Interpolate vertex values to this face
             call LINEAR_PROP (f, Phi_v, Phi_f)

             ! Accumulate the dot product
             dPhi_dx = dPhi_dx + Cell%Face_Normal(1,f)*Cell%Face_Area(f)*Phi_f
             dPhi_dy = dPhi_dy + Cell%Face_Normal(2,f)*Cell%Face_Area(f)*Phi_f
             dPhi_dz = dPhi_dz + Cell%Face_Normal(3,f)*Cell%Face_Area(f)*Phi_f
          end do

          ! Normalize by cell volume
          Phi_f = 1.0_r8/Cell%Volume
          dPhi_dx = dPhi_dx*Phi_f
          dPhi_dy = dPhi_dy*Phi_f
          dPhi_dz = dPhi_dz*Phi_f

       case ('lslr')
          ! Allocate arrays for neighboring coordinates and Phi values
          call ARRAYCREATE (X_e, 1, nfc, 1, ncells, 'Array X_e(nfc,ncells)')
          call ARRAYCREATE (Y_e, 1, nfc, 1, ncells, 'Array Y_e(nfc,ncells)')
          call ARRAYCREATE (Z_e, 1, nfc, 1, ncells, 'Array Z_e(nfc,ncells)')
          call ARRAYCREATE (Phi_e, 1, nfc, 1, ncells, 'Array Phi_e(nfc,ncells)')

          ! Gather neighboring coordinates and Phi values
          call EE_GATHER (X_e, Cell%Centroid(1))
          call EE_GATHER (Y_e, Cell%Centroid(2))
          call EE_GATHER (Z_e, Cell%Centroid(3))
          call EE_GATHER (Phi_e, Phi)

!!$          ! Compute the least-squares linear reconstruction gradient
!!$          if (dirichlet_bc) then
!!$             ! Correct neighbor coordinates and Phi
!!$             ! values on faces having a Dirichlet BC
!!$             do f = 1,nfc
!!$                ! Set the mask for this face
!!$                if (velocity_bc) then
!!$                   Mask = In_Flow (BC%Flag, Phi_Bit(f))
!!$                else
!!$                   Mask = Dirichlet (BC%Flag, Phi_Bit(f))
!!$                end if
!!$                if (ANY(Mask)) then
!!$                   ! Get the face coordinates
!!$                   call LINEAR_PROP (f, X_v, dPhi_dx)
!!$                   call LINEAR_PROP (f, Y_v, dPhi_dy)
!!$                   call LINEAR_PROP (f, Z_v, dPhi_dz)
!!$                   where (Mask)
!!$                      Phi_e(f,:) = Phi_BC(f,:)
!!$                      X_e(f,:) = dPhi_dx
!!$                      Y_e(f,:) = dPhi_dy
!!$                      Z_e(f,:) = dPhi_dz
!!$                   end where
!!$                end if
!!$             end do
!!$
!!$             ! Dirichlet BCs present
!!$             if (PRESENT(Weight_Data)) then
!!$                call LSLR (X_v, Y_v, Z_v, Phi_v, X_e, Y_e, Z_e, Phi_e,                 &
!!$                           Cell%Centroid(1), Cell%Centroid(2), Cell%Centroid(3),       &
!!$                           Phi, dPhi_dx, dPhi_dy, dPhi_dz, Phi_Bit = Phi_Bit,          &
!!$                           dirichlet_flag = dirichlet_bc, velocity_flag = velocity_bc, &
!!$                           Weight_Data = Weight_Data)
!!$             else
!!$                call LSLR (X_v, Y_v, Z_v, Phi_v, X_e, Y_e, Z_e, Phi_e, &
!!$                           Cell%Centroid(1), Cell%Centroid(2), Cell%Centroid(3), &
!!$                           Phi, dPhi_dx, dPhi_dy, dPhi_dz, Phi_Bit = Phi_Bit, &
!!$                           dirichlet_flag = dirichlet_bc, velocity_flag = velocity_bc)
!!$             end if
!!$          else
             ! No Dirichlet BCs present
             if (PRESENT(Weight_Data)) then
                call LSLR (X_v, Y_v, Z_v, Phi_v, X_e, Y_e, Z_e, Phi_e,               &
                           Cell%Centroid(1), Cell%Centroid(2), Cell%Centroid(3),     &
                           Phi, dPhi_dx, dPhi_dy, dPhi_dz, Weight_Data = Weight_Data)
             else
                call LSLR (X_v, Y_v, Z_v, Phi_v, X_e, Y_e, Z_e, Phi_e,               &
                           Cell%Centroid(1), Cell%Centroid(2), Cell%Centroid(3),     &
                           Phi, dPhi_dx, dPhi_dy, dPhi_dz)
             end if
!!$          end if

          ! Deallocate arrays
          call ARRAYDESTROY (X_e, 'Array X_e(nfc,ncells)')
          call ARRAYDESTROY (Y_e, 'Array Y_e(nfc,ncells)')
          call ARRAYDESTROY (Z_e, 'Array Z_e(nfc,ncells)')
          call ARRAYDESTROY (Phi_e, 'Array Phi_e(nfc,ncells)')
    end select

    ! Deallocate arrays
    select case (gradient_method)
       case default

       case ('volume-average', 'lslr')
          call ARRAYDESTROY (X_v, 'Array X_v(nvc,ncells)')
          call ARRAYDESTROY (Y_v, 'Array Y_v(nvc,ncells)')
          call ARRAYDESTROY (Z_v, 'Array Z_v(nvc,ncells)')
    end select

    ! Eliminate noise
    dPhi_dx = MERGE(0.0_r8, dPhi_dx, ABS(dPhi_dx) <= alittle)
    dPhi_dy = MERGE(0.0_r8, dPhi_dy, ABS(dPhi_dy) <= alittle)
    dPhi_dz = MERGE(0.0_r8, dPhi_dz, ABS(dPhi_dz) <= alittle)

  END SUBROUTINE GRADIENT

  SUBROUTINE LSLR (X_v, Y_v, Z_v, Phi_v, X_e, Y_e, Z_e, Phi_e, X, Y, Z, Phi, &
                   dPhi_dx, dPhi_dy, dPhi_dz, Phi_c, face, Phi_Bit,          &
                   dirichlet_flag, velocity_flag, Weight_Data)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the least-squares linear reconstruction gradient of Phi
    !   at either cell centers or cell faces. (X_v,Y_v,Z_v) are the
    !   cell vertex coordinates, Phi_v are the vertex values of Phi,
    !   (X_e,Y_e,Z_e) are the cell neighbor centroid coordinates, Phi_e
    !   are the cell neighbor values of Phi, and (X,Y,Z) is the point where
    !   the gradient of Phi is desired.
    !
    !=======================================================================
    use ArrayAllocate_Module, only: ARRAYCREATE, ARRAYDESTROY
    use cutoffs_module,       only: alittle
    use mesh_module,          only: Cell, Face_Vrtx, Mesh, Vrtx_Face
    use parameter_module,     only: ncells, nfc, nnodes, nvc

    ! Arguments
    real(r8), dimension(nvc,ncells), intent(IN)  :: X_v, Y_v, Z_v, Phi_v
    real(r8), dimension(nfc,ncells), intent(IN)  :: X_e, Y_e, Z_e, Phi_e
    real(r8), dimension(ncells), intent(IN)  :: X, Y, Z, Phi
    real(r8), dimension(ncells), optional, intent(IN) :: Phi_c
    integer, dimension(nfc), optional, intent(IN) :: Phi_Bit
    integer, optional, intent(IN) :: face
    logical, optional, intent(IN) :: dirichlet_flag
    logical, optional, intent(IN) :: velocity_flag
    real(r8), dimension(ncells), optional, intent(IN) :: Weight_Data
    real(r8), dimension(ncells), intent(OUT) :: dPhi_dx
    real(r8), dimension(ncells), intent(OUT) :: dPhi_dy
    real(r8), dimension(ncells), intent(OUT) :: dPhi_dz

    ! Local Variables
    logical :: dirichlet_bc, gradient_face
    logical, dimension(ncells) :: Mask
    integer :: f, ff = 0, v

    real(r8), dimension(ncells) :: Axx, Axy, Axz
    real(r8), dimension(ncells) :: Ayy, Ayz, Azz
    real(r8), dimension(ncells) :: Bx, By, Bz
    real(r8), dimension(ncells) :: dX, dY, dZ
    real(r8), dimension(ncells) :: W
    real(r8), pointer, dimension(:)   :: WD_Vrtx
    real(r8), pointer, dimension(:,:) :: WD_Vrtx_Cell, Weight_Data_Nbr
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! See if this is to be a face or cell-centered gradient
    gradient_face = PRESENT(face)
    if (gradient_face) ff = face

    ! See if there are Dirichlet BCs to be applied
    dirichlet_bc = PRESENT(dirichlet_flag) .or. PRESENT(velocity_flag)

    ! Construct the normal equation matrix A and RHS vector B
    ! Add in those vertices that are not on a boundary face
    Axx = 0.0_r8; Axy = 0.0_r8; Axz = 0.0_r8
    Ayy = 0.0_r8; Ayz = 0.0_r8; Azz = 0.0_r8
    Bx  = 0.0_r8; By  = 0.0_r8; Bz  = 0.0_r8
    Mask = .false.

    if (PRESENT(Weight_Data)) then
       call ARRAYCREATE (WD_Vrtx, 1, nnodes, 'WD_Vrtx(nnodes)')
       call ARRAYCREATE (WD_Vrtx_Cell, 1, nvc, 1, ncells, 'WD_Vrtx_Cell(nvc,ncells)')
       call ARRAYCREATE (Weight_Data_Nbr, 1, nfc, 1, ncells, 'Weight_Data_Nbr(nfc,ncells)')
       ! Cell-Vertex Quantity ...................................
       ! Vertex Average
       call VERTEX_AVG (Weight_Data, WD_Vrtx)
       ! Gather Quantity
       call EN_GATHER  (WD_Vrtx_Cell, WD_Vrtx)
       call EE_GATHER  (Weight_Data_Nbr, Weight_Data)
    end if

    do v = 1,nvc
       ! Construct deltas
       dX =  X_v(v,:) - X
       dY =  Y_v(v,:) - Y
       dZ =  Z_v(v,:) - Z
       W = dX*dX + dY*dY + dZ*dZ
       if (PRESENT(Weight_Data)) W = W * WD_Vrtx_Cell(v,:)
       W = MERGE (0.0_r8, 1.0_r8/W, ABS(W) <= alittle)

       ! Mask off the boundary vertices for cell-centered gradients
       ! and those vertices not on this face for face-centered gradients
       if (gradient_face) then
          Mask = v == Face_Vrtx(face,1) .or. v == Face_Vrtx(face,2) .or. &
                 v == Face_Vrtx(face,3) .or. v == Face_Vrtx(face,4)
       else
          Mask = Mesh%Ngbr_cell(Vrtx_Face(v,1)) /= 0 .and. &
                 Mesh%Ngbr_cell(Vrtx_Face(v,2)) /= 0 .and. &
                 Mesh%Ngbr_cell(Vrtx_Face(v,3)) /= 0
       end if

       ! Do the sums
       where (Mask)
          Axx = Axx + dX*dX*W
          Axy = Axy + dX*dY*W
          Axz = Axz + dX*dZ*W
          Ayy = Ayy + dY*dY*W
          Ayz = Ayz + dY*dZ*W
          Azz = Azz + dZ*dZ*W

          Bx = Bx + (Phi_v(v,:) - Phi)*dX*W
          By = By + (Phi_v(v,:) - Phi)*dY*W
          Bz = Bz + (Phi_v(v,:) - Phi)*dZ*W
       end where
    end do

    ! Add in neighboring elements
    do f = 1,nfc
       ! If this is a face gradient and not the
       ! face of concern, go to the next face
       if (gradient_face .and. f /= ff) cycle

       ! Construct deltas
       dX  = X_e(f,:) - X
       dY  = Y_e(f,:) - Y
       dZ  = Z_e(f,:) - Z
       W = dX*dX + dY*dY + dZ*dZ
       if (PRESENT(Weight_Data)) W = W * Weight_Data_Nbr(f,:)
       W = MERGE (0.0_r8, 1.0_r8/W, ABS(W) <= alittle)

       ! Mask off boundaries
       Mask = Mesh%Ngbr_cell(f) /= 0

!!$       ! If this is a cell-centered gradient, correct the
!!$       ! Mask if we have Dirichlet BCs by resetting it to
!!$       ! true on Dirichlet boundary faces
!!$       if (.not. gradient_face) then
!!$          if (dirichlet_bc) then
!!$             if (velocity_flag) then
!!$                Mask = Mask .or. In_Flow (BC%Flag, Phi_Bit(f))
!!$             else if (dirichlet_flag) then
!!$                Mask = Mask .or. Dirichlet (BC%Flag, Phi_Bit(f))
!!$             end if
!!$          end if
!!$       end if

       ! Do the sums
       where (Mask)
          Axx = Axx + dX*dX*W
          Axy = Axy + dX*dY*W
          Axz = Axz + dX*dZ*W
          Ayy = Ayy + dY*dY*W
          Ayz = Ayz + dY*dZ*W
          Azz = Azz + dZ*dZ*W

          Bx = Bx + (Phi_e(f,:) - Phi)*dX*W
          By = By + (Phi_e(f,:) - Phi)*dY*W
          Bz = Bz + (Phi_e(f,:) - Phi)*dZ*W
       end where
    end do

    ! Add in cell itself if this is a face gradient
    if (gradient_face) then
       ! Construct deltas
       dX =  Cell%Centroid(1) - X
       dY =  Cell%Centroid(2) - Y
       dZ =  Cell%Centroid(3) - Z

       W = dX*dX + dY*dY + dZ*dZ
       if (PRESENT(Weight_Data)) W = W * Weight_Data
       W = MERGE (0.0_r8, 1.0_r8/W, ABS(W) <= alittle)

       ! Do the sums
       Axx = Axx + dX*dX*W
       Axy = Axy + dX*dY*W
       Axz = Axz + dX*dZ*W
       Ayy = Ayy + dY*dY*W
       Ayz = Ayz + dY*dZ*W
       Azz = Azz + dZ*dZ*W

       Bx = Bx + (Phi_c - Phi)*dX*W
       By = By + (Phi_c - Phi)*dY*W
       Bz = Bz + (Phi_c - Phi)*dZ*W
    end if

    ! Solve the 3X3 system with Cramer's rule
    ! Denominator (determinant of coefficient matrix)
    dX = Axx*(Ayy*Azz - Ayz*Ayz) + &
         Axy*(Axz*Ayz - Axy*Azz) + &
         Axz*(Axy*Ayz - Axz*Ayy)
    dX = MERGE (0.0_r8, 1.0_r8/dX, ABS(dX) <= alittle)

    ! Numerators (determinant of matrix with B inserted)
    dPhi_dx =  Bx*(Ayy*Azz - Ayz*Ayz) + &
              Axy*(Bz *Ayz - By *Azz) + &
              Axz*(By *Ayz - Bz *Ayy)
    dPhi_dy = Axx*(By *Azz - Bz *Ayz) + &
               Bx*(Axz*Ayz - Axy*Azz) + &
              Axz*(Bz *Axy - By *Axz)
    dPhi_dz = Axx*(Bz *Ayy - By *Ayz) + &
              Axy*(By *Axz - Bz *Axy) + &
               Bx*(Axy*Ayz - Axz*Ayy)

    ! Multiply by inverse determinant
    dPhi_dx = dX*dPhi_dx
    dPhi_dy = dX*dPhi_dy
    dPhi_dz = dX*dPhi_dz

    ! Destroy unneeded arrays.
    if (PRESENT(Weight_Data)) then
       call ARRAYDESTROY (WD_Vrtx, 'WD_Vrtx(nnodes)')
       call ARRAYDESTROY (WD_Vrtx_Cell, 'WD_Vrtx_Cell(nvc,ncells)')
       call ARRAYDESTROY (Weight_Data_Nbr, 'Weight_Data_Nbr(nfc,ncells)')
    end if

  END SUBROUTINE LSLR

  SUBROUTINE VERTEX_AVG (X_cell, X_vertex, BOUNDARY)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute a vertex-averaged quantity X_vertex from a cell-centered
    !   quantity X_cell.  X_vertex is an inverse-volume-weighted average
    !   of X_cell: X_vertex = SUM(X_cell/Vol)/SUM(1/Vol).
    !
    !=======================================================================
    use mesh_module,      only: Cell, Vertex
    use parameter_module, only: ncells, nnodes

    ! Arguments
    real(r8), dimension(ncells), intent(IN)  :: X_cell
    real(r8), dimension(nnodes), intent(OUT) :: X_vertex
    real(r8), dimension(:), pointer, optional :: BOUNDARY

    ! Local Variables
    real(r8), dimension(ncells) :: Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Scatter inverse-volume-weighted value of X_cell
    Tmp = X_cell/Cell%Volume
    call EN_SUM_SCATTER (X_vertex, Tmp, BOUNDARY=BOUNDARY)

    ! X_vertex: Vertex-averaged X_cell
    X_vertex = X_vertex*Vertex%Rsum_rvol

  END SUBROUTINE VERTEX_AVG


  SUBROUTINE GRADIENT_CELL (Gradient_Phi, Phi)
    !=======================================================================
    ! Purpose(s):
    ! 
    !   Compute the cell-centered gradient (Gradient_Phi) of a
    !   cell-centered scalar quantity Phi with one of the following methods:
    !   'green-gauss':
    !       With a discrete approximation to Gauss's theorem, whereby
    !       the integral of (Gradient_Phi) over the cell
    !       volume is converted to an integral of Phi over the cell
    !       surface area vector.  The area integral is approximated
    !       as a discrete sum over cell faces.  This control volume
    !       formulation is discretely conservative, i.e., adjacent
    !       face contributions will telescope, leaving only boundary
    !       contributions.
    !=======================================================================
    use cutoffs_module,       only: alittle
    use linear_module,        only: LINEAR_PROP
    use mesh_module,          only: Cell
    use parameter_module,     only: ncells, nfc, nnodes, nvc, ndim

    ! Arguments
    real(r8), dimension(ndim,ncells), intent(OUT) :: Gradient_Phi
    real(r8), dimension(ncells),      intent(IN)  :: Phi

    ! Local Variables
    integer :: i, f, n

    real(r8), dimension(ncells) :: Phi_f
    real(r8) :: Ptmp
    real(r8) :: Gtmp
    real(r8), dimension(nvc,ncells) :: Phi_v
    real(r8), dimension(nnodes) :: Phi_vtx

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Compute Phi_vtx, a vertex value of Phi
    call VERTEX_AVG (Phi, Phi_vtx)

    ! Gather vertex Phi_vtx into cell vector Phi_v
    call EN_GATHER (Phi_v, Phi_vtx)

    ! compute gradient according to Green-Gauss method
    ! Loop over faces, accumulating the product
    ! Phi_f*Face_Normal for each area vector component

    Gradient_Phi = 0.0_r8
    do f = 1,nfc
       ! Interpolate vertex values to this face
       call LINEAR_PROP (f, Phi_v, Phi_f)
       ! Accumulate the dot product
       do n = 1,ndim
          Gradient_Phi(n,:) = Gradient_Phi(n,:) + Cell%Face_Normal(n,f)*Cell%Face_Area(f)*Phi_f   
       end do
    end do
    ! Normalize by cell volume
    do i = 1,ncells
    Ptmp = 1.0_r8/Cell(i)%Volume
    do n = 1,ndim
       Gtmp = Gradient_Phi(n,i)*Ptmp
       ! Eliminate noise
       if(ABS(Gtmp) <= alittle)then
         Gradient_Phi(n,i) = 0.0_r8
       else
         Gradient_Phi(n,i) = Gtmp
       endif
    end do
    end do
  END SUBROUTINE GRADIENT_CELL

  SUBROUTINE DYNAMIC_PRESSURE_FACE_GRADIENT (facegrad, pc, ghc, ghn)
    !=======================================================================
    ! Purpose(s):
    ! 
    !   Compute the cell-centered gradient (Gradient_Phi) of a
    !   cell-centered scalar quantity Phi with one of the following methods:
    !   'green-gauss':
    !       With a discrete approximation to Gauss's theorem, whereby
    !       the integral of (Gradient_Phi) over the cell
    !       volume is converted to an integral of Phi over the cell
    !       surface area vector.  The area integral is approximated
    !       as a discrete sum over cell faces.  This control volume
    !       formulation is discretely conservative, i.e., adjacent
    !       face contributions will telescope, leaving only boundary
    !       contributions.
    !=======================================================================
    use mesh_module,          only: Cell, Mesh
    use parameter_module,     only: ncells, nfc, ndim
    use gs_module,            only: EE_GATHER 

    ! Arguments
    real(r8), dimension(ndim,nfc,ncells), intent(OUT) :: facegrad
    real(r8), dimension(ncells),  intent(IN) :: pc
    real(r8), dimension(nfc,ncells),  intent(IN) :: ghc
    real(r8), dimension(nfc,ncells),  intent(IN) :: ghn

    ! Local Variables
    integer :: i, f, n, d
    real(r8) :: dl2, dl2i, plr1, plr2
    real(r8), dimension(nfc,ncells) :: pn
    real(r8), dimension(ndim) :: dx, cn
    real(r8), dimension(ndim,nfc,ncells) :: Cell_Ngbr_Coord

    !***************************************************************************

    facegrad = 0.0_r8
    
    do i=1,ndim
       call EE_GATHER(Cell_Ngbr_Coord(i,:,:), Cell(:)%Centroid(i))
    end do
    
    call EE_GATHER(pn,pc)
    
    do n=1,ncells
       do f=1,nfc

          ! for this cell value...
          plr1 = pc(n)   + ghc(f,n)
      
          ! for neighbor cell value...
          plr2 = pn(f,n) + ghn(f,n)
          
          if(Mesh(n)%Ngbr_Cell(f) == 0) then
             do d=1,ndim
                cn(d) = Cell(n)%Face_Centroid(d,f)
             end do
          else
             do d=1,ndim
                cn(d) = Cell_Ngbr_Coord(d,f,n)
             end do
          end if
        
          do d=1,ndim
             dx(d) = (Cell(n)%Centroid(d) - cn(d))
          end do
          
          dl2 = 0.0_r8
          do d=1,ndim
             dl2 = dl2 + dx(d)*dx(d)
          end do 
          dl2i = 1.0/dl2
    
          facegrad(1,f,n) = (plr1 - plr2)*dx(1)*dl2i
          facegrad(2,f,n) = (plr1 - plr2)*dx(2)*dl2i
          facegrad(3,f,n) = (plr1 - plr2)*dx(3)*dl2i
          
          if (Mesh(n)%Ngbr_Cell(f) == 0) then
             facegrad(:,f,n) = 0.0_r8
          end if
          
       end do
    end do

  END SUBROUTINE DYNAMIC_PRESSURE_FACE_GRADIENT

  !-----------------------------------------------------------------------------

END MODULE DISCRETE_OP_MODULE
