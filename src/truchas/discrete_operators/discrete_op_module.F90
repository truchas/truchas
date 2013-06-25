MODULE DISCRETE_OP_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define various procedures to perform discrete operations
  !   such as gradient, divergence, etc.
  !
  !   Public Interface(s):
  !     * DIVERGENCE_OLD (D, Vc)
  !
  !         Returns the cell-centered divergence D of the
  !         cell-centered vector Vc(ndim) or Zone%Vc(ndim)
  !
  !     * FACE_CENTER (Cell, Mesh, Vertex, Q_Cell, Q_Face)
  !
  !         Evaluate a cell-face centroid quantity given a cell-center quantity.
  !
  !     * FACE_GRADIENT (f, dPhi_dx, dPhi_dy, dPhi_dz, Phi, Phi_n, &
  !                      X_n, Y_n, Z_n, Phi_e, X_e, Y_e, Z_e, &
  !                      Phi_Bit, Phi_BC, node_extrap, method)
  !
  !         Returns the face-f-centered gradient (dPhi_dx,dPhi_dy,dPhi_dz)
  !         of a cell-centered quantity Phi. Phi_e and Phi_n are vectors
  !         containing values of Phi at element and node neighbors, 
  !         respectively.
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
  !           DIVERGENCE_OLD
  !           FACE_CENTER
  !           FACE_GRADIENT
  !           FACE_GRADIENT_T
  !           FACE_GRADIENT_SAHOTA
  !           GRADIENT
  !           LSLR
  !           LSLR_T
  !           NODE_EXTRAPOLATION
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

  use gs_module,            only: EN_GATHER, EE_GATHER, EN_SUM_SCATTER

  implicit none

  ! Private Module
  private

  ! Public Subroutines
  public :: DETERMINANT_VOL_AVG, &
            DIVERGENCE_OLD,          &
            FACE_CENTER,         &
            FACE_GRADIENT,       &
            FACE_GRADIENT_T,     &
            GRADIENT,            &
            GRADIENT_CELL,       &
            VERTEX_AVG,          &
            GRADIENT_CELL_FROM_FACE_VALUES, &
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
    use constants_module, only: one_twelfth, zero
    use kind_module,      only: int_kind, real_kind
    use parameter_module, only: ncells, nfc, nvc

    implicit none

    ! Arguments
    real(KIND = real_kind), dimension(nvc,ncells), intent(IN)  :: Xv, Yv, Zv
    real(KIND = real_kind), dimension(ncells),     intent(OUT) :: Avg

    ! Local Variables
    integer(KIND = int_kind) :: f
    integer(KIND = int_kind) :: v1, v2, v3, v4, v5, v6

    real(KIND = real_kind), dimension(ncells) :: X1, Y1, Z1
    real(KIND = real_kind), dimension(ncells) :: X2, Y2, Z2
    real(KIND = real_kind), dimension(ncells) :: X3, Y3, Z3

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant quantities
    Avg = zero

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

    Avg = one_twelfth*Avg

    return

  END SUBROUTINE DETERMINANT_VOL_AVG

  SUBROUTINE DIVERGENCE_OLD (D, Vc, Weight_Data)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the cell-centered divergence (D) of a cell-centered
    !   vector Vc with a discrete approximation
    !   to Gauss's theorem, whereby the integral of Div*V over the
    !   cell volume is converted to an integral of V over the cell
    !   surface area vector.  The area integral is approximated as a
    !   discrete sum over cell faces. This control volume formulation
    !   is discretely conservative, i.e., adjacent face contributions
    !   telescope, leaving only boundary contributions.
    !
    !=======================================================================
    use bc_module,         only: BC, DIRICHLET, SET_VELOCITY_BC, Vel, Prs
    use constants_module,  only: zero
    use cutoffs_module,    only: alittle
    use fluid_data_module, only: Fluxing_Velocity
    use kind_module,       only: int_kind, log_kind, real_kind
    use linear_module,     only: LINEAR_PROP
    use mesh_module,       only: Cell, Mesh, Vrtx_Face, Vertex
    use parameter_module,  only: ncells, ndim, nfc, nfv, nnodes, nvc
    use zone_module,       only: Zone

    implicit none

    ! Arguments
    real(KIND = real_kind),   dimension(ndim,ncells), optional, intent(IN)  :: Vc
    real(KIND = real_kind),   dimension(ncells),      optional, intent(IN)  :: Weight_Data
    real(KIND = real_kind),   dimension(ncells),                intent(OUT) :: D

    ! Local Variables
    logical(KIND = log_kind), dimension(ncells)                             :: Mask
    integer(KIND = int_kind)                                                :: f, face, n, v, i
    real(KIND = real_kind),   dimension(ncells)                             :: Tmp, Velocity
    real(KIND = real_kind),   dimension(ndim,nfc,ncells)                    :: V_Face
    real(KIND = real_kind),   dimension(nvc,ncells)                         :: V_v
    real(KIND = real_kind),   dimension(nnodes)                             :: V_vtx

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Loop over each velocity component, computing a face-centered
    ! velocity component from cell-centered velocities
    do n = 1,ndim
       ! Loop over each vertex, computing the portion of the
       ! cell-centered velocity component apportioned to each vertex
       V_v = zero
       do v = 1,nvc
          ! Select the appropriate component
          if (PRESENT(Vc)) then
             Velocity = Vc(n,:)
          else
             Velocity = Zone%Vc(n)
          end if

          ! Loop over those faces associated with this vertex,
          ! and subtract out the face-normal part of the velocity
          ! component for those faces that are boundary faces. This
          ! operation mocks up the presence of a ghost cell across
          ! from each boundary face that has an equal and opposite
          ! cell-centered velocity. Note: nfv - number of faces/vertex,
          ! Vrtx_Face(v,f) = face number associated with vertex f
          do f = 1,nfv
             face = Vrtx_Face(v,f)

             ! Set the Mask to true on boundary faces, 
             ! except for those faces with a specified pressure,
             ! where flow can go in/out
             Mask = Mesh%Ngbr_cell(face) == 0 .and. &
                    .not. (DIRICHLET (BC%Flag, Prs%Face_Bit(face)))

             ! Normal part of the velocity
             Tmp = zero
             do i = 1,ndim
                if (PRESENT(Vc)) then
                   Tmp = Tmp + Vc(i,:)*Cell%Face_Normal(i,face)
                else
                   Tmp = Tmp + Zone%Vc(i)*Cell%Face_Normal(i,face)
                end if
             end do

             ! Subtract out the normal part
             where (Mask)
                Velocity = Velocity - Tmp*Cell%Face_Normal(n,face)
             end where
          end do

          ! Weight the contribution of the velocity to this
          ! vertex with the inverse cell volume
          V_v(v,:) = Velocity/Cell%Volume
       end do

       ! Sum-Scatter the velocity vector to each vertex
       call EN_SUM_SCATTER (V_vtx, V_v)

       ! Multiply by the inverse volume sum
       V_Vtx = V_Vtx * Vertex%Rsum_Rvol

       ! Gather V_Vtx into cell vector V_v
       call EN_GATHER (V_v, V_Vtx)

       ! Loop over faces, computing a face-centered
       ! vector from the vertex-averaged velocities
       do f = 1,nfc
          ! Interpolate vertex values to this face
          call LINEAR_PROP (f, V_v, V_face(n,f,:))
       end do
    end do

    ! Loop over faces and accumulate the divergence
    ! from these face-centered velocities
    D = zero
    do f = 1,nfc

       Tmp = zero

       ! Apply BC on the face velocities
       call SET_VELOCITY_BC (V_Face(:,f,:), Vel%Face_bit(f), f)

       ! Mask off dirichlet pressure faces
       Mask = DIRICHLET (BC%Flag, Prs%Face_Bit(f))

       ! Dot product of face velocities and unit normals
       do n = 1,ndim
          where (Mask) V_Face(n,f,:) = Fluxing_Velocity(f,:)*Cell%Face_Normal(n,f)
          Tmp = Tmp + V_Face(n,f,:)*Cell%Face_Normal(n,f)
       end do

       ! Multiply by area
       D = D + Tmp*Cell%Face_Area(f)
    end do

    ! Normalize by cell volume
    D = D/Cell%Volume

    ! Eliminate noise
    D = MERGE(zero, D, ABS(D) <= alittle)

    return

  END SUBROUTINE DIVERGENCE_OLD

  SUBROUTINE FACE_GRADIENT (f, dPhi_dx, dPhi_dy, dPhi_dz, Phi, &
                            Phi_n, X_n, Y_n, Z_n,              &
                            Phi_e, X_e, Y_e, Z_e,              &
                            node_extrap, Weight_Data,          &
                            method, Phi_Bit, Phi_BC)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face-centered gradient (dPhi_dx,dPhi_dy,dPhi_dz)
    !   of a cell-centered scalar quantity Phi on face f. The face
    !   gradient is computed by solving a 3x3 linear system based on
    !   cell-centered values of Phi across face f, node-averaged
    !   values of Phi (Phi_n) at each of the nvf nodes on face f,
    !   node coordinates (X_n,Y_n,Z_n), and element coordinates
    !   (X_e,Y_e,Z_e) across face f, and Phi across face f (Phi_e).
    !
    !=======================================================================
    use bc_module,        only: Vel
    use constants_module, only: zero
    use cutoffs_module,   only: alittle
    use kind_module,      only: int_kind, log_kind, real_kind
    use linear_module,    only: LINEAR_PROP
    use parameter_module, only: ncells, nfc, nvc


    implicit none

    ! Arguments
    integer(KIND = int_kind),                                  intent(IN)    :: f
    real(KIND = real_kind),   dimension(ncells),               intent(IN)    :: Phi
    real(KIND = real_kind),   dimension(nfc,ncells),           intent(IN)    :: X_e
    real(KIND = real_kind),   dimension(nfc,ncells),           intent(IN)    :: Y_e
    real(KIND = real_kind),   dimension(nfc,ncells),           intent(IN)    :: Z_e
    real(KIND = real_kind),   dimension(nfc,ncells),           intent(IN)    :: Phi_e
    real(KIND = real_kind),   dimension(nvc,ncells),           intent(IN)    :: X_n
    real(KIND = real_kind),   dimension(nvc,ncells),           intent(IN)    :: Y_n
    real(KIND = real_kind),   dimension(nvc,ncells),           intent(IN)    :: Z_n
    logical(KIND = log_kind),                                  intent(IN)    :: node_extrap
    character(LEN = *),                              optional, intent(IN)    :: method
    real(KIND = real_kind),   dimension(ncells),     optional, intent(IN)    :: Weight_Data
    integer(KIND = int_kind), dimension(nfc),        optional, intent(IN)    :: Phi_Bit
    real(KIND = real_kind),   dimension(nfc,ncells), optional, intent(IN)    :: Phi_BC
    real(KIND = real_kind),   dimension(nvc,ncells),           intent(INOUT) :: Phi_n
    real(KIND = real_kind),   dimension(ncells),               intent(OUT)   :: dPhi_dx
    real(KIND = real_kind),   dimension(ncells),               intent(OUT)   :: dPhi_dy
    real(KIND = real_kind),   dimension(ncells),               intent(OUT)   :: dPhi_dz

    ! Local Variables
    character(LEN = 80)      :: gradient_method
    logical(KIND = log_kind) :: specified_method

    logical(KIND = log_kind) :: applied_bc
    logical(KIND = log_kind) :: fixed_bc
    logical(KIND = log_kind) :: velocity_bc
    logical(KIND = log_kind) :: zero_bc
    real(KIND = real_kind), dimension(ncells) :: Phi_f
    real(KIND = real_kind), dimension(ncells) :: X_f, Y_f, Z_f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Specify Face Gradient Method ...........................
    ! Default Gradient Method
    gradient_method = 'lslr'

    ! Other Gradient Method
    specified_method = PRESENT(method)
    if (specified_method) gradient_method = method

    ! Boundary Corrections ...................................
    applied_bc = PRESENT(Phi_Bit) ! Any Applied BCs

    ! Type of Applied BCs
    if (applied_bc) then
       ! Dirichlet with Zero Constant
       ! (This is No-Slip for Velocity BCs)
       zero_bc = .not. PRESENT(Phi_BC)

       ! Dirichlet with Non-Zero Constant
       ! (Tradtional Dirichlet BC. Dirichlet and InFlow for Velocity BCs)
       fixed_bc = PRESENT(Phi_BC)

       ! Velocity Boundary Conditions
       ! (Equivalency in First Face Bit Location)
       velocity_bc = Vel%Face_Bit(1) - Phi_Bit(1) == 0
    end if

    ! Extrapolate Boundary Node Values
    if (node_extrap) then
       if (.not. applied_bc) then
          call NODE_EXTRAPOLATION (Phi, Phi_n, X_n, Y_n, Z_n, dirichlet_bc = .false.)
       else if (zero_bc) then
          call NODE_EXTRAPOLATION (Phi, Phi_n, X_n, Y_n, Z_n, dirichlet_bc = zero_bc, &
                                   Phi_Bit = Phi_Bit) 
       else if (fixed_bc) then
          call NODE_EXTRAPOLATION (Phi, Phi_n, X_n, Y_n, Z_n, dirichlet_bc = fixed_bc, &
                                   Phi_Bit = Phi_Bit, Phi_BC  = Phi_BC)
       end if
    end if

    ! Select Face Gradient Method ............................
    select case (gradient_method)

       ! Sahota Method ..................................
       case ('Sahota')
          call FACE_GRADIENT_SAHOTA (f, Phi, Phi_e, X_e, Y_e, Z_e, &
                                             Phi_n, X_n, Y_n, Z_n, &
                                             dPhi_dx, dPhi_dy, dPhi_dz)

       ! LSLR Method
       case ('lslr')

          ! Cell-Face Coordinates and Phi ..........
          call LINEAR_PROP (f, X_n, X_f)
          call LINEAR_PROP (f, Y_n, Y_f)
          call LINEAR_PROP (f, Z_n, Z_f)
          call LINEAR_PROP (f, Phi_n, Phi_f)

!!$          ! Boundary Face Correction ...............
!!$          if (applied_bc) then
!!$             ! Velocity BCs
!!$             if (velocity_bc) then
!!$                if (fixed_bc) then ! Dirichlet and InFlow
!!$!                  Mask = DIRICHLET (BC%Flag, Phi_Bit(f)) ! Not enough bit locations
!!$!                  where (Mask) Phi_f(:) = Phi_BC(f,:)    ! Considered as InFlow BC
!!$                   Mask = IN_FLOW (BC%Flag, Phi_Bit(f))
!!$                   where (Mask) Phi_f(:) = Phi_BC(f,:)
!!$                   Mask = NO_SLIP (BC%Flag, Phi_Bit(f)) ! Consider No-Slip BC
!!$                   where (Mask) Phi_f(:) = zero         ! to be 'Fixed Velocity BC'
!!$                end if
!!$                if (zero_bc) then  ! No-Slip
!!$                   Mask = NO_SLIP (BC%Flag, Phi_Bit(f))
!!$                   where (Mask) Phi_f(:) = zero
!!$                end if
!!$             end if
!!$
!!$             ! Pressure BCs
!!$             if (.not. velocity_bc) then
!!$                if (fixed_bc) then ! Dirichlet
!!$                   Mask = DIRICHLET (BC%Flag, Phi_Bit(f))
!!$                   where (Mask) Phi_f(:) = Phi_BC(f,:)
!!$                end if
!!$             end if
!!$          end if

          ! LSLR Phi Gradient ......................
          if (applied_bc) then
             ! Applied BCs Present
             if (PRESENT(Weight_Data)) then
                call LSLR (X_n, Y_n, Z_n, Phi_n, X_e, Y_e, Z_e, Phi_e,             &
                           X_f, Y_f, Z_f, Phi_f, dPhi_dx, dPhi_dy, dPhi_dz,        &
                           Phi_c = Phi, face = f, Phi_Bit = Phi_Bit,               &
                           dirichlet_flag = fixed_bc, velocity_flag = velocity_bc, &
                           weight_data = Weight_Data)
             else
                call LSLR (X_n, Y_n, Z_n, Phi_n, X_e, Y_e, Z_e, Phi_e,    &
                        X_f, Y_f, Z_f, Phi_f, dPhi_dx, dPhi_dy, dPhi_dz,  &
                        Phi_c = Phi, face = f, Phi_Bit = Phi_Bit,         &
                        dirichlet_flag = fixed_bc, velocity_flag = velocity_bc)
             end if
          else
             ! No Applied BCs Present
             if (PRESENT(Weight_Data)) then
                call LSLR (X_n, Y_n, Z_n, Phi_n, X_e, Y_e, Z_e, Phi_e,      &
                           X_f, Y_f, Z_f, Phi_f, dPhi_dx, dPhi_dy, dPhi_dz, &
                           Phi_c = Phi, face = f, weight_data = Weight_Data)
             else
                call LSLR (X_n, Y_n, Z_n, Phi_n, X_e, Y_e, Z_e, Phi_e,   &
                           X_f, Y_f, Z_f, Phi_f, dPhi_dx, dPhi_dy, dPhi_dz, &
                           Phi_c = Phi, face = f)
             end if
          end if
    end select

    ! Eliminate Face Gradient Noise ..........................
    dPhi_dx = MERGE(zero, dPhi_dx, ABS(dPhi_dx) <= alittle)
    dPhi_dy = MERGE(zero, dPhi_dy, ABS(dPhi_dy) <= alittle)
    dPhi_dz = MERGE(zero, dPhi_dz, ABS(dPhi_dz) <= alittle)

    return

  END SUBROUTINE FACE_GRADIENT

  SUBROUTINE FACE_GRADIENT_SAHOTA (f, Phi, Phi_e, X_e, Y_e, Z_e, &
                                           Phi_n, X_n, Y_n, Z_n, &
                                           dPhi_dx, dPhi_dy, dPhi_dz)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face-centered gradient (dPhi_dx,dPhi_dy,dPhi_dz) of
    !   a cell-centered scalar quantity Phi on cell-face f using the
    !   'Sahota' method. The face gradient is computed by solving a 3x3
    !   linear system based on cell-centered values of Phi across face f
    !   (Phi_e) with element coordinates (X_e, Y_e, Z_e), node-averaged
    !   values of Phi at each of the nvf nodes on face f (Phi_n), and
    !   node coordinates (X_n,Y_n,Z_n).
    !
    !=======================================================================
    use ArrayAllocate_Module, only: ARRAYCREATE, ARRAYDESTROY
    use constants_module,     only: one, zero
    use cutoffs_module,       only: alittle
    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Cell, Mesh
    use parameter_module,     only: ncells, nfc, nvc

    implicit none

    ! Argument List

    integer(KIND = int_kind), intent(IN) :: f

    real(KIND = real_kind), dimension(ncells),     intent(IN)    :: Phi
    real(KIND = real_kind), dimension(nfc,ncells), intent(IN)    :: Phi_e
    real(KIND = real_kind), dimension(nfc,ncells), intent(IN)    :: X_e, Y_e, Z_e
    real(KIND = real_kind), dimension(nvc,ncells), intent(IN)    :: X_n, Y_n, Z_n
    real(KIND = real_kind), dimension(nvc,ncells), intent(INOUT) :: Phi_n
    real(KIND = real_kind), dimension(ncells),     intent(INOUT)   :: dPhi_dx, dPhi_dy, dPhi_dz

    ! Local Variables
    integer(KIND = int_kind)                      :: v1, v2, v3, v4,i
    real (KIND = real_kind)                       :: dist, dot
    real(KIND = real_kind), dimension(:), pointer :: D, dPhi
    real(KIND = real_kind), dimension(:), pointer :: dX, dY, dZ

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Allocate Arrays
    call ARRAYCREATE (D,    1, ncells, 'Array D(ncells)')
    call ARRAYCREATE (dPhi, 1, ncells, 'Array dPhi(ncells)')

    call ARRAYCREATE (dX, 1, ncells, 'Array dX(ncells)')
    call ARRAYCREATE (dY, 1, ncells, 'Array dY(ncells)')
    call ARRAYCREATE (dZ, 1, ncells, 'Array dZ(ncells)')

    ! Compute Deltas (Cell to Face Neighbors)
    ! Physical Coordinate Deltas
    dX(:) = Cell(:)%Centroid(1) - X_e(f,:)
    dY(:) = Cell(:)%Centroid(2) - Y_e(f,:)
    dZ(:) = Cell(:)%Centroid(3) - Z_e(f,:)

    ! Cell-Center Phi Deltas
    where (Mesh%Ngbr_Cell(f) == 0)
       dPhi = 0
    elsewhere
       dPhi = Phi - Phi_e(f,:)
    end where

    ! Cell-Center Phi Deltas
    !dPhi(:) = Phi(:) - Phi_e(f,:)

    if(.true.)then
      !EDM  Compute a simple cell face gradient a la CFDLIB
      do i = 1,ncells
       dist = dX(i)**2 + dY(i)**2 + dZ(i)**2 
       dist = sqrt(dist) 
       ! dPhi_dx(i) = -dPhi(i)/(dist + 1.d-30)*Cell(i)%Face_Normal(1,f)
       ! dPhi_dy(i) = -dPhi(i)/(dist + 1.d-30)*Cell(i)%Face_Normal(2,f)
       ! dPhi_dz(i) = -dPhi(i)/(dist + 1.d-30)*Cell(i)%Face_Normal(3,f)
        dPhi_dx(i) = dPhi(i)/(dist + 1.d-30)**2*dX(i)/2.d0
        dPhi_dy(i) = dPhi(i)/(dist + 1.d-30)**2*dY(i)/2.d0
        dPhi_dz(i) = dPhi(i)/(dist + 1.d-30)**2*dZ(i)/2.d0
      end do

      call ARRAYDESTROY (D,    'Array D(ncells)')
      call ARRAYDESTROY (dPhi, 'Array dPhi(ncells)')
  
      call ARRAYDESTROY (dX, 'Array dX(ncells)')
      call ARRAYDESTROY (dY, 'Array dX(ncells)')
      call ARRAYDESTROY (dZ, 'Array dX(ncells)')
   
      return
    end if

    if(.false.)then
    ! Cell-Face Vertex Numbers
    select case(f)
       case default
       case (1) ! Face One
          v1 = 4; v2 = 8; v3 = 7; v4 = 3
       case (2) ! Face Two
          v1 = 1; v2 = 2; v3 = 6; v4 = 5
       case (3) ! Face Three
          v1 = 4; v2 = 1; v3 = 5; v4 = 8
       case (4) ! Face Four
          v1 = 3; v2 = 7; v3 = 6; v4 = 2
       case (5) ! Face Five
          v1 = 4; v2 = 3; v3 = 2; v4 = 1
       case (6) ! Face Six
          v1 = 8; v2 = 5; v3 = 6; v4 = 7
    end select

    ! Denominator
    D(:) = dX(:)*((Y_n(v1,:) - Y_n(v3,:)) * (Z_n(v2,:) - Z_n(v4,:))  - &
                  (Y_n(v2,:) - Y_n(v4,:)) * (Z_n(v1,:) - Z_n(v3,:))) + &
           dY(:)*((X_n(v2,:) - X_n(v4,:)) * (Z_n(v1,:) - Z_n(v3,:))  - &
                  (X_n(v1,:) - X_n(v3,:)) * (Z_n(v2,:) - Z_n(v4,:))) + &
           dZ(:)*((X_n(v1,:) - X_n(v3,:)) * (Y_n(v2,:) - Y_n(v4,:))  - &
                  (X_n(v2,:) - X_n(v4,:)) * (Y_n(v1,:) - Y_n(v3,:)))
    D = MERGE(zero, one/D, ABS(D) <= alittle)

    ! X-Derivative: dPhi_dx
    dPhi_dx(:) = dPhi(:) * ((Y_n(v1,:) - Y_n(v3,:)) * (Z_n(v2,:) - Z_n(v4,:))      - &
                            (Y_n(v2,:) - Y_n(v4,:)) * (Z_n(v1,:) - Z_n(v3,:)))     + &
                   dY(:) * ((Phi_n(v2,:) - Phi_n(v4,:)) * (Z_n(v1,:) - Z_n(v3,:))  - &
                            (Phi_n(v1,:) - Phi_n(v3,:)) * (Z_n(v2,:) - Z_n(v4,:))) + &
                   dZ(:) * ((Phi_n(v1,:) - Phi_n(v3,:)) * (Y_n(v2,:) - Y_n(v4,:))  - &
                            (Phi_n(v2,:) - Phi_n(v4,:)) * (Y_n(v1,:) - Y_n(v3,:)))

    ! Y-Derivative: dPhi_dy
    dPhi_dy(:) =   dX(:) * ((Phi_n(v1,:) - Phi_n(v3,:)) * (Z_n(v2,:) - Z_n(v4,:))  - &
                            (Phi_n(v2,:) - Phi_n(v4,:)) * (Z_n(v1,:) - Z_n(v3,:))) + &
                 dPhi(:) * ((X_n(v2,:) - X_n(v4,:)) * (Z_n(v1,:) - Z_n(v3,:))      - &
                            (X_n(v1,:) - X_n(v3,:)) * (Z_n(v2,:) - Z_n(v4,:)))     + &
                   dZ(:) * ((X_n(v1,:) - X_n(v3,:)) * (Phi_n(v2,:) - Phi_n(v4,:))  - &
                            (X_n(v2,:) - X_n(v4,:)) * (Phi_n(v1,:) - Phi_n(v3,:)))

    ! Z-Derivative: dPhi_dz
    dPhi_dz(:) =   dX(:) * ((Y_n(v1,:) - Y_n(v3,:)) * (Phi_n(v2,:) - Phi_n(v4,:))  - &
                            (Y_n(v2,:) - Y_n(v4,:)) * (Phi_n(v1,:) - Phi_n(v3,:))) + &
                   dY(:) * ((X_n(v2,:) - X_n(v4,:)) * (Phi_n(v1,:) - Phi_n(v3,:))  - &
                            (X_n(v1,:) - X_n(v3,:)) * (Phi_n(v2,:) - Phi_n(v4,:))) + &
                 dPhi(:) * ((X_n(v1,:) - X_n(v3,:)) * (Y_n(v2,:) - Y_n(v4,:))      - &
                            (X_n(v2,:) - X_n(v4,:)) * (Y_n(v1,:) - Y_n(v3,:)))

    ! Denominator Correction
    dPhi_dx(:) = D(:) * dPhi_dx(:)
    dPhi_dy(:) = D(:) * dPhi_dy(:)
    dPhi_dz(:) = D(:) * dPhi_dz(:)
    endif


    !subtract out component along cell centers and replace with a compact treatment
    do i = 1,ncells
     dist = dX(i)**2 + dY(i)**2 + dZ(i)**2 
     dot =  dPhi_dx(i)*dX(i) + dPhi_dy(i)*dY(i) + dPhi_dZ(i)*dZ(i)
     dPhi_dX(i) = dPhi_dX(i) +  (dPhi(i) + dot)/(dist + 1.d-30)*dX(i)  
     dPhi_dY(i) = dPhi_dY(i) +  (dPhi(i) + dot)/(dist + 1.d-30)*dY(i)  
     dPhi_dZ(i) = dPhi_dZ(i) +  (dPhi(i) + dot)/(dist + 1.d-30)*dZ(i)  

     dPhi_dX(i) = dPhi(i)/(dist + 1.d-30)*dX(i)  
     dPhi_dY(i) = dPhi(i)/(dist + 1.d-30)*dY(i)  
     dPhi_dZ(i) = dPhi(i)/(dist + 1.d-30)*dZ(i)  
    end do
    
    ! Deallocate Arrays
    call ARRAYDESTROY (D,    'Array D(ncells)')
    call ARRAYDESTROY (dPhi, 'Array dPhi(ncells)')

    call ARRAYDESTROY (dX, 'Array dX(ncells)')
    call ARRAYDESTROY (dY, 'Array dX(ncells)')
    call ARRAYDESTROY (dZ, 'Array dX(ncells)')

    return

  END SUBROUTINE FACE_GRADIENT_SAHOTA

  SUBROUTINE FACE_CENTER (Q_Cell, Q_Face)
    !=======================================================================
    ! Purpose(s):
    !   Evaluate a cell-face centroid quantity given a cell-center quantity.
    !=======================================================================
    use constants_module,   only: zero
    use kind_module,        only: int_kind, real_kind
    use linear_module,      only: LINEAR_PROP
    use parameter_module,   only: ncells, nfc, nnodes, nvc

    implicit none

    ! Argument List

    real(KIND = real_kind), dimension(ncells),     intent(IN)    :: Q_Cell
    real(KIND = real_kind), dimension(nfc,ncells), intent(OUT)   :: Q_Face

    ! Local Variables
    integer(KIND = int_kind) :: f
    real(KIND = real_kind), dimension(nnodes)     :: Q_Vrtx
    real(KIND = real_kind), dimension(nvc,ncells) :: Q_Vrtx_Cell

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Cell-Vertex Quantity ...................................
    ! Vertex Average
    call VERTEX_AVG (Q_Cell, Q_Vrtx)

    ! Gather Quantity
    call EN_GATHER  (Q_Vrtx_Cell, Q_Vrtx)

    ! Cell-Face Quantity .....................................
    ! Initialize Quantity
    Q_Face = zero

    ! Evaluate Quantity
    do f = 1, nfc
       call LINEAR_PROP (f, Q_Vrtx_Cell, Q_Face(f,:))
    end do

    return

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
    use constants_module,     only: one, zero
    use cutoffs_module,       only: alittle
    use kind_module,          only: int_kind, log_kind, real_kind
    use linear_module,        only: LINEAR_PROP
    use mesh_module,          only: Cell, Vertex, Vrtx_Bdy
    use parameter_module,     only: ncells, nfc, nnodes, nvc
    implicit none

    ! Arguments
    character(LEN = *), intent(IN), optional :: method
    real(KIND = real_kind), dimension(ncells),       intent(OUT)          :: dPhi_dx
    real(KIND = real_kind), dimension(ncells),       intent(OUT)          :: dPhi_dy
    real(KIND = real_kind), dimension(ncells),       intent(OUT)          :: dPhi_dz
    real(KIND = real_kind),   dimension(ncells),     intent(IN)           :: Phi
    real(KIND = real_kind),   dimension(ncells),     intent(IN), optional :: Weight_Data
    real(KIND = real_kind),   dimension(nfc,ncells), intent(IN), optional :: Phi_BC
    integer(KIND = int_kind), dimension(nfc),        intent(IN), optional :: Phi_Bit

    ! Local Variables
    character(LEN = 80) :: gradient_method = 'volume-average'
    logical(KIND = log_kind) :: specified_method
    logical(KIND = log_kind) :: dirichlet_bc, velocity_bc
    integer(KIND = int_kind)                    :: f

    real(KIND = real_kind), dimension(ncells)              :: Phi_f
    real(KIND = real_kind), dimension(:,:),        pointer :: Phi_e
    real(KIND = real_kind), dimension(nvc,ncells)          :: Phi_v
    real(KIND = real_kind), dimension(nnodes)              :: Phi_vtx

    real(KIND = real_kind), dimension(:,:),        pointer :: X_v, Y_v, Z_v
    real(KIND = real_kind), dimension(:,:),        pointer :: X_e, Y_e, Z_e

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
          Phi_f = one/Cell%Volume
          dPhi_dx = dPhi_dx*Phi_f
          dPhi_dy = dPhi_dy*Phi_f
          dPhi_dz = dPhi_dz*Phi_f

       case ('Green-Gauss')
          ! Loop over faces, accumulating the product
          ! Phi_f*Face_Normal for each area vector component
          dPhi_dx = zero; dPhi_dy = zero; dPhi_dz = zero
          do f = 1,nfc
             ! Interpolate vertex values to this face
             call LINEAR_PROP (f, Phi_v, Phi_f)

             ! Accumulate the dot product
             dPhi_dx = dPhi_dx + Cell%Face_Normal(1,f)*Cell%Face_Area(f)*Phi_f
             dPhi_dy = dPhi_dy + Cell%Face_Normal(2,f)*Cell%Face_Area(f)*Phi_f
             dPhi_dz = dPhi_dz + Cell%Face_Normal(3,f)*Cell%Face_Area(f)*Phi_f
          end do

          ! Normalize by cell volume
          Phi_f = one/Cell%Volume
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
    dPhi_dx = MERGE(zero, dPhi_dx, ABS(dPhi_dx) <= alittle)
    dPhi_dy = MERGE(zero, dPhi_dy, ABS(dPhi_dy) <= alittle)
    dPhi_dz = MERGE(zero, dPhi_dz, ABS(dPhi_dz) <= alittle)

    return

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
    use constants_module,     only: one, zero
    use cutoffs_module,       only: alittle
    use kind_module,          only: int_kind, log_kind, real_kind
    use mesh_module,          only: Cell, Face_Vrtx, Mesh, Vrtx_Face
    use parameter_module,     only: ncells, nfc, nnodes, nvc

    implicit none

    ! Arguments
    real(KIND = real_kind),   dimension(nvc,ncells),         intent(IN)  :: X_v, Y_v, Z_v, Phi_v
    real(KIND = real_kind),   dimension(nfc,ncells),         intent(IN)  :: X_e, Y_e, Z_e, Phi_e
    real(KIND = real_kind),   dimension(ncells),             intent(IN)  :: X, Y, Z, Phi
    real(KIND = real_kind),   dimension(ncells),   optional, intent(IN)  :: Phi_c
    integer(KIND = int_kind), dimension(nfc),      optional, intent(IN)  :: Phi_Bit
    integer(KIND = int_kind),                      optional, intent(IN)  :: face
    logical(KIND = log_kind),                      optional, intent(IN)  :: dirichlet_flag
    logical(KIND = log_kind),                      optional, intent(IN)  :: velocity_flag
    real(KIND = real_kind),   dimension(ncells),   optional, intent(IN)  :: Weight_Data
    real(KIND = real_kind),   dimension(ncells),             intent(OUT) :: dPhi_dx
    real(KIND = real_kind),   dimension(ncells),             intent(OUT) :: dPhi_dy
    real(KIND = real_kind),   dimension(ncells),             intent(OUT) :: dPhi_dz

    ! Local Variables
    logical(KIND = log_kind)                    :: dirichlet_bc
    logical(KIND = log_kind)                    :: gradient_face
    logical(KIND = log_kind), dimension(ncells) :: Mask
    integer(KIND = int_kind)                    :: f, ff = 0, v

    real(KIND = real_kind), dimension(ncells)       :: Axx, Axy, Axz
    real(KIND = real_kind), dimension(ncells)       :: Ayy, Ayz, Azz
    real(KIND = real_kind), dimension(ncells)       :: Bx, By, Bz
    real(KIND = real_kind), dimension(ncells)       :: dX, dY, dZ
    real(KIND = real_kind), dimension(ncells)       :: W
    real(KIND = real_kind), pointer, dimension(:)   :: WD_Vrtx
    real(KIND = real_kind), pointer, dimension(:,:) :: WD_Vrtx_Cell, Weight_Data_Nbr
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! See if this is to be a face or cell-centered gradient
    gradient_face = PRESENT(face)
    if (gradient_face) ff = face

    ! See if there are Dirichlet BCs to be applied
    dirichlet_bc = PRESENT(dirichlet_flag) .or. PRESENT(velocity_flag)

    ! Construct the normal equation matrix A and RHS vector B
    ! Add in those vertices that are not on a boundary face
    Axx = zero; Axy = zero; Axz = zero
    Ayy = zero; Ayz = zero; Azz = zero
    Bx  = zero; By  = zero; Bz  = zero
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
       W = MERGE (zero, one/W, ABS(W) <= alittle)

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
       W = MERGE (zero, one/W, ABS(W) <= alittle)

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
       W = MERGE (zero, one/W, ABS(W) <= alittle)

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
    dX = MERGE (zero, one/dX, ABS(dX) <= alittle)

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

    return

  END SUBROUTINE LSLR

  SUBROUTINE NODE_EXTRAPOLATION (Phi, Phi_n, X_n, Y_n, Z_n, dirichlet_bc, &
                                 Phi_Bit, Phi_BC)
    !=======================================================================
    ! Purpose(s):
    !
    !   Extrapolate a cell-centered quantity Phi to boundary vertices
    !   by performing a distance-weighted Taylor series expansion from
    !   cell centers.
    !
    !=======================================================================
    use constants_module, only: one, zero
    use cutoffs_module,   only: alittle
    use kind_module,      only: int_kind, log_kind, real_kind
    use mesh_module,      only: Cell, Mesh, Vrtx_Face
    use parameter_module, only: ncells, ndim, nfc, nnodes, nvc

    implicit none

    ! Argument List

    real(KIND = real_kind), dimension(ncells),                 intent(IN)    :: Phi
    real(KIND = real_kind), dimension(nvc,ncells),             intent(IN)    :: X_n, Y_n, Z_n
    logical(KIND = log_kind),                                  intent(IN)    :: dirichlet_bc
    integer(KIND = int_kind), dimension(nfc),        optional, intent(IN)    :: Phi_Bit
    real(KIND = real_kind),   dimension(nfc,ncells), optional, intent(IN)    :: Phi_BC
    real(KIND = real_kind), dimension(nvc,ncells),             intent(OUT)   :: Phi_n

    ! Local Variables
    integer(KIND = int_kind)                         :: v, n
    logical(KIND = log_kind), dimension(ncells)      :: Mask
    real(KIND = real_kind),   dimension(ndim,ncells) :: Coord, Grad
    real(KIND = real_kind),   dimension(nvc,ncells)  :: Tmp1_n, Tmp2_n
    real(KIND = real_kind),   dimension(nnodes)      :: Vtx1, Vtx2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get the cell-centered gradient.
    if (dirichlet_bc) then
       call GRADIENT (Grad(1,:), Grad(2,:), Grad(3,:), Phi, &
                      Phi_Bit = Phi_Bit, Phi_BC = Phi_BC, &
                      method = 'lslr')
    else
       call GRADIENT (Grad(1,:), Grad(2,:), Grad(3,:), Phi, &
                      method = 'lslr')
    end if

    ! Loop over vertices, and use the cell-centered gradient
    ! to extrapolate a value of Phi to boundary vertices
    Vtx1 = zero; Vtx2 = zero
    do v = 1,nvc

       ! Distance from cell centroid to vertex
       Coord(1,:) = X_n(v,:) - Cell%Centroid(1)
       Coord(2,:) = Y_n(v,:) - Cell%Centroid(2)
       Coord(3,:) = Z_n(v,:) - Cell%Centroid(3)

       Tmp1_n(v,:) = zero
       do n = 1,ndim
          Tmp1_n(v,:) = Tmp1_n(v,:) + Coord(n,:)**2
       end do
       Tmp1_n(v,:) = MERGE (zero, one/SQRT(Tmp1_n(v,:)), Tmp1_n(v,:) <= alittle)

       ! Taylor series expansion to this vertex.
       Tmp2_n(v,:) = Phi
       do n = 1,ndim
          Tmp2_n(v,:) = Tmp2_n(v,:) + Coord(n,:)*Grad(n,:)
       end do

       ! Distance-weight the expansion
       Tmp2_n(v,:) = Tmp1_n(v,:)*Tmp2_n(v,:)

    end do

    ! Scatter the weighted expansions and accumulate them.
    call EN_SUM_SCATTER (Vtx2, Tmp2_n)

    ! Scatter the inverse distances and accumulate them.
    call EN_SUM_SCATTER (Vtx1, Tmp1_n)

    ! Compute the average expansion at vertices; gather it.
    Vtx2 = MERGE (zero, Vtx2/Vtx1, ABS(Vtx1) <= alittle)
    call EN_GATHER (Tmp1_n, Vtx2)

    ! Only take the expansion at boundary vertices
    do v = 1,nvc
       Mask = .false.
       where (Mesh%Ngbr_cell(Vrtx_Face(v,1)) == 0 .or. &
              Mesh%Ngbr_cell(Vrtx_Face(v,2)) == 0 .or. &
              Mesh%Ngbr_cell(Vrtx_Face(v,3)) == 0)
              Phi_n(v,:) = Tmp1_n(v,:)
              Mask = .true.
       end where
    end do

    return

  END SUBROUTINE NODE_EXTRAPOLATION

  SUBROUTINE VERTEX_AVG (X_cell, X_vertex, BOUNDARY)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute a vertex-averaged quantity X_vertex from a cell-centered
    !   quantity X_cell.  X_vertex is an inverse-volume-weighted average
    !   of X_cell: X_vertex = SUM(X_cell/Vol)/SUM(1/Vol).
    !
    !=======================================================================
    use kind_module,      only: real_kind
    use mesh_module,      only: Cell, Vertex
    use parameter_module, only: ncells, nnodes

    implicit none

    ! Arguments
    real(KIND = real_kind),  dimension(ncells),  intent(IN)  :: X_cell
    real(KIND = real_kind),  dimension(nnodes),  intent(OUT) :: X_vertex
    real(KIND = real_kind),  dimension(:), pointer, optional :: BOUNDARY

    ! Local Variables
    real(KIND = real_kind), dimension(ncells) :: Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Scatter inverse-volume-weighted value of X_cell
    Tmp = X_cell/Cell%Volume
    call EN_SUM_SCATTER (X_vertex, Tmp, BOUNDARY=BOUNDARY)

    ! X_vertex: Vertex-averaged X_cell
    X_vertex = X_vertex*Vertex%Rsum_rvol

    return

  END SUBROUTINE VERTEX_AVG

  !-----------------------------------------------------------------------------

  Subroutine FACE_GRADIENT_T (face, dPhi_dx, dPhi_dy, dPhi_dz, Phi, &
                              Phi_n, X_n, Y_n, Z_n,                 &
                              Phi_e, X_e, Y_e, Z_e,                 &
                              node_extrap, WD,                      &
                              method_arg)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   *** special case - temperature only ***
    !
    !   Compute the face-centered gradient (dPhi_dx,dPhi_dy,dPhi_dz)
    !   of a cell-centered scalar quantity Phi on face f. The face
    !   gradient is computed by solving a 3x3 linear system based on
    !   cell-centered values of Phi across face f, node-averaged
    !   values of Phi (Phi_n) at each of the nvf nodes on face f,
    !   node coordinates (X_n,Y_n,Z_n), and element coordinates
    !   (X_e,Y_e,Z_e) across face f, and Phi across face f (Phi_e).
    !---------------------------------------------------------------------------

    use bc_module,        only: BCMATCH, BC_T_DIRICHLET, BC_T_HNEUMANN, BC_T_NEUMANN,&
                                BC_T_NO_BC, BC_T, BC_T_RADIATION, BC_T_VFRADIATION, BC_T_HTC
    use kind_module,      only: log_kind, int_kind, real_kind
    use linear_module,    only: LINEAR_PROP
    use parameter_module, only: ncells, nfc, nvc
    use truchas_logging_services

    implicit none

    ! arguments
    integer(int_kind),                                  intent(IN)    :: face
    real(real_kind),   dimension(ncells),               intent(IN)    :: Phi
    real(real_kind),   dimension(nfc,ncells),           intent(IN)    :: X_e
    real(real_kind),   dimension(nfc,ncells),           intent(IN)    :: Y_e
    real(real_kind),   dimension(nfc,ncells),           intent(IN)    :: Z_e
    real(real_kind),   dimension(nfc,ncells),           intent(IN)    :: Phi_e
    real(real_kind),   dimension(nvc,ncells),           intent(IN)    :: X_n
    real(real_kind),   dimension(nvc,ncells),           intent(IN)    :: Y_n
    real(real_kind),   dimension(nvc,ncells),           intent(IN)    :: Z_n
    logical(log_kind),                                  intent(IN)    :: node_extrap
    character(LEN=*),                         optional, intent(IN)    :: method_arg
    real(real_kind),   dimension(ncells),     optional, intent(IN)    :: WD
    real(real_kind),   dimension(nvc,ncells),           intent(IN)    :: Phi_n
    real(real_kind),   dimension(ncells),               intent(OUT)   :: dPhi_dx
    real(real_kind),   dimension(ncells),               intent(OUT)   :: dPhi_dy
    real(real_kind),   dimension(ncells),               intent(OUT)   :: dPhi_dz

    ! local variables
    character(LEN=80) :: method
    integer(int_kind) :: cell
    real(real_kind),   dimension(ncells) :: Phi_f
    real(real_kind),   dimension(ncells) :: X_f
    real(real_kind),   dimension(ncells) :: Y_f
    real(real_kind),   dimension(ncells) :: Z_f

    !---------------------------------------------------------------------------

    ! extrapolate boundary node values
    If (node_extrap) Then
       Call TLS_fatal ('FACE_GRADIENT_T: node extrapolation not implemented')
    End If

    ! select face gradient method - if not specified, use LSLR
    If (Present(method_arg)) Then
       method = method_arg
    Else
       method = 'lslr'
    End If

    ! calculate face gradients using the appropriate method
    Select Case (method)
    Case ('Sahota')
       Call TLS_fatal ('FACE_GRADIENT_T: sahota method not implemented')

    Case ('lslr')
       ! calculate coordinates of cell face center from surrounding nodal values
       Call LINEAR_PROP (face, X_n, X_f)
       Call LINEAR_PROP (face, Y_n, Y_f)
       Call LINEAR_PROP (face, Z_n, Z_f)

       ! calculate phi at cell face center from surrounding nodal values
       Call LINEAR_PROP (face, Phi_n, Phi_f)

       ! correct for Dirichlet boundary conditions (only case considered)
       Do cell = 1, ncells
          Select Case (BCMATCH(BC_T, face, cell))
          Case (BC_T_NO_BC)
          Case (BC_T_DIRICHLET)
             phi_f(cell) = BC_T%Value1(face, cell)
          Case (BC_T_HNEUMANN)
             continue                   ! do _NOTHING_, and assume Dana knows what he's doing
          Case (BC_T_RADIATION)
             continue                   ! do _NOTHING_, and assume Dana knows what he's doing
          Case (BC_T_VFRADIATION)
             continue                   ! do _NOTHING_, and assume Dana knows what he's doing
          Case (BC_T_HTC)
             continue                   ! do _NOTHING_, and assume Dana knows what he's doing
          Case (BC_T_NEUMANN)
             continue                   ! do _NOTHING_, and assume Matt knows what he's doing
          Case Default
             Call TLS_fatal ('FACE_GRADIENT_T: bad BC type - only Dirichlet conditions implemented')
          End Select
       End Do

       ! calculate the gradient of phi by LSLR
       If (Present(WD)) Then
          Call LSLR_T (X_n, Y_n, Z_n, Phi_n,      &
                       X_e, Y_e, Z_e, Phi_e,      &
                       X_f, Y_f, Z_f, Phi_f,      &
                       dPhi_dx, dPhi_dy, dPhi_dz, &
                       Phi_c=Phi, face=face, WD=WD)
       Else
          Call LSLR_T (X_n, Y_n, Z_n, Phi_n,      &
                       X_e, Y_e, Z_e, Phi_e,      &
                       X_f, Y_f, Z_f, Phi_f,      &
                       dPhi_dx, dPhi_dy, dPhi_dz, &
                       Phi_c=Phi, face=face)
       End If

    Case Default
       Call TLS_fatal ('FACE_GRADIENT_T: invalid face gradient method specified')

    End Select

    Return

  End Subroutine FACE_GRADIENT_T

  !-----------------------------------------------------------------------------

  Subroutine LSLR_T (X_n, Y_n, Z_n, Phi_n,      &
                     X_e, Y_e, Z_e, Phi_e,      &
                     X,   Y,   Z,   Phi,        &
                     dPhi_dx, dPhi_dy, dPhi_dz, &
                     Phi_c, face, WD)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   *** special case - temperature only ***
    !
    !   Compute the least-squares linear reconstruction gradient of Phi
    !   at either cell centers or cell faces. (X_v,Y_v,Z_v) are the
    !   cell vertex coordinates, Phi_v are the vertex values of Phi,
    !   (X_e,Y_e,Z_e) are the cell neighbor centroid coordinates, Phi_e
    !   are the cell neighbor values of Phi, and (X,Y,Z) is the point where
    !   the gradient of Phi is desired.
    !
    !   can be used to calculate face gradients (FG) or cell gradients (CG)
    !---------------------------------------------------------------------------

    Use cutoffs_module,       Only: alittle
    Use kind_module,          Only: log_kind, int_kind, real_kind
    Use mesh_module,          Only: Cell, Face_Vrtx, Mesh, Vrtx_Face
    Use parameter_module,     Only: ncells, nfc, nnodes, nvc
    Use bc_module,            Only: BCMATCH, BC_T_DIRICHLET, BC_T
    use truchas_logging_services, only: TLS_fatal_if_any

    Implicit None

    ! arguments
    Real(real_kind),   Dimension(nvc,ncells),           Intent(IN)  :: X_n, Y_n, Z_n ! node coordinates
    Real(real_kind),   Dimension(nvc,ncells),           Intent(IN)  :: Phi_n         ! node values
    Real(real_kind),   Dimension(nfc,ncells),           Intent(IN)  :: X_e, Y_e, Z_e ! neighbor cell coordinates
    Real(real_kind),   Dimension(nfc,ncells),           Intent(IN)  :: Phi_e         ! neighbor cell values
    Real(real_kind),   Dimension(ncells),               Intent(IN)  :: X, Y, Z       ! desired point coordinates
    Real(real_kind),   Dimension(ncells),               Intent(IN)  :: Phi           ! desired point value
    Real(real_kind),   Dimension(ncells),     Optional, Intent(IN)  :: Phi_c         ! cell center value (not always needed)
    Integer(int_kind),                        Optional, Intent(IN)  :: face          ! face number (not always needed)
    Real(real_kind),   Dimension(ncells),     Optional, Intent(IN)  :: WD            ! weight data (not always wanted)
    Real(real_kind),   Dimension(ncells),               Intent(OUT) :: dPhi_dx       ! gradient - x component
    Real(real_kind),   Dimension(ncells),               Intent(OUT) :: dPhi_dy       ! gradient - y component
    Real(real_kind),   Dimension(ncells),               Intent(OUT) :: dPhi_dz       ! gradient - z component

    ! local variables
    Logical(log_kind), Dimension(ncells) :: Mask
    Integer(int_kind)                    :: f
    Integer(int_kind)                    :: v
    Integer(int_kind)                    :: c
    Integer                              :: status
    Real(real_kind), Dimension(ncells)       :: Axx
    Real(real_kind), Dimension(ncells)       :: Axy
    Real(real_kind), Dimension(ncells)       :: Axz
    Real(real_kind), Dimension(ncells)       :: Ayy
    Real(real_kind), Dimension(ncells)       :: Ayz
    Real(real_kind), Dimension(ncells)       :: Azz
    Real(real_kind), Dimension(ncells)       :: Bx
    Real(real_kind), Dimension(ncells)       :: By
    Real(real_kind), Dimension(ncells)       :: Bz
    Real(real_kind), Dimension(ncells)       :: dX
    Real(real_kind), Dimension(ncells)       :: dY
    Real(real_kind), Dimension(ncells)       :: dZ
    Real(real_kind), Dimension(ncells)       :: W
    Real(real_kind), Allocatable, Dimension(:)   :: WD_Vrtx
    Real(real_kind), Allocatable, Dimension(:,:) :: WD_Vrtx_Cell
    Real(real_kind), Allocatable, Dimension(:,:) :: WD_Nbr

    !---------------------------------------------------------------------------

    ! See if this is to be a face or cell-centered gradient
    ! Construct the normal equation matrix A and RHS vector B
    ! Add in those vertices that are not on a boundary face
    Axx = 0.0d0
    Axy = 0.0d0
    Axz = 0.0d0
    Ayy = 0.0d0
    Ayz = 0.0d0
    Azz = 0.0d0
    Bx  = 0.0d0
    By  = 0.0d0
    Bz  = 0.0d0

    If (Present(WD)) Then
       ! allocate some working arrays for the weights
       Allocate (WD_Vrtx(nnodes), STAT=status)
       Call TLS_fatal_if_any ((status /= 0), 'LSLR_T: error allocating WD_Vrtx(nnodes) array')
       Allocate (WD_Vrtx_Cell(nvc,ncells), STAT=status)
       Call TLS_fatal_if_any ((status /= 0), 'LSLR_T: error allocating WD_Vrtx_Cell(nvc,ncells) array')
       Allocate (WD_Nbr(nfc,ncells), STAT=status)
       Call TLS_fatal_if_any ((status /= 0), 'LSLR_T: error allocating WD_NBR(nfc,ncells) array')

       ! average the supplied weight data to the vertices
       Call VERTEX_AVG (WD, WD_Vrtx)
       ! gather the weight data from the neighbors
       Call EN_GATHER (WD_Vrtx_Cell, WD_Vrtx)
       Call EE_GATHER (WD_Nbr, WD)
    End If

    Do v = 1, nvc
       ! construct distance deltas
       dX =  X_n(v,:) - X
       dY =  Y_n(v,:) - Y
       dZ =  Z_n(v,:) - Z

       ! inverse distance weighting
       W = dX*dX + dY*dY + dZ*dZ
       If (Present(WD)) Then
          W = W * WD_Vrtx_Cell(v,:)
       End If
       W = Merge (0.0d0, 1.0d0/W, Abs(W) <= alittle)

       ! Mask off the boundary vertices for cell-centered gradients
       ! and those vertices not on this face for face-centered gradients
       If (Present(face)) Then
          Mask = v == Face_Vrtx(face,1) .Or. v == Face_Vrtx(face,2) .Or. &
                 v == Face_Vrtx(face,3) .Or. v == Face_Vrtx(face,4)
       Else
          ! note - this is broken for internal boundaries
          Mask = Mesh%Ngbr_cell(Vrtx_Face(v,1)) /= 0 .And. &
                 Mesh%Ngbr_cell(Vrtx_Face(v,2)) /= 0 .and. &
                 Mesh%Ngbr_cell(Vrtx_Face(v,3)) /= 0
       End If

       ! do the sums
       Where (Mask)
          Axx = Axx + dX*dX*W
          Axy = Axy + dX*dY*W
          Axz = Axz + dX*dZ*W
          Ayy = Ayy + dY*dY*W
          Ayz = Ayz + dY*dZ*W
          Azz = Azz + dZ*dZ*W
          Bx = Bx + (Phi_n(v,:) - Phi)*dX*W
          By = By + (Phi_n(v,:) - Phi)*dY*W
          Bz = Bz + (Phi_n(v,:) - Phi)*dZ*W
       End Where
    End Do

    ! add in neighboring elements
    Do f = 1, nfc
       ! if calculating a face gradient and f isn't the face of concern, continue
       If (Present(face)) Then
          If (f /= face) Then
             Cycle
          End If
       End If

       ! construct deltas
       dX  = X_e(f,:) - X
       dY  = Y_e(f,:) - Y
       dZ  = Z_e(f,:) - Z
       W = dX*dX + dY*dY + dZ*dZ
       if (PRESENT(WD)) W = W * WD_Nbr(f,:)
       W = MERGE (0.0d0, 1.0d0/W, ABS(W) <= alittle)

       ! mask off boundaries
       ! note - this is broken for internal boundaries
       Mask = Mesh%Ngbr_cell(f) /= 0

       ! if this is a cell-centered gradient, correct the mask for Dirichlet boundary faces
       If (.Not. Present(face)) Then
          Do c = 1, ncells
             mask(c) = mask(c) .Or. BCMATCH(BC_T, f, c) == BC_T_DIRICHLET
          End Do
       End If

       ! do the sums
       Where (Mask)
          Axx = Axx + dX*dX*W
          Axy = Axy + dX*dY*W
          Axz = Axz + dX*dZ*W
          Ayy = Ayy + dY*dY*W
          Ayz = Ayz + dY*dZ*W
          Azz = Azz + dZ*dZ*W
          Bx = Bx + (Phi_e(f,:) - Phi)*dX*W
          By = By + (Phi_e(f,:) - Phi)*dY*W
          Bz = Bz + (Phi_e(f,:) - Phi)*dZ*W
       End Where
    End Do

    ! add in cell itself if this is a face gradient
    If (Present(face)) Then
       ! construct deltas
       dX =  Cell%Centroid(1) - X
       dY =  Cell%Centroid(2) - Y
       dZ =  Cell%Centroid(3) - Z

       W = dX*dX + dY*dY + dZ*dZ
       If (Present(WD)) W = W * WD
       W = MERGE (0.0d0, 1.0d0/W, ABS(W) <= alittle)

       ! do the sums
       Axx = Axx + dX*dX*W
       Axy = Axy + dX*dY*W
       Axz = Axz + dX*dZ*W
       Ayy = Ayy + dY*dY*W
       Ayz = Ayz + dY*dZ*W
       Azz = Azz + dZ*dZ*W
       Bx = Bx + (Phi_c - Phi)*dX*W
       By = By + (Phi_c - Phi)*dY*W
       Bz = Bz + (Phi_c - Phi)*dZ*W
    End If

    ! solve the 3X3 system with Cramer's rule
    ! denominator (determinant of coefficient matrix)
    dX = Axx*(Ayy*Azz - Ayz*Ayz) + &
         Axy*(Axz*Ayz - Axy*Azz) + &
         Axz*(Axy*Ayz - Axz*Ayy)
    dX = MERGE (0.0d0, 1.0d0/dX, ABS(dX) <= alittle)

    ! numerators (determinant of matrix with B inserted)
    dPhi_dx =  Bx*(Ayy*Azz - Ayz*Ayz) + &
              Axy*(Bz *Ayz - By *Azz) + &
              Axz*(By *Ayz - Bz *Ayy)
    dPhi_dy = Axx*(By *Azz - Bz *Ayz) + &
               Bx*(Axz*Ayz - Axy*Azz) + &
              Axz*(Bz *Axy - By *Axz)
    dPhi_dz = Axx*(Bz *Ayy - By *Ayz) + &
              Axy*(By *Axz - Bz *Axy) + &
               Bx*(Axy*Ayz - Axz*Ayy)

    ! multiply by inverse determinant
    dPhi_dx = dX*dPhi_dx
    dPhi_dy = dX*dPhi_dy
    dPhi_dz = dX*dPhi_dz

    ! free allocated space
    If (Present(WD)) Then
       Deallocate (WD_Vrtx)
       Deallocate (WD_Vrtx_Cell)
       Deallocate (WD_Nbr)
    End If

    Return

  End Subroutine LSLR_T


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
    use constants_module,     only: one, zero
    use cutoffs_module,       only: alittle
    use kind_module,          only: int_kind, real_kind
    use linear_module,        only: LINEAR_PROP
!   use matl_module,          only: GATHER_VOF,MATL
!   use mesh_module,          only: Cell, Mesh, Face_Vrtx, Vertex, Vrtx_Bdy
    use mesh_module,          only: Cell
    use parameter_module,     only: ncells, nfc, nnodes, nvc, ndim

    implicit none

    ! Arguments
    real(KIND = real_kind),   dimension(ndim,ncells), intent(OUT)          :: Gradient_Phi
    real(KIND = real_kind),   dimension(ncells),      intent(IN)           :: Phi

    ! Local Variables
    integer(KIND = int_kind) :: i, f, n

!   real(KIND = real_kind),   dimension(ncells)             :: Phi_f, Weight, Vof
    real(KIND = real_kind),   dimension(ncells)             :: Phi_f
    real(KIND = real_kind)                                  :: Ptmp
    real(KIND = real_kind)                                  :: Gtmp
    real(KIND = real_kind),   dimension(nvc,ncells)         :: Phi_v
    real(KIND = real_kind),   dimension(nnodes)             :: Phi_vtx
!   real(KIND = real_kind),   dimension(:,:,:),     pointer :: X_v

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Compute Phi_vtx, a vertex value of Phi
    call VERTEX_AVG (Phi, Phi_vtx)

    ! Gather vertex Phi_vtx into cell vector Phi_v
    call EN_GATHER (Phi_v, Phi_vtx)

    ! compute gradient according to Green-Gauss method
    ! Loop over faces, accumulating the product
    ! Phi_f*Face_Normal for each area vector component

    Gradient_Phi = zero
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
    Ptmp = one/Cell(i)%Volume
    do n = 1,ndim
       Gtmp = Gradient_Phi(n,i)*Ptmp
       ! Eliminate noise
       if(ABS(Gtmp) <= alittle)then
         Gradient_Phi(n,i) = zero
       else
         Gradient_Phi(n,i) = Gtmp
       endif
    end do
    end do
    return
  END SUBROUTINE GRADIENT_CELL

  SUBROUTINE GRADIENT_CELL_FROM_FACE_VALUES (Gradient_Phi, Phi_f)
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
    use constants_module,     only: one, zero
    use cutoffs_module,       only: alittle
    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Cell
    use parameter_module,     only: ncells, nfc, ndim

    implicit none

    ! Arguments
    real(KIND = real_kind),   dimension(ndim,ncells), intent(OUT)          :: Gradient_Phi
    real(KIND = real_kind),   dimension(nfc,ncells),  intent(IN)           :: Phi_f

    ! Local Variables
    integer(KIND = int_kind) :: i, f, n

    real(KIND = real_kind)                                  :: Ptmp
    real(KIND = real_kind)                                  :: Gtmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! compute gradient according to Green-Gauss method
    ! Loop over faces, accumulating the product
    ! Phi_f*Face_Normal for each area vector component

    Gradient_Phi = zero
    do f = 1,nfc
       ! Accumulate the dot product
       do n = 1,ndim  
          Gradient_Phi(n,:) = Gradient_Phi(n,:) + Cell%Face_Normal(n,f)*Cell%Face_Area(f)*Phi_f(f,:)
       end do
    end do
    ! Normalize by cell volume
    do i = 1,ncells
    Ptmp = one/Cell(i)%Volume
    do n = 1,ndim
       Gtmp = Gradient_Phi(n,i)*Ptmp
       ! Eliminate noise
       if(ABS(Gtmp) <= alittle)then
         Gradient_Phi(n,i) = zero
       else
         Gradient_Phi(n,i) = Gtmp
       endif
    end do
    end do
    return
  END SUBROUTINE GRADIENT_CELL_FROM_FACE_VALUES

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
    use constants_module,     only: zero
    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Cell, Mesh
    use parameter_module,     only: ncells, nfc, ndim

    use gs_module,            only: EE_GATHER 

    implicit none

    ! Arguments
    real(KIND = real_kind),   dimension(ndim,nfc,ncells), intent(OUT)  :: facegrad
    real(KIND = real_kind),   dimension(ncells),  intent(IN)           :: pc
    real(KIND = real_kind),   dimension(nfc,ncells),  intent(IN)       :: ghc
    real(KIND = real_kind),   dimension(nfc,ncells),  intent(IN)       :: ghn

    ! Local Variables
    integer(KIND = int_kind) :: i, f, n, d

    real(KIND = real_kind)                             :: dl2, dl2i, plr1, plr2

    real(KIND = real_kind), dimension(nfc,ncells)      :: pn
    real(KIND = real_kind), dimension(ndim)            :: dx, cn
    real(KIND = real_kind), dimension(ndim,nfc,ncells) :: Cell_Ngbr_Coord

    !***************************************************************************

    facegrad = zero
    
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
          
          dl2 = zero
          do d=1,ndim
             dl2 = dl2 + dx(d)*dx(d)
          end do 
          dl2i = 1.0/dl2
    
          facegrad(1,f,n) = (plr1 - plr2)*dx(1)*dl2i
          facegrad(2,f,n) = (plr1 - plr2)*dx(2)*dl2i
          facegrad(3,f,n) = (plr1 - plr2)*dx(3)*dl2i
          
          if (Mesh(n)%Ngbr_Cell(f) == 0) then
             facegrad(:,f,n) = zero
          end if
          
       end do
    end do
    
    return

  END SUBROUTINE DYNAMIC_PRESSURE_FACE_GRADIENT

  !-----------------------------------------------------------------------------

END MODULE DISCRETE_OP_MODULE
