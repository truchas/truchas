!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NODE_OP_SETUP_MODULE
  !=============================================================================
  ! Purpose:
  !
  !  Allocate data structures and set up control volume geometries, integration
  !  points and jacobians for a node based operator.
  !
  !  Public procedures:
  !       call ALLOCATE_CONTROL_VOLUME()
  !       call CELL_CV_FACE()
  !
  ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
  !
  !=============================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: ALLOCATE_CONTROL_VOLUME, CELL_CV_FACE, BOUNDARY_CV_FACE
  ! Private variables
  logical, pointer, dimension(:,:)  :: Node_Mask
  integer, pointer, dimension(:)    :: GNode
  real(r8), pointer, dimension(:,:)    :: Node_Value
  real(r8), parameter           :: small_angle = 0.02 ! Largest angle treated as a corner

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>
  SUBROUTINE ALLOCATE_CONTROL_VOLUME()
    !=============================================================================
    !
    !  Allocate data structures for internal control volume faces
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use node_operator_module,     only: CV_Internal, Nodal_Volume, nipc
    use legacy_mesh_api,          only: ncells, nnodes
    use solid_mechanics_mesh,     only: ndim
    use solid_mechanics_input,    only: stress_reduced_integration

    ! Local variables
    integer :: status
    ! 
    ! Control volume faces internal to the cells
    allocate(CV_Internal%Face_Normal(ndim,nipc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Internal%Face_Normal')
    allocate(CV_Internal%Face_Area(nipc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Internal%Face_Area')
    allocate(Nodal_Volume(nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Nodal_Volume')
    allocate(CV_Internal%Face_Coord(ndim,nipc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Internal%Face_Coord')
    if (.not. stress_reduced_integration) then
       allocate(CV_Internal%Face_Ijac(ndim,ndim,nipc,ncells), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Internal%Face_Ijac')
    end if

    call MECH_BC_ALLOCATE

  END SUBROUTINE ALLOCATE_CONTROL_VOLUME

  !---------------------------------------------------------------------------------
  SUBROUTINE MECH_BC_ALLOCATE
    !=============================================================================
    !
    !  Get BC data, find connectivity for gap nodes and allocate data structures 
    !  for boundary control volume faces.
    !
    ! This routine collect BC data to set up arrays used to transfer data from the
    ! BC atlas structures to the solid mech BC structures.
    !
    !  Node_Mask lists which nodes are included in each BC_OPERATOR
    !  Node_Value collects the BC values for each node for each operator
    !  GNode is allocated and used later to collect gap neighbor info
    !
    ! The nbface and nbnode arrays list the number of nodes and faces for each BC
    ! operator
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use node_operator_module,     only: nipbf, CV_Boundary, nbface, nbnode
    use bc_data_types
    use mech_bc_data_module
    use legacy_mesh_api,          only: ncells, nnodes, EN_OR_Scatter, EN_MAX_Scatter
    use solid_mechanics_mesh,     only: ndim, nvc
    use var_vector_module

    ! Local variables
    integer :: status, i, j, atlas_size, c, inode, nbnodes_tot
    logical, pointer, dimension(:,:) :: Mask_EN
    real(r8), pointer, dimension(:,:) :: Value_EN
    integer, dimension(nipbf) :: v
    type (BC_Operator), POINTER :: Operator
    type (BC_Atlas), POINTER :: Atlas
    integer, pointer, dimension(:) :: Cell_List
    integer, pointer, dimension(:) :: Face_List
    real(r8), pointer, dimension(:,:) :: Value_List

    ! Control volume faces on boundaries  - using "new" bc stuff
    ! Number of bc atlases with data
    allocate(Node_Mask(BC_MAX_OPERATORS,nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Mask')
    allocate(Node_Value(BC_MAX_OPERATORS,nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Value')
    allocate(GNode(nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: GNode')
    allocate(Mask_EN(nvc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Mask_EN')
    allocate(Value_EN(nvc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Value_EN')
    Node_Mask = .false.
    Node_Value = -1.0e10
    GNode = 0
    nbface = 0.0
    nbnode = 0.0

    ! Loop over all valid bc operators.  Only operate on the solid mech BCs applied at nodes.
    BC_OP_LOOP: do i = 1, BC_MAX_OPERATORS
       select case (i)
          ! Traction BCs
       case(BC_X_TRACTION_OP, BC_Y_TRACTION_OP, BC_Z_TRACTION_OP, &
            BC_NORM_TRACTION_OP)
          Operator => BC_Spec_Get_Operator(Displacement_BC,i)
          Atlas => BC_Op_Get_Atlas(Operator)
          atlas_size = DATA_SIZE(Atlas)
          nbface(i) = atlas_size
          ! Displacement BCs
       case(BC_X_DISPLACEMENT_OP, BC_Y_DISPLACEMENT_OP, BC_Z_DISPLACEMENT_OP, &
            BC_NORM_DISPLACEMENT_OP, BC_NORMAL_CONSTRAINT_OP, BC_CONTACT_OP,&
            BC_FREE_INTERFACE_OP)
          Operator => BC_Spec_Get_Operator(Displacement_BC,i)
          Atlas => BC_Op_Get_Atlas(Operator)
          atlas_size = DATA_SIZE(Atlas)
          nbface(i) = atlas_size
          Cell_List => BC_Get_Cell(Atlas)
          Face_List => BC_Get_Face(Atlas)
          Value_List => BC_Get_Values(Atlas)
          Mask_EN = .false.
          Value_EN = -1.0e10
          do c = 1,atlas_size
             select case (Face_List(c))
             case (1)
                v(1)=3; v(2)=4; v(3)=8; v(4)=7
             case (2)
                v(1)=1; v(2)=2; v(3)=6; v(4)=5
             case (3)
                v(1)=4; v(2)=1; v(3)=5; v(4)=8
             case (4)
                v(1)=2; v(2)=3; v(3)=7; v(4)=6
             case (5)
                v(1)=1; v(2)=2; v(3)=3; v(4)=4
             case (6)
                v(1)=7; v(2)=8; v(3)=5; v(4)=6
             end select
             do j =1, nipbf
                Mask_EN(v(j),Cell_List(c)) = .true.
                Value_EN(v(j),Cell_List(c)) = Value_List(1,c)
             end do
          end do
          call EN_OR_Scatter(Node_Mask(i,:),Mask_EN)
          ! Take the most positive value for each node.  This may not be what is
          ! wanted, but if the input specifies two values (node is on a corner between 
          ! two surfaces) then the user should overwrite with individual nodal BCs.
          call EN_MAX_Scatter(Node_Value(i,:),Value_EN)
          nbnode(i) = COUNT(Node_Mask(i,:))
       end select
    end do BC_OP_LOOP

    ! Add individual nodal displacement constraints.  Currently it is not possible to specify normal 
    ! constraints or contact at a single node.  Node_Disp_BC_Temp only has data in the first index of 
    ! its components.

    do i = 1,SIZE(Node_Disp_BC_Temp%Node)
       inode = Node_Disp_BC_Temp%Node(i)
       select case(Node_Disp_BC_Temp%BC_Type(i,1))
       case(X_DISPLACEMENT)
          if (Node_Mask(BC_X_DISPLACEMENT_OP,inode)) then
             ! issue warning here
          end if
          Node_Mask(BC_X_DISPLACEMENT_OP,inode) = .true.
          Node_Value(BC_X_DISPLACEMENT_OP,inode) = Node_Disp_BC_Temp%Value(i,1)
          nbnode(BC_X_DISPLACEMENT_OP) = COUNT(Node_Mask(BC_X_DISPLACEMENT_OP,:))
       case(Y_DISPLACEMENT)
          if (Node_Mask(BC_Y_DISPLACEMENT_OP,inode)) then
             ! issue warning here
          end if
          Node_Mask(BC_Y_DISPLACEMENT_OP,inode) = .true.
          Node_Value(BC_Y_DISPLACEMENT_OP,inode) = Node_Disp_BC_Temp%Value(i,1)
          nbnode(BC_Y_DISPLACEMENT_OP) = COUNT(Node_Mask(BC_Y_DISPLACEMENT_OP,:))
       case(Z_DISPLACEMENT)
          if (Node_Mask(BC_Z_DISPLACEMENT_OP,inode)) then
             ! issue warning here
          end if
          Node_Mask(BC_Z_DISPLACEMENT_OP,inode) = .true.
          Node_Value(BC_Z_DISPLACEMENT_OP,inode) = Node_Disp_BC_Temp%Value(i,1)
          nbnode(BC_Z_DISPLACEMENT_OP) = COUNT(Node_Mask(BC_Z_DISPLACEMENT_OP,:))
       case(NORMAL_DISPLACEMENT)
          if (Node_Mask(BC_NORM_DISPLACEMENT_OP,inode)) then
             ! issue warning here
          end if
          Node_Mask(BC_NORM_DISPLACEMENT_OP,inode) = .true.
          Node_Value(BC_NORM_DISPLACEMENT_OP,inode) = Node_Disp_BC_Temp%Value(i,1)
          nbnode(BC_NORM_DISPLACEMENT_OP) = COUNT(Node_Mask(BC_NORM_DISPLACEMENT_OP,:))
       case(CONTACT)
          if (Node_Mask(BC_CONTACT_OP,inode)) then
             ! issue warning here
          end if
          Node_Mask(BC_CONTACT_OP,inode) = .true.
          Node_Value(BC_CONTACT_OP,inode) = Node_Disp_BC_Temp%Value(i,1)
          nbnode(BC_CONTACT_OP) = COUNT(Node_Mask(BC_CONTACT_OP,:))
       end select
    end do

    nbnodes_tot=0
    do i = 1,nnodes
       if (ANY(Node_Mask(:,i))) nbnodes_tot = nbnodes_tot + 1
    end do

    NULLIFY (Operator, Atlas, Cell_List, Face_List)

    ! Allocate for boundary
    allocate(CV_Boundary(BC_MAX_OPERATORS), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary')
    do i=1,BC_MAX_OPERATORS
       ! Cell numbers
       allocate(CV_Boundary(i)%Cell(nbface(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Cell')
       ! Node numbers (array size will be zero for traction BCs)
       allocate(CV_Boundary(i)%Node(nbnode(i),1+ndim), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Node')
!       ! Interface ID for each of up to ndim normals, values, etc.
!       allocate(CV_Boundary(i)%Interface(ndim,nbnode(i)), stat=status)
!       if (status /= 0) call punt((/'allocation error: CV_Boundary%Interface'/), 'ALLOCATE_CONTROL_VOLUME')
       ! BC nodal values (displacements or friction coefficients or tangential tractions or ??)
       allocate(CV_Boundary(i)%Node_Value(nbnode(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Node_Value')
       ! Vector normal to the surface at the node, pointing towards the interface or free surface 
       allocate(CV_Boundary(i)%Node_Normal(ndim,nbnode(i),ndim), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Node_Normal')
       ! Number of interfaces (and normals) associated with each node
       allocate(CV_Boundary(i)%NN_Count(nbnode(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%NN_Count')
       ! Cell face numbers 
       allocate(CV_Boundary(i)%Face(nbface(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Face')
       ! Node numbers for each cell face - these also correspond to individual subfaces on the cell face
       allocate(CV_Boundary(i)%Face_Node(nipbf,nbface(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Face_Node')
       ! Subface normals at subface centroids, pointing outward from the cell
       allocate(CV_Boundary(i)%Face_Normal(ndim,nipbf,nbface(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Face_Normal')
       ! Subface areas
       allocate(CV_Boundary(i)%Face_Area(nipbf,nbface(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Face_Area')
       ! Subface centroid coordinates
       allocate(CV_Boundary(i)%Face_Coord(ndim,nipbf,nbface(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Face_Coord')
       ! Inverse jacobians at subface centroids
       allocate(CV_Boundary(i)%Face_Ijac(ndim,ndim,nipbf,nbface(i)), stat=status)
       if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: CV_Boundary%Face_Ijac')
       ! Face gap displacement

       ! Initialize to known values
       CV_Boundary(i)%Cell = 0
       CV_Boundary(i)%Node = 0
       CV_Boundary(i)%Node_Value = 0.0
       CV_Boundary(i)%Node_Normal = 0.0
       CV_Boundary(i)%NN_Count = 0
       CV_Boundary(i)%Face = 0
       CV_Boundary(i)%Face_Node = 0
       CV_Boundary(i)%Face_Normal = 0.0
       CV_Boundary(i)%Face_Area = 0.0
       CV_Boundary(i)%Face_Coord = -1.0e10
       CV_Boundary(i)%Face_Ijac = 0.0
    end do

    ! Allocate nodal BC structure
    allocate(Node_Displacement_BC%Node(nbnodes_tot), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Node')
    allocate(Node_Displacement_BC%Gap_Node(nbnodes_tot,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Gap_Node')
    allocate(Node_Displacement_BC%BC_Type(nbnodes_tot,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%BC_Type')
    allocate(Node_Displacement_BC%Interface(ndim,nbnodes_tot), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Interface')
    allocate(Node_Displacement_BC%Normal(nbnodes_tot,ndim,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Normal')
    allocate(Node_Displacement_BC%Combination(nbnodes_tot), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Combination')
    allocate(Node_Displacement_BC%Value(nbnodes_tot,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Value')
    allocate(Node_Displacement_BC%Vector(nbnodes_tot,ndim,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Vector')
    allocate(Node_Displacement_BC%Lambda(nbnodes_tot,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Lambda')
    allocate(Node_Displacement_BC%Scalar(nbnodes_tot,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Scalar')
    allocate(Node_Displacement_BC%Gap_Disp(nbnodes_tot,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Gap_Disp')
    allocate(Node_Displacement_BC%Normal_Traction(nbnodes_tot,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Normal_Traction')
    allocate(Node_Displacement_BC%Area(nbnodes_tot,ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: Node_Displacement_BC%Area')

    Node_Displacement_BC%Node = 0
    Node_Displacement_BC%Gap_Node = 0
    Node_Displacement_BC%BC_Type = 0
    Node_Displacement_BC%Interface = 0
    Node_Displacement_BC%Normal = 0.0
    Node_Displacement_BC%Combination = TRACTION_ONLY
    Node_Displacement_BC%Value = 0.0
    Node_Displacement_BC%Vector = 0.0
    Node_Displacement_BC%Lambda = 0.0
    Node_Displacement_BC%Scalar = 0.0
    Node_Displacement_BC%Gap_Disp = 0.0
    Node_Displacement_BC%Normal_Traction = 0.0
    Node_Displacement_BC%Area = 0.0

    ! Allocate array to relate mesh node numbers to position in Node_Displacement_BC
    allocate(NBC_index(nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'ALLOCATE_CONTROL_VOLUME: allocation error: NBC_index')

    deallocate(Mask_EN)
    deallocate(Value_EN)

    return
  END SUBROUTINE MECH_BC_ALLOCATE

  !---------------------------------------------------------------------------------
  SUBROUTINE CELL_CV_FACE()
    !
    !  Calculate all internal control volume face areas, and normals and inverse 
    !  jacobians for a cell
    !
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use legacy_mesh_api, only: ncells, Mesh, gather_vertex_coord, EN_SUM_Scatter
    use legacy_mesh_api, only: CELL_TET, CELL_PYRAMID, CELL_PRISM, CELL_EDGE, GAP_ELEMENT_1
    use cutoffs_module,        only: alittle
    use lu_solve_module,       only: LU_SOLVE, factsolve, backsolve
    use node_operator_module,  only: CV_Internal, nipc, Nodal_Volume, LINEAR_PROP, LINEAR_GRAD
    use solid_mechanics_input, only: stress_reduced_integration
    use solid_mechanics_mesh,  only: ndim, nvc, nec, nvf

    ! Local Variables
    integer :: status
    real(r8), pointer, dimension(:,:,:) :: Xv
    integer :: edge, i, j, v1, v2, v11, v12, v13, v14, v21, v22, v23, v24, &
         v31, v32, v33, v34, n, nn, i1, i2, solve_flag
    integer, Dimension(ndim)            :: indx
    real(r8), Dimension(:,:,:),pointer  :: Xn
    real(r8), Dimension(nvc,ndim)       :: Xlc
    real(r8), Dimension(:,:,:),pointer  :: Xl
    real(r8), Dimension(ndim)           :: Coef
    real(r8), Dimension(:,:), pointer   :: X1, X2
    real(r8), Dimension(:,:), pointer   :: cell_cen_l
    real(r8)                            :: area_l, nsign
    real(r8), Dimension(:,:,:), pointer :: ip_crd_lp
    real(r8), Dimension(:,:), pointer   :: Cell_Vertex_Volume
    real(r8), Dimension(:,:), pointer   :: Projection
    real(r8), Dimension(:,:,:),pointer  :: Xlt
    real(r8), Dimension(:,:,:),pointer  :: Xg
    ! jac only need to be allocated if we are using full integration
    real(r8), Dimension(:,:,:,:), pointer  :: jac
    logical                               :: pivot = .true.

    ! Inform the user of control volume face calculation.
    call TLS_info ('')
    call TLS_info (' Finding control volume faces internal to cells ... ')
    ! Define vertex coordinates for logical cell
    Xlc(1,1) = 1.0_r8; Xlc(1,2) = 0.0_r8; Xlc(1,3) = 0.0_r8
    Xlc(2,1) = 1.0_r8; Xlc(2,2) = 1.0_r8; Xlc(2,3) = 0.0_r8
    Xlc(3,1) = 0.0_r8; Xlc(3,2) = 1.0_r8; Xlc(3,3) = 0.0_r8
    Xlc(4,1) = 0.0_r8; Xlc(4,2) = 0.0_r8; Xlc(4,3) = 0.0_r8
    Xlc(5,1) = 1.0_r8; Xlc(5,2) = 0.0_r8; Xlc(5,3) = 1.0_r8
    Xlc(6,1) = 1.0_r8; Xlc(6,2) = 1.0_r8; Xlc(6,3) = 1.0_r8
    Xlc(7,1) = 0.0_r8; Xlc(7,2) = 1.0_r8; Xlc(7,3) = 1.0_r8
    Xlc(8,1) = 0.0_r8; Xlc(8,2) = 0.0_r8; Xlc(8,3) = 1.0_r8
    !
    ! Create temporary arrays
    allocate(cell_cen_l(ndim,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: cell_cen_l')
    allocate(ip_crd_lp(ndim,nipc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: ip_crd_lp')
    allocate(Xl(nvc,ndim,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: Xl')
    allocate(Xg(ndim,nipc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: Xg')
    allocate(Xn(ndim,nvf,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: Xn')
    allocate(X1(ndim,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: X1')
    allocate(X2(ndim,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: X2')
    allocate(Xv(ndim,nvc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: Xv')
    allocate(Cell_Vertex_Volume(nvc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: Cell_Vertex_Volume')
    allocate(Projection(2,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: Projection')
    ! Xlt is used for more than one purpose as a temporary array
    allocate(Xlt(nvc,nipc,ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: Xlt')
    if (.not. stress_reduced_integration) then
       allocate(jac(ndim,ndim,nipc,ncells), stat=status)
       if (status /= 0) call TLS_panic ( 'CELL_CV_FACES: allocation error: jac')
    end if
    ! Gather vertex coordinates
    call gather_vertex_coord (Xv)
    ! Get logical centroid
    call CELL_LOGICAL_CENTROID(Xv,cell_cen_l)
    ! write(*,*) ' Cell logical centroid', cell_cen_l(1),cell_cen_l(2),cell_cen_l(3)

    
    ! Loop over edges - each edge corresponds to a control volume face which is used for the
    ! two nodes on that edge
    EDGE_LOOP: do edge=1,nec

       call CV_FACE_COORDS(edge, ndim, Xlc, Xv, cell_cen_l, nn, Xl, Xn, Coef, nsign,  ip_crd_lp,&
                           v11, v12, v13, v14, v21, v22, v23, v24, v31, v32, v33, v34)


       ! Calculate CV face area and normal
       do i = 1,ndim
          X1(i,:) = Xn(i,1,:) - Xn(i,2,:)
          X2(i,:) = Xn(i,3,:) - Xn(i,4,:)
       end do
       ! These are face area vectors for now.
       do i = 1,ndim
          select case (i)
          case (1)
             v1 = 2; v2 = 3
             if (ndim == 2) CV_Internal%Face_Normal(i,edge,:) = X1(v1,:)
          case (2)
             v1 = 3; v2 = 1
             if (ndim == 2) CV_Internal%Face_Normal(i,edge,:) = X2(v2,:)
          case (3)
             v1 = 1; v2 = 2
          end select
          if (ndim == 3) CV_Internal%Face_Normal(i,edge,:) = 0.5_r8*(X1(v1,:)*X2(v2,:) - X2(v1,:)*X1(v2,:))
       end do
       ! Set components to zero if they're small; accumulate areas.
       CV_Internal%Face_Area(edge,:) = 0.0_r8
       do i = 1,ndim
          CV_Internal%Face_Normal(i,edge,:) = MERGE(0.0_r8, CV_Internal%Face_Normal(i,edge,:), &
               ABS(CV_Internal%Face_Normal(i,edge,:)) < alittle)
          CV_Internal%Face_Area(edge,:) = CV_Internal%Face_Area(edge,:) + CV_Internal%Face_Normal(i,edge,:)**2
       end do
       ! Compute areas.
!       if (PGSLib_GLOBAL_ANY(CV_Internal%Face_Area(edge,:) == 0.0)) &
!            call PUNT((/'FATAL: Internal control volume faces with zero area found'/), 'CELL_CV_FACES')

       CV_Internal%Face_Area(edge,:) = SQRT(CV_Internal%Face_Area(edge,:))
       ! Convert the face vectors to unit normals.
       do i = 1,ndim
          where (CV_Internal%Face_Area(edge,:) >= alittle**2)
             CV_Internal%Face_Normal(i,edge,:) = CV_Internal%Face_Normal(i,edge,:) / CV_Internal%Face_Area(edge,:)
          elsewhere
             CV_Internal%Face_Normal(i,edge,:) = 0.0_r8
          end where
       end do
       ! Find the logical centroid of the CV face this is in a new (xi', eta', zeta') logical
       ! coordinate system.
       ! If we are using the one-point operator this is not needed, except for the preconditioner
       ! volumes.
       !
       ! Need to fix
       ! if (ndim == 2) cycle FACE_LOOP
       !
       ! Since this is in logical coordinates for the original cell, initial area=(1/12)/(1/4)
       ! Face_Area = one_twelfth / (CV_Internal(:)%Face_Area + alittle)
       area_l = 1.0_r8/3.0_r8
       NDIM_LOOP: do i = 1,ndim
          select case (i)
          case (1)
             i1 = 2; i2 = 3
          case (2)
             i1 = 3; i2 = 1
          case (3)
             i1 = 1; i2 = 2
          end select
          ip_crd_lp(1,edge,:) = ip_crd_lp(1,edge,:) + &
               Coef(1)*area_l*nsign*CV_Internal%Face_Normal(i,edge,:)* &
               ((Xn(i1,v11,:) - Xn(i1,v12,:))*(Xn(i2,v13,:) - Xn(i2,v14,:)) &
               - (Xn(i2,v11,:) - Xn(i2,v12,:))*(Xn(i1,v13,:) - Xn(i1,v14,:)))
          ip_crd_lp(2,edge,:) = ip_crd_lp(2,edge,:) + &
               Coef(2)*area_l*nsign*CV_Internal%Face_Normal(i,edge,:)* &
               ((Xn(i1,v21,:) - Xn(i1,v22,:))*(Xn(i2,v23,:) - Xn(i2,v24,:)) &
               - (Xn(i2,v21,:) - Xn(i2,v22,:))*(Xn(i1,v23,:) - Xn(i1,v24,:)))
          ip_crd_lp(3,edge,:) = ip_crd_lp(3,edge,:) + &
               Coef(3)*area_l*nsign*CV_Internal%Face_Normal(i,edge,:)* &
               ((Xn(i1,v31,:) - Xn(i1,v32,:))*(Xn(i2,v33,:) - Xn(i2,v34,:)) &
               - (Xn(i2,v31,:) - Xn(i2,v32,:))*(Xn(i1,v33,:) - Xn(i1,v34,:)))
       end do NDIM_LOOP
       ! Convert logical prime coordinates to logical coordinates for the original cell
       do i = 1,ndim
          Xlt(:,edge,:)=Xl(:,i,:)
          call LINEAR_PROP (ncells,ip_crd_lp(:,edge,:),Xlt(:,edge,:),CV_Internal%Face_Coord(i,edge,:))
       end do
       ! Convert logical coordinates of the CV face centroid to global coordinates for the 
       !nodal volume calculation.  Temporary array Xg is used; global coordinates are not saved.
       do i = 1,ndim
          Xlt(:,edge,:)=Xv(i,:,:)
          call LINEAR_PROP (ncells,CV_Internal%Face_Coord(:,edge,:),Xlt(:,edge,:),Xg(i,edge,:))
       end do
       !
    end do EDGE_LOOP
    !
    ! Set areas of missing edges of degenerate cells to zero.  (Leave faces for redundant edges.)
    where (Mesh%Cell_Shape <= CELL_PRISM)
       CV_Internal%Face_Area(10,:) = 0.0_r8
       CV_Internal%Face_Area(12,:) = 0.0_r8
    end where
    where (Mesh%Cell_Shape <= CELL_PYRAMID)
       CV_Internal%Face_Area(9,:) = 0.0_r8
       CV_Internal%Face_Area(11,:) = 0.0_r8
    end where
    where (Mesh%Cell_Shape == CELL_TET)
       CV_Internal%Face_Area(1,:) = 0.0_r8
    end where
    !
    ! Calculate volume of the control volumes for nodes:
    ! Accumulate volumes for tetrahedrons formed by the control volume faces and the two nodes associated
    ! with each face.  These will be sum scattered to the nodes.
    ! Initialize nodal volume array
    Cell_Vertex_Volume(:,:) = 0.0_r8
    VOLUME_EDGE_LOOP: do edge = 1,nec
       Projection(:,:) = 0.0_r8
       do i = 1,ndim
          Projection(1,:) = Projection(1,:) + (Xg(i,edge,:)-Xv(i,Cell_Edge(1,edge),:)) * &
               CV_Internal%Face_Normal(i,edge,:)
          Projection(2,:) = Projection(2,:) - (Xg(i,edge,:)-Xv(i,Cell_Edge(2,edge),:)) * &
               CV_Internal%Face_Normal(i,edge,:)
       end do
       Cell_Vertex_Volume(Cell_Edge(1,edge),:) = Cell_Vertex_Volume(Cell_Edge(1,edge),:) + Projection(1,:) * &
            CV_Internal%Face_Area(edge,:) / 3.0_r8
       Cell_Vertex_Volume(Cell_Edge(2,edge),:) = Cell_Vertex_Volume(Cell_Edge(2,edge),:) + Projection(2,:) * &
            CV_Internal%Face_Area(edge,:) / 3.0_r8
    end do VOLUME_EDGE_LOOP
    call EN_SUM_Scatter(Nodal_Volume,Cell_Vertex_Volume)
    !

    call TLS_info ('')
    call TLS_info (' Finding jacobians ... ')
    do n=1,ncells
       ! How are other derived types initialized??
       CV_Internal%Face_Ijac(:,:,:,n)=0.0_r8
       do i = 1,ndim
          do edge = 1,nipc
             Xlt(:,edge,n)=Xv(i,:,n)
             do j=1,ndim
                CV_Internal%Face_Ijac(j,j,edge,n)=1.0_r8
             end do
          end do
          ! Leave the jacobians for gap elements set to the identity matrix
          if (Mesh(n)%Cell_Shape < GAP_ELEMENT_1)  &
               call LINEAR_GRAD(nipc,CV_Internal%Face_Coord(:,:,n),Xlt(:,:,n),jac(i,:,:,n))
       end do
    end do
    call TLS_info ('')
    call TLS_info (' Inverting jacobians ... ')
    do n = 1,ncells
       ! Leave gap elements set to [I]
       if (Mesh(n)%Cell_Shape < GAP_ELEMENT_1) then
          do edge=1,nipc
             do i=1,ndim
                if(i==1) then
                   solve_flag=factsolve
                else
                   solve_flag=backsolve
                end if
                call LU_SOLVE (jac(:,:,edge,n),CV_Internal%Face_Ijac(:,i,edge,n),indx,solve_flag,pivot)
             end do
          end do
       end if
    end do
    !
    ! Destroy temporary arrays
    deallocate(cell_cen_l)
    deallocate(ip_crd_lp)
    deallocate(Xl)
    deallocate(Xg)
    deallocate(Xn)
    deallocate(X1)
    deallocate(X2)
    deallocate(Xv)
    deallocate(Cell_Vertex_Volume)
    deallocate(Projection)
    deallocate(Xlt)
    if (.not. stress_reduced_integration) then
       deallocate(jac)
    end if
    ! Inform the user of successful geometry computation.
    call TLS_info ('done.')
  END SUBROUTINE CELL_CV_FACE

  !---------------------------------------------------------------------------------
  SUBROUTINE BOUNDARY_CV_FACE
    !
    !  Calculate all boundary control volume face areas, and normals and for a cell
    !  This routine also calls routines to put solid mechanics boundary condition
    ! data into the Node_Displacement_BC data structure.
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use legacy_mesh_api,      only: Cell, ncells, nnodes, gather_vertex_coord
    use cutoffs_module,       only: alittle
    use pgslib_module,        only: PGSLib_GLOBAL_ANY
    use node_operator_module, only: CV_Boundary, nbnode, nbface
    use solid_mechanics_mesh, only: ndim, nvc, nvf
    use bc_data_types
    use mech_bc_data_module
    use var_vector_module

    ! Local Variables
    integer :: status
    real(r8), pointer, dimension(:,:,:) :: Xv
    real(r8), Dimension(:,:,:), pointer :: Xn
    real(r8), Dimension(:,:), pointer :: X1, X2
    integer :: atlas_size, i, j, v1, v2, nodecount
    type (BC_Operator), POINTER :: Operator
    type (BC_Atlas),    POINTER :: Atlas
    integer, pointer, dimension(:) :: Cell_List
    integer, pointer, dimension(:) :: Face_List
    ! Boundary stuff 
    real(r8), pointer, Dimension(:,:) :: F_Cent
    integer :: ICell, Face, F_Node
    integer :: ibcop, inode

    !
    ! Inform the user of control volume face calculation.
    call TLS_info ('')
    call TLS_info (' Finding control volume faces on external boundaries ...')
    call TLS_info ('')
    ! Create temporary arrays
    allocate(Xv(ndim,nvc,ncells), stat=status)
    if (status /= 0) call TLS_panic ('BOUNDARY_CV_FACE: allocation error: Xv')
    !
    ! Gather vertex coordinates
    call gather_vertex_coord (Xv)

    ! Get node numbers for the node BC structure
    nodecount = 0
    NBC_index = 0
    do i = 1,nnodes
       if (ANY(Node_Mask(:,i))) then
          nodecount = nodecount + 1
          Node_Displacement_BC%Node(nodecount) = i
          NBC_index(i) = nodecount
       end if
    end do

    ! Loop over all mech bc operators
    BOUNDARY_LOOP: do ibcop = 1, BC_MAX_OPERATORS
       ! Transfer node displacement BCs to control volume data 
       if(PGSLib_GLOBAL_ANY(nbnode(ibcop) /= 0)) then
          ! Displacement BCs
          nodecount = 0
          do inode = 1, nnodes
             if (Node_Mask(ibcop,inode)) then
                nodecount = nodecount + 1
                CV_Boundary(ibcop)%Node(nodecount,1) = inode 
             end if
          end do
       end if
       ! Get (point to) data in "new" bc structures
       Operator => BC_Spec_Get_Operator(Displacement_BC,ibcop)
       Atlas => BC_Op_Get_Atlas(Operator)
       atlas_size = DATA_SIZE(Atlas)
       ! If atlas has face data in it continue
       ATLAS_HAS_DATA: if (atlas_size > 0) then
          Cell_List => BC_Get_Cell(Atlas)
          Face_List => BC_Get_Face(Atlas)
          !
          CV_Boundary(ibcop)%Cell = Cell_List
          CV_Boundary(ibcop)%Face = Face_List
          !
          ! Allocate temp arrays for this boundary
          allocate(Xn(ndim,nvf,nbface(ibcop)), stat=status)
          if (status /= 0) call TLS_panic ('BOUNDARY_CV_FACE: allocation error: Xn')
          allocate(X1(ndim,nbface(ibcop)), stat=status)
          if (status /= 0) call TLS_panic ('BOUNDARY_CV_FACE: allocation error: X1')
          allocate(X2(ndim,nbface(ibcop)), stat=status)
          if (status /= 0) call TLS_panic ('BOUNDARY_CV_FACE: allocation error: X2')
          !
          ! Loop over nodes on a cell face
          INODE_LOOP: do inode = 1, nvf
             ! Loop over faces on a boundary
             NBDFACE_LOOP: do j=1,nbface(ibcop)
                ! Get physical coordinates for each CV face on a cell face, and get list
                ! of nodes for each cell face
                ICell = CV_Boundary(ibcop)%Cell(j)
                Face = CV_Boundary(ibcop)%Face(j)
                F_Cent => Cell(CV_Boundary(ibcop)%Cell(j))%Face_Centroid

                CALL BOUNDARY_FACE_COORDS(inode, j, Xv, ICell, Face, F_Cent, Xn, F_Node)
                
                CV_Boundary(ibcop)%Face_Node(inode,j) = F_Node
  
             end do NBDFACE_LOOP
             ! This is lifted from CELL_CV_FACE, which is adapted from CELL_GEOMETRY
             do i = 1,ndim
                X1(i,:) = Xn(i,1,:) - Xn(i,2,:)
                X2(i,:) = Xn(i,3,:) - Xn(i,4,:)
             end do
             ! These are face area vectors for now.
             do i = 1,ndim
                select case (i)
                case (1)
                   v1 = 2; v2 = 3
                   if (ndim == 2) CV_Boundary(ibcop)%Face_Normal(i,inode,:) = X1(v1,:)
                case (2)
                   v1 = 3; v2 = 1
                   if (ndim == 2) CV_Boundary(ibcop)%Face_Normal(i,inode,:) = X2(v2,:)
                case (3)
                   v1 = 1; v2 = 2
                end select
                if (ndim == 3) CV_Boundary(ibcop)%Face_Normal(i,inode,:) = 0.5_r8*(X1(v1,:)*X2(v2,:) - X2(v1,:)*X1(v2,:))
             end do
             ! Set components to zero if they're small; accumulate areas.
             CV_Boundary(ibcop)%Face_Area(inode,:) = 0.0_r8
             do i = 1,ndim
                CV_Boundary(ibcop)%Face_Normal(i,inode,:) = MERGE(0.0_r8, CV_Boundary(ibcop)%Face_Normal(i,inode,:), &
                     ABS(CV_Boundary(ibcop)%Face_Normal(i,inode,:)) < alittle)
                CV_Boundary(ibcop)%Face_Area(inode,:) = CV_Boundary(ibcop)%Face_Area(inode,:) + &
                     CV_Boundary(ibcop)%Face_Normal(i,inode,:)**2
             end do
             ! Compute areas.
             CV_Boundary(ibcop)%Face_Area(inode,:) = SQRT(CV_Boundary(ibcop)%Face_Area(inode,:))
             ! Convert the face vectors to unit normals.
             do i = 1,ndim
                where (CV_Boundary(ibcop)%Face_Area(inode,:) >= alittle**2)
                   CV_Boundary(ibcop)%Face_Normal(i,inode,:) = CV_Boundary(ibcop)%Face_Normal(i,inode,:) / &
                        CV_Boundary(ibcop)%Face_Area(inode,:)
                elsewhere
                   CV_Boundary(ibcop)%Face_Normal(i,inode,:) = 0.0_r8
                end where
             end do
          end do INODE_LOOP
          ! Deallocate temps because each boundary can have a different number of faces and nodes
          deallocate(Xn)
          deallocate(X1)
          deallocate(X2)
       end if ATLAS_HAS_DATA

       call NODE_BC_DATA(ibcop, nbnode(ibcop), nbface(ibcop))


    end do BOUNDARY_LOOP
    !
    deallocate(Xv)

    ! Calculate new vectors for combined normal constraint nodel BCs
    CALL DISPLACEMENT_CONSTRAINT_VECTORS()

    ! Inform the user of successful geometry computation.
    call TLS_info ('done.')

  END SUBROUTINE BOUNDARY_CV_FACE

  !---------------------------------------------------------------------------------
  SUBROUTINE CV_FACE_COORDS(edge, ndim, Xlc, Xv, cell_cen_l, nn, Xl, Xn, Coef, nsign,  ip_crd_lp,&
                           v11, v12, v13, v14, v21, v22, v23, v24, v31, v32, v33, v34)
    !
    !  Find the coordinates of each nodal control volume and its integration points,
    !  along with other data needed for the CV_Internal data structure
    ! 
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use legacy_mesh_api, only: Cell
    
    !Arguments
    ! Inputs
    integer, intent(IN) :: edge, ndim
    real(r8), Dimension(:,:) , intent(IN) :: Xlc
    real(r8), pointer, dimension(:,:,:) :: Xv
    real(r8), Dimension(:,:), pointer :: cell_cen_l
    ! Outputs
    integer, intent(OUT) :: nn
    real(r8), Dimension(:,:,:),pointer :: Xl
    real(r8), Dimension(:,:,:),pointer :: Xn
    real(r8), Dimension(:), intent(OUT) :: Coef
    real(r8), intent(OUT) :: nsign
    real(r8), Dimension(:,:,:), pointer    :: ip_crd_lp
    integer, intent(OUT) :: v11, v12, v13, v14, v21, v22, v23, v24, v31, v32, v33, v34

    ! Local variables
    integer :: i, n
    
    ! Preinitialize the face centroid
    do n = 1,ndim
       ip_crd_lp(n,edge,:) = 0.5_r8
       Coef(n)          = 1.0_r8
    end do

    select case (edge)
       ! Set up logical sub-cells in four corners of the original cell for coordinate
       ! transforms to map control volume face centroids to logical cell coordinates. 
       ! Each sub-cell has three control volume faces associated with it that correspond
       ! to edges in the original cell.
    case (1,4,8)
       ! Get logical cell vertices of the corner volume associated with these 3 edges.
       do i=1,ndim
          Xl(1,i,:)= Xlc(1,i)
          Xl(2,i,:)= (Xlc(1,i)+Xlc(2,i))/2.0_r8
          Xl(3,i,:)= Cell%Face_centroid_L(i,5)
          Xl(4,i,:)= (Xlc(1,i)+Xlc(4,i))/2.0_r8
          Xl(5,i,:)= (Xlc(1,i)+Xlc(5,i))/2.0_r8
          Xl(6,i,:)= Cell%Face_centroid_L(i,2)
          Xl(7,i,:)= cell_cen_l(i,:)
          Xl(8,i,:)= Cell%Face_centroid_L(i,3)
       end do
       select case (edge)
          ! Get coordinates of control volume face vertices, one per edge.
       case(1)
          ! Get physical coordinates for face and put in order for area/normal calculation
          do i=1,ndim
             Xn(i,1,:)= Cell%Centroid(i)
             Xn(i,2,:)= 0.5_r8*(Xv(i,1,:)+Xv(i,2,:))
             Xn(i,3,:)= Cell%Face_Centroid(i,2)
             Xn(i,4,:)= Cell%Face_Centroid(i,5)
          end do
          ! Get indices of coordinates for logical face centroid calculation
          v11 = 3; v12 = 1; v13 = 2; v14 = 4
          v21 = 1; v22 = 1; v23 = 1; v24 = 1
          v31 = 1; v32 = 4; v33 = 3; v34 = 2
          nn = 2; ip_crd_lp(nn, edge,:) = 1.0_r8; Coef(nn)=0.0_r8; nsign = 1.0_r8
       case (4)
          do i=1,ndim
             Xn(i,1,:)= Cell%Centroid(i)
             Xn(i,2,:)= 0.5_r8*(Xv(i,4,:)+Xv(i,1,:))
             Xn(i,3,:)= Cell%Face_Centroid(i,3)
             Xn(i,4,:)= Cell%Face_Centroid(i,5)
          end do
          ! Get indices of coordinates for logical face centroid calculation
          v11 = 1; v12 = 1; v13 = 1; v14 = 1
          v21 = 1; v22 = 3; v23 = 4; v24 = 2
          v31 = 3; v32 = 2; v33 = 1; v34 = 4
          nn = 1; ip_crd_lp(nn, edge,:) = 0.0_r8; Coef(nn)=0.0_r8; nsign = -1.0_r8
       case (8)
          do i=1,ndim
             Xn(i,1,:)= Cell%Face_Centroid(i,2)
             Xn(i,2,:)= Cell%Face_Centroid(i,3)
             Xn(i,3,:)= Cell%Centroid(i)
             Xn(i,4,:)= 0.5_r8*(Xv(i,1,:)+Xv(i,5,:))
          end do
          v11 = 4; v12 = 2; v13 = 1; v14 = 3
          v21 = 1; v22 = 4; v23 = 3; v24 = 2
          v31 = 1; v32 = 1; v33 = 1; v34 = 1
          nn = 3; ip_crd_lp(nn, edge,:) = 1.0_r8; Coef(nn)=0.0_r8; nsign = 1.0_r8
       end select
    case(7,11,12)
       do i=1,ndim
          Xl(1,i,:)= Cell%Face_centroid_L(i,3)
          Xl(2,i,:)= cell_cen_l(i,:)
          Xl(3,i,:)= Cell%Face_centroid_L(i,1)
          Xl(4,i,:)= (Xlc(8,i)+Xlc(4,i))/2.0_r8
          Xl(5,i,:)= (Xlc(8,i)+Xlc(5,i))/2.0_r8
          Xl(6,i,:)= Cell%Face_centroid_L(i,6)
          Xl(7,i,:)= (Xlc(8,i)+Xlc(7,i))/2.0_r8
          Xl(8,i,:)= Xlc(8,i)
       end do
       select case (edge)
       case (7)
          do i=1,ndim
             Xn(i,1,:)= Cell%Centroid(i)
             Xn(i,2,:)= 0.5_r8*(Xv(i,4,:)+Xv(i,8,:))
             Xn(i,3,:)= Cell%Face_Centroid(i,1)
             Xn(i,4,:)= Cell%Face_Centroid(i,3)
          end do
          v11 = 1; v12 = 3; v13 = 4; v14 = 2
          v21 = 3; v22 = 2; v23 = 1; v24 = 4
          v31 = 1; v32 = 1; v33 = 1; v34 = 1
          nn = 3; ip_crd_lp(nn, edge,:) = 0.0_r8; Coef(nn)=0.0_r8; nsign = -1.0_r8
       case (11)
          do i=1,ndim
             Xn(i,1,:)= Cell%Face_Centroid(i,6)
             Xn(i,2,:)= Cell%Face_Centroid(i,1)
             Xn(i,3,:)= 0.5_r8*(Xv(i,7,:)+Xv(i,8,:))
             Xn(i,4,:)= Cell%Centroid(i)
          end do
          v11 = 1; v12 = 3; v13 = 4; v14 = 2
          v21 = 1; v22 = 1; v23 = 1; v24 = 1
          v31 = 3; v32 = 2; v33 = 1; v34 = 4
          nn = 2; ip_crd_lp(nn, edge,:) = 1.0_r8; Coef(nn)=0.0_r8; nsign = -1.0_r8
       case (12)
          do i=1,ndim
             Xn(i,1,:)= Cell%Face_Centroid(i,6)
             Xn(i,2,:)= Cell%Face_Centroid(i,3)
             Xn(i,3,:)= 0.5_r8*(Xv(i,5,:)+Xv(i,8,:))
             Xn(i,4,:)= Cell%Centroid(i)
          end do
          v11 = 1; v12 = 1; v13 = 1; v14 = 1
          v21 = 4; v22 = 2; v23 = 1; v24 = 3
          v31 = 1; v32 = 4; v33 = 3; v34 = 2
          nn = 1; ip_crd_lp(nn, edge,:) = 1.0_r8; Coef(nn)=0.0_r8; nsign = 1.0_r8
       end select
    case(2,3,6)
       do i=1,ndim
          Xl(1,i,:)= Cell%Face_centroid_L(i,5)
          Xl(2,i,:)= (Xlc(3,i)+Xlc(2,i))/2.0_r8
          Xl(3,i,:)= Xlc(3,i)
          Xl(4,i,:)= (Xlc(3,i)+Xlc(4,i))/2.0_r8
          Xl(5,i,:)= cell_cen_l(i,:)
          Xl(6,i,:)= Cell%Face_centroid_L(i,4)
          Xl(7,i,:)= (Xlc(3,i)+Xlc(7,i))/2.0_r8
          Xl(8,i,:)= Cell%Face_centroid_L(i,1)
       end do
       select case (edge)
       case (2)
          do i=1,ndim
             Xn(i,1,:)= Cell%Centroid(i)
             Xn(i,2,:)= 0.5_r8*(Xv(i,2,:)+Xv(i,3,:))
             Xn(i,3,:)= Cell%Face_Centroid(i,4)
             Xn(i,4,:)= Cell%Face_Centroid(i,5)
          end do
          v11 = 1; v12 = 1; v13 = 1; v14 = 1
          v21 = 2; v22 = 4; v23 = 3; v24 = 1
          v31 = 3; v32 = 2; v33 = 1; v34 = 4
          nn = 1; ip_crd_lp(nn, edge,:) = 1.0_r8; Coef(nn)=0.0_r8; nsign = -1.0_r8
       case (3)
          do i=1,ndim
             Xn(i,1,:)= Cell%Centroid(i)
             Xn(i,2,:)= 0.5_r8*(Xv(i,3,:)+Xv(i,4,:))
             Xn(i,3,:)= Cell%Face_Centroid(i,1)
             Xn(i,4,:)= Cell%Face_Centroid(i,5)
          end do
          v11 = 4; v12 = 2; v13 = 1; v14 = 3
          v21 = 1; v22 = 1; v23 = 1; v24 = 1
          v31 = 1; v32 = 4; v33 = 3; v34 = 2
          nn = 2; ip_crd_lp(nn, edge,:) = 0.0_r8; Coef(nn)=0.0_r8; nsign = 1.0_r8
       case (6)
          do i=1,ndim
             Xn(i,1,:)= Cell%Face_Centroid(i,4)
             Xn(i,2,:)= Cell%Face_Centroid(i,1)
             Xn(i,3,:)= 0.5_r8*(Xv(i,3,:)+Xv(i,7,:))
             Xn(i,4,:)= Cell%Centroid(i)
          end do
          v11 = 4; v12 = 2; v13 = 1; v14 = 3
          v21 = 1; v22 = 4; v23 = 3; v24 = 2
          v31 = 1; v32 = 1; v33 = 1; v34 = 1
          nn = 3; ip_crd_lp(nn, edge,:) = 1.0_r8; Coef(nn)=0.0_r8; nsign = 1.0_r8
       end select
    case(5,9,10)
       do i=1,ndim
          Xl(1,i,:)= Cell%Face_centroid_L(i,2)
          Xl(2,i,:)= (Xlc(6,i)+Xlc(2,i))/2.0_r8
          Xl(3,i,:)= Cell%Face_centroid_L(i,4)
          Xl(4,i,:)= cell_cen_l(i,:)
          Xl(5,i,:)= (Xlc(6,i)+Xlc(5,i))/2.0_r8
          Xl(6,i,:)= Xlc(6,i)
          Xl(7,i,:)= (Xlc(6,i)+Xlc(7,i))/2.0_r8
          Xl(8,i,:)= Cell%Face_centroid_L(i,6)
       end do
       select case (edge)
       case (5)
          do i=1,ndim
             Xn(i,1,:)= 0.5_r8*(Xv(i,2,:)+Xv(i,6,:))
             Xn(i,2,:)= Cell%Centroid(i)
             Xn(i,3,:)= Cell%Face_Centroid(i,4)
             Xn(i,4,:)= Cell%Face_Centroid(i,2)
          end do
          v11 = 1; v12 = 3; v13 = 4; v14 = 2
          v21 = 3; v22 = 2; v23 = 1; v24 = 4
          v31 = 1; v32 = 1; v33 = 1; v34 = 1
          nn = 3; ip_crd_lp(nn, edge,:) = 0.0_r8; Coef(nn)=0.0_r8; nsign = -1.0_r8
       case (9)
          do i=1,ndim
             Xn(i,1,:)= Cell%Face_Centroid(i,6)
             Xn(i,2,:)= Cell%Face_Centroid(i,2)
             Xn(i,3,:)= 0.5_r8*(Xv(i,5,:)+Xv(i,6,:))
             Xn(i,4,:)= Cell%Centroid(i)
          end do
          v11 = 2; v12 = 4; v13 = 3; v14 = 1
          v21 = 1; v22 = 1; v23 = 1; v24 = 1
          v31 = 3; v32 = 2; v33 = 1; v34 = 4
          nn = 2; ip_crd_lp(nn, edge,:) = 0.0_r8; Coef(nn)=0.0_r8; nsign = -1.0_r8
       case (10)
          do i=1,ndim
             Xn(i,1,:)= Cell%Face_Centroid(i,6)
             Xn(i,2,:)= Cell%Face_Centroid(i,4)
             Xn(i,3,:)= 0.5_r8*(Xv(i,6,:)+Xv(i,7,:))
             Xn(i,4,:)= Cell%Centroid(i)
          end do
          v11 = 1; v12 = 1; v13 = 1; v14 = 1
          v21 = 3; v22 = 1; v23 = 2; v24 = 4
          v31 = 1; v32 = 4; v33 = 3; v34 = 2
          nn = 1; ip_crd_lp(nn, edge,:) = 0.0_r8; Coef(nn)=0.0_r8; nsign = 1.0_r8
       end select
    end select

  END SUBROUTINE CV_FACE_COORDS

  !---------------------------------------------------------------------------------
  SUBROUTINE BOUNDARY_FACE_COORDS(inode, jface, Xv, ICell, Face, F_Cent, Xn, F_Node)
    !
    !  Find the coordinates of each nodal control volume face that is on a boundary
    !  and its integration points,  along with other data needed for the CV_Boundary
    !  data structure
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use solid_mechanics_mesh,     only: ndim

    !Arguments
    integer, intent(IN) :: inode, jface, ICell, Face
    real(r8), dimension(:,:), pointer   :: F_Cent
    real(r8), dimension(:,:,:), pointer :: Xv
    real(r8), dimension(:,:,:), pointer :: Xn
    integer, intent(OUT) :: F_Node

    ! Local variables
    integer :: i

    select case (Face)
    case (1)
       select case (inode)
       case (1)
          do i=1,ndim
             Xn(i,1,jface)= Xv(i,8,ICell)
                         Xn(i,2,jface)= F_Cent(i,1)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,7,ICell)+Xv(i,8,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,4,ICell)+Xv(i,8,ICell))
                      end do
                      F_Node = 8
                   case (2)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,7,ICell)+Xv(i,8,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,3,ICell)+Xv(i,7,ICell))
                         Xn(i,3,jface)= Xv(i,7,ICell)
                         Xn(i,4,jface)= F_Cent(i,1)
                      end do
                      F_Node = 7
                   case (3)
                      do i=1,ndim
                         Xn(i,1,jface)= F_Cent(i,1)
                         Xn(i,2,jface)= Xv(i,3,ICell)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,3,ICell)+Xv(i,7,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,3,ICell)+Xv(i,4,ICell))
                      end do
                      F_Node = 3
                   case (4)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,8,ICell)+Xv(i,4,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,3,ICell)+Xv(i,4,ICell))
                         Xn(i,3,jface)= F_Cent(i,1)
                         Xn(i,4,jface)= Xv(i,4,ICell)
                      end do
                      F_Node = 4
                   end select
                case (2)
                   select case (inode)
                   case (1)
                      do i=1,ndim
                         Xn(i,1,jface)= Xv(i,6,ICell)
                         Xn(i,2,jface)= F_Cent(i,2)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,6,ICell)+Xv(i,5,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,6,ICell)+Xv(i,2,ICell))
                      end do
                      F_Node = 6
                   case (2)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,6,ICell)+Xv(i,5,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,1,ICell)+Xv(i,5,ICell))
                         Xn(i,3,jface)= Xv(i,5,ICell)
                         Xn(i,4,jface)= F_Cent(i,2)
                      end do
                      F_Node = 5
                   case (3)
                      do i=1,ndim
                         Xn(i,1,jface)= F_Cent(i,2)
                         Xn(i,2,jface)= Xv(i,1,ICell)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,1,ICell)+Xv(i,5,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,1,ICell)+Xv(i,2,ICell))
                      end do
                      F_Node = 1
                   case (4)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,2,ICell)+Xv(i,6,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,2,ICell)+Xv(i,1,ICell))
                         Xn(i,3,jface)= F_Cent(i,2)
                         Xn(i,4,jface)= Xv(i,2,ICell)
                      end do
                      F_Node = 2
                   end select
                case (3)
                   select case (inode)
                   case (1)
                      do i=1,ndim
                         Xn(i,1,jface)= Xv(i,5,ICell)
                         Xn(i,2,jface)= F_Cent(i,3)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,5,ICell)+Xv(i,8,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,5,ICell)+Xv(i,1,ICell))
                      end do
                      F_Node = 5
                   case (2)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,8,ICell)+Xv(i,5,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,8,ICell)+Xv(i,4,ICell))
                         Xn(i,3,jface)= Xv(i,8,ICell)
                         Xn(i,4,jface)= F_Cent(i,3)
                      end do
                      F_Node = 8
                   case (3)
                      do i=1,ndim
                         Xn(i,1,jface)= F_Cent(i,3)
                         Xn(i,2,jface)= Xv(i,4,ICell)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,4,ICell)+Xv(i,8,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,4,ICell)+Xv(i,1,ICell))
                      end do
                      F_Node = 4
                   case (4)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,1,ICell)+Xv(i,5,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,1,ICell)+Xv(i,4,ICell))
                         Xn(i,3,jface)= F_Cent(i,3)
                         Xn(i,4,jface)= Xv(i,1,ICell)
                      end do
                      F_Node = 1
                   end select
                case (4)
                   select case (inode)
                   case (1)
                      do i=1,ndim
                         Xn(i,1,jface)= Xv(i,7,ICell)
                         Xn(i,2,jface)= F_Cent(i,4)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,7,ICell)+Xv(i,6,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,7,ICell)+Xv(i,3,ICell))
                      end do
                      F_Node = 7
                   case (2)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,6,ICell)+Xv(i,7,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,6,ICell)+Xv(i,2,ICell))
                         Xn(i,3,jface)= Xv(i,6,ICell)
                         Xn(i,4,jface)= F_Cent(i,4)
                      end do
                      F_Node = 6
                   case (3)
                      do i=1,ndim
                         Xn(i,1,jface)= F_Cent(i,4)
                         Xn(i,2,jface)= Xv(i,2,ICell)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,2,ICell)+Xv(i,6,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,2,ICell)+Xv(i,3,ICell))
                      end do
                      F_Node = 2
                   case (4)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,3,ICell)+Xv(i,7,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,3,ICell)+Xv(i,2,ICell))
                         Xn(i,3,jface)= F_Cent(i,4)
                         Xn(i,4,jface)= Xv(i,3,ICell)
                      end do
                      F_Node = 3
                   end select
                case (5)
                   select case (inode)
                   case (1)
                      do i=1,ndim
                         Xn(i,1,jface)= Xv(i,4,ICell)
                         Xn(i,2,jface)= F_Cent(i,5)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,4,ICell)+Xv(i,3,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,4,ICell)+Xv(i,1,ICell))
                      end do
                      F_Node = 4
                   case (2)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,3,ICell)+Xv(i,4,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,3,ICell)+Xv(i,2,ICell))
                         Xn(i,3,jface)= Xv(i,3,ICell)
                         Xn(i,4,jface)= F_Cent(i,5)
                      end do
                      F_Node = 3
                   case (3)
                      do i=1,ndim
                         Xn(i,1,jface)= F_Cent(i,5)
                         Xn(i,2,jface)= Xv(i,2,ICell)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,2,ICell)+Xv(i,3,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,2,ICell)+Xv(i,1,ICell))
                      end do
                      F_Node = 2
                   case (4)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,1,ICell)+Xv(i,4,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,1,ICell)+Xv(i,2,ICell))
                         Xn(i,3,jface)= F_Cent(i,5)
                         Xn(i,4,jface)= Xv(i,1,ICell)
                      end do
                      F_Node = 1
                   end select
                case (6)
                   select case (inode)
                   case (1)
                      do i=1,ndim
                         Xn(i,1,jface)= Xv(i,5,ICell)
                         Xn(i,2,jface)= F_Cent(i,6)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,5,ICell)+Xv(i,6,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,5,ICell)+Xv(i,8,ICell))
                      end do
                      F_Node = 5
                   case (2)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,6,ICell)+Xv(i,5,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,6,ICell)+Xv(i,7,ICell))
                         Xn(i,3,jface)= Xv(i,6,ICell)
                         Xn(i,4,jface)= F_Cent(i,6)
                      end do
                      F_Node = 6
                   case (3)
                      do i=1,ndim
                         Xn(i,1,jface)= F_Cent(i,6)
                         Xn(i,2,jface)= Xv(i,7,ICell)
                         Xn(i,3,jface)= 0.5_r8*(Xv(i,7,ICell)+Xv(i,6,ICell))
                         Xn(i,4,jface)= 0.5_r8*(Xv(i,7,ICell)+Xv(i,8,ICell))
                      end do
                      F_Node = 7
                   case (4)
                      do i=1,ndim
                         Xn(i,1,jface)= 0.5_r8*(Xv(i,8,ICell)+Xv(i,5,ICell))
                         Xn(i,2,jface)= 0.5_r8*(Xv(i,8,ICell)+Xv(i,7,ICell))
                         Xn(i,3,jface)= F_Cent(i,6)
                         Xn(i,4,jface)= Xv(i,8,ICell)
                      end do
                      F_Node = 8
                   end select
                end select

  END SUBROUTINE BOUNDARY_FACE_COORDS

  !---------------------------------------------------------------------------------
  SUBROUTINE NODE_BC_DATA(ibcop, ibcnodes, ibcfaces)

    ! Compute the node BC data
    !
    ! This code depends on the BC_??_DISPLACEMENT_OP integers being in specific order, such that 
    ! all single node displacement constraints are processed before the constraints that involve
    ! interfaces (normal-constraint and contact).  We are looping over the BC operators (ibcop) 
    ! in the order that they are specified in bc_enum_types.F90
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use node_operator_module, only: CV_Boundary
    use mech_bc_data_module
    use legacy_mesh_api,      only: ncells, nnodes, Cell, Mesh, GAP_ELEMENT_1
    use legacy_mesh_api,      only: EN_SUM_SCATTER, NN_Gather_BoundaryData
    use legacy_mesh_api,      only: Vertex, Vertex_Ngbr_All
    use pgslib_module,        only: PGSLib_GLOBAL_ANY, PGSLib_Global_MAXVAL
    use solid_mechanics_mesh, only: ndim, nvc, nvf, nfc
    use var_vector_module

    ! Arguments
    integer, intent(IN) :: ibcop, ibcnodes, ibcfaces

    ! Local variables
    integer :: status
    real(r8), allocatable, dimension(:,:) :: Fnorm, Fcount, Farea
    real(r8), allocatable, dimension(:)   :: Ncount, Narea
    real(r8), allocatable, dimension(:,:) :: Nnorm
    integer :: bctype, inode, nnum, i, j, k, n, nmax, nnbr, ncnt
    logical :: found_node, dup1, dup2, fatal
    real(r8) :: vmag, dist, costheta, sintheta, mindist
    real(r8), pointer, dimension(:,:) :: Bcoords
    real(r8), pointer, dimension(:) :: Btemp
    real(r8), pointer, dimension(:) :: BNcount 
    real(r8), pointer, dimension(:,:) :: BNnorm 
    integer, pointer, dimension(:) :: NN_Vec
    real(r8), dimension(ndim) :: Nbrvec
    real(r8) :: dotprod
    real(r8), dimension(nfc) :: htemp

    if (PGSLib_GLOBAL_ANY(ibcnodes /= 0)) then
       allocate(Fnorm(nvc, ncells), stat= status)
       if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: Fnorm')
       allocate(Farea(nvc, ncells), stat= status)
       if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: Farea')
       allocate(Nnorm(nnodes,ndim), stat= status)
       if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: Nnorm')
       allocate(Fcount(nvc, ncells), stat= status)
       if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: Fcount')
       allocate(Ncount(nnodes), stat= status)           
       if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: Ncount')
       allocate(Narea(nnodes), stat= status)           
       if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: Narea')
       !

       Nnorm = 0.0
       select case (ibcop)
          ! First collect normals, bc types and values for displacement constraints, including 
          ! x,y,z displacements
       case(BC_X_DISPLACEMENT_OP, BC_Y_DISPLACEMENT_OP, BC_Z_DISPLACEMENT_OP)
          select case (ibcop)
          case(BC_X_DISPLACEMENT_OP)
             bctype = X_DISPLACEMENT
             Nnorm(:,1) = 1.0
          case(BC_Y_DISPLACEMENT_OP)
             bctype = Y_DISPLACEMENT
             Nnorm(:,2) = 1.0
          case(BC_Z_DISPLACEMENT_OP)
             bctype = Z_DISPLACEMENT
             Nnorm(:,3) = 1.0
          end select
          
          ! For all nodes in this bc operator
          NODE_DISP_LOOP:do inode = 1, ibcnodes
             nnum = CV_Boundary(ibcop)%Node(inode,1)
             ! Find the first empty normal vector in the CV_Boundary structure
             CV_B_LOOP: do j = 1,ndim
                if (ALL(CV_Boundary(ibcop)%Node_Normal(:,inode,j) == 0.0)) then
                   CV_Boundary(ibcop)%Node_Normal(:,inode,j) = Nnorm(nnum,:)
                   CV_Boundary(ibcop)%Node_Value(inode) = Node_Value(ibcop,nnum)
                   ! Find the first available vector in the node bc structure
                   NODE_BC_LOOP: do k = 1,ndim
                      if (ALL(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) == 0.0)) then
                         Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = Nnorm(nnum,:)
                         Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = bctype
                         Node_Displacement_BC%Value(NBC_index(nnum),k) = Node_Value(ibcop,nnum)
                         select case(k)
                         case(1)
                            Node_Displacement_BC%Combination(NBC_index(nnum)) =  ONE_DISPLACEMENT
                         case(2)
                            ! Check to see if we have the same normal vector
                            ! This is easy since we only have X, Y, or Z directions
                            if (Node_Displacement_BC%BC_Type(NBC_index(nnum),k) /= &
                                 Node_Displacement_BC%BC_Type(NBC_index(nnum),1)) then
                               Node_Displacement_BC%Combination(NBC_index(nnum)) =  TWO_DISPLACEMENTS
                            else
                               ! The two vectors are the same - combine them if the values are the same
                               if (Node_Displacement_BC%Value(NBC_index(nnum),k) == &
                                    Node_Displacement_BC%Value(NBC_index(nnum),1)) then
                                  Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = 0.0
                                  Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = 0
                                  Node_Displacement_BC%Value(NBC_index(nnum),k) = 0.0
                               else
                                  call TLS_panic ('NODE_BC_DATA: two different displacement values specified in the same direction')
                               end if
                            end if
                            
                         case(3)
                            ! Check to see if we have the same normal vector
                            ! This is easy since we only have X, Y, or Z directions
                            fatal = .false.
                            ! Check for duplicate vectors
                            dup1 = (Node_Displacement_BC%BC_Type(NBC_index(nnum),k) == &
                                 Node_Displacement_BC%BC_Type(NBC_index(nnum),1))
                            dup2 = (Node_Displacement_BC%BC_Type(NBC_index(nnum),k) == &
                                 Node_Displacement_BC%BC_Type(NBC_index(nnum),2))
                            if (dup1) then
                               ! The first and third vectors are the same - ignore the third vector
                               ! if the values are the same
                               if (Node_Displacement_BC%Value(NBC_index(nnum),k) == &
                                    Node_Displacement_BC%Value(NBC_index(nnum),1)) then
                                  Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = 0.0
                                  Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = 0
                                  Node_Displacement_BC%Value(NBC_index(nnum),k) = 0.0
                               else
                                  fatal = .true.
                               end if
                            else if (dup2) then
                               ! The second and third vectors are the same - ignore the third vector &
                               ! if the values are the same
                               if (Node_Displacement_BC%Value(NBC_index(nnum),k) == &
                                    Node_Displacement_BC%Value(NBC_index(nnum),2)) then
                                  Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = 0.0
                                  Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = 0
                                  Node_Displacement_BC%Value(NBC_index(nnum),k) = 0.0
                               else
                                  fatal = .true.
                               end if
                            else 
                               Node_Displacement_BC%Combination(NBC_index(nnum)) =  THREE_DISPLACEMENTS
                            end if
                            if (fatal) call TLS_panic ('NODE_BC_DATA: two different displacement values specified in the same direction')
                         end select
                         EXIT CV_B_LOOP
                      end if
                      if(k==ndim) then
                         ! We already have more than ndim normal vectors 
                         ! Issue warning that node is already fully constrained
                         
                         ! Need to check for duplicates at boundaries between element blocks
                         ! This should be OK and handled gracefully
                      end if
                   end do NODE_BC_LOOP
                   EXIT NODE_DISP_LOOP
                end if
             end do CV_B_LOOP
             if(j==ndim) then
                ! We already have more than ndim normal vectors 
                ! Issue warning that node is already fully constrained
             end if
          end do NODE_DISP_LOOP

       case(BC_NORM_DISPLACEMENT_OP)
          bctype = NORMAL_DISPLACEMENT
          INTERFACE_LOOP: do n = 1, SIZE(Interface_List)
             if (Interface_List(n) /= 0) then
                Narea = 0.0
                Farea = 0.0
                do i = 1,ndim
                   Fnorm = 0.0
                   Fcount = 0.0
                   Ncount = 0.0
                   ! Transfer the subface normal component to an array for sum scattering 
                   do j = 1, nvf
                      do k = 1,ibcfaces
                         if ((Interface_Id(CV_Boundary(ibcop)%Face(k),CV_Boundary(ibcop)%Cell(k)) == Interface_List(n)) .and. &
                              (Mesh(CV_Boundary(ibcop)%Cell(k))%Cell_Shape <= GAP_ELEMENT_1)) then
                            Fnorm(CV_Boundary(ibcop)%Face_Node(j,k),CV_Boundary(ibcop)%Cell(k)) = &
                                 CV_Boundary(ibcop)%Face_Normal(i,j,k)
                            ! Face area ia a scalar
                            if (i == 1) Farea(CV_Boundary(ibcop)%Face_Node(j,k),CV_Boundary(ibcop)%Cell(k)) = &
                                 CV_Boundary(ibcop)%Face_Area(j,k)
                            Fcount(CV_Boundary(ibcop)%Face_Node(j,k),CV_Boundary(ibcop)%Cell(k)) = 1.0
                         end if
                      end do
                   end do
                   call EN_SUM_SCATTER(Nnorm(:,i), Fnorm)
                   call EN_SUM_SCATTER(Ncount, Fcount)
                   if (i == 1) call EN_SUM_SCATTER(Narea, Farea)
                   where (Ncount /= 0.0) 
                      Nnorm(:,i) = Nnorm(:,i) / Ncount
                   end where
                end do
                ! Normalize the normal vectors
                do inode = 1, nnodes
                   vmag = sqrt(SUM(Nnorm(inode,:)**2))
                   if (vmag > 0.0) &
                        Nnorm(inode,:) = Nnorm(inode,:)/vmag
                end do

                ! Fill in the mech bc data structures
                NODE_DISP_LOOP2: do inode = 1, ibcnodes
                   nnum = CV_Boundary(ibcop)%Node(inode,1)
                   ! If on this interface (and not a gap element)
                   if (Ncount(nnum) > 0) then
                      CV_B_LOOP2: do j = 1,ndim
                         ! Find the first empty normal vector in the CV_Boundary structure
                         if (ALL(CV_Boundary(ibcop)%Node_Normal(:,inode,j) == 0.0)) then
                            CV_Boundary(ibcop)%Node_Normal(:,inode,j) = Nnorm(nnum,:)
                            CV_Boundary(ibcop)%Node_Value(inode) = Node_Value(ibcop,nnum)
                            ! Find the first available vector in the node bc structure
                            NODE_BC_LOOP2: do k = 1,ndim
                               if (ALL(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) == 0.0)) then
                                  Node_Displacement_BC%Interface(k,NBC_index(nnum)) = n
                                  Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = Nnorm(nnum,:)
                                  Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = bctype
                                  Node_Displacement_BC%Value(NBC_index(nnum),k) = Node_Value(ibcop,nnum)
                                  Node_Displacement_BC%Area(NBC_index(nnum),k) = Narea(nnum)
                                  select case(k)
                                  case(1)
                                     Node_Displacement_BC%Combination(NBC_index(nnum)) =  ONE_DISPLACEMENT
                                  case(2)
                                     ! Check to see if we have the same normal vector
                                     costheta = SUM(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) * &
                                          Node_Displacement_BC%Normal(NBC_index(nnum),1,:))
                                     sintheta = sqrt(1.0-costheta*costheta)
                                     ! Check to see if Nvec1 and Nvec2 are equal or nearly equal
                                     if (sintheta > small_angle) then 
                                        Node_Displacement_BC%Combination(NBC_index(nnum)) =  TWO_DISPLACEMENTS
                                     else
                                        ! The two vectors are the same - combine them if the values are the same,
                                        ! averaging the normal vectors
                                        if (Node_Displacement_BC%Value(NBC_index(nnum),k) == &
                                             Node_Displacement_BC%Value(NBC_index(nnum),1)) then
                                           Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = 0.0
                                           Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = 0
                                           Node_Displacement_BC%Value(NBC_index(nnum),k) = 0.0
                                           ! Average the two vectors
                                           Node_Displacement_BC%Normal(NBC_index(nnum),1,:) = (Nnorm(nnum,:) + &
                                                Node_Displacement_BC%Normal(NBC_index(nnum),1,:)) / 2.0
                                        else
                                           call TLS_panic ('NODE_BC_DATA: two different displacement values specified in the same direction')
                                        end if
                                     end if
                                  case(3)
                                     ! Check to see if we have the same normal vector
                                     fatal = .false.
                                     ! Check for duplicate vectors
                                     costheta = SUM(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) * &
                                          Node_Displacement_BC%Normal(NBC_index(nnum),1,:))
                                     sintheta = sqrt(1.0-costheta*costheta)
                                     ! Check to see if Nvec1 and Nvec3 are equal or nearly equal
                                     dup1 = (sintheta <= small_angle) 
                                     costheta = SUM(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) * &
                                          Node_Displacement_BC%Normal(NBC_index(nnum),2,:))
                                     sintheta = sqrt(1.0-costheta*costheta)
                                     ! Check to see if Nvec2 and Nvec3 are equal or nearly equal
                                     dup2 = (sintheta <= small_angle) 
                                     if (dup1) then
                                        ! The first and third vectors are the same - ignore the third vector
                                        ! if the values are the same
                                        if (Node_Displacement_BC%Value(NBC_index(nnum),k) == &
                                             Node_Displacement_BC%Value(NBC_index(nnum),1)) then
                                           Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = 0.0
                                           Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = 0
                                           Node_Displacement_BC%Value(NBC_index(nnum),k) = 0.0
                                        else
                                           fatal = .true.
                                        end if
                                     else if (dup2) then
                                        ! The second and third vectors are the same - ignore the third vector &
                                        ! if the values are the same
                                        if (Node_Displacement_BC%Value(NBC_index(nnum),k) == &
                                             Node_Displacement_BC%Value(NBC_index(nnum),2)) then
                                           Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = 0.0
                                           Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = 0
                                           Node_Displacement_BC%Value(NBC_index(nnum),k) = 0.0
                                        else
                                           fatal = .true.
                                        end if
                                     else 
                                        Node_Displacement_BC%Combination(NBC_index(nnum)) =  THREE_DISPLACEMENTS
                                     end if
                                     if (fatal) call TLS_panic ('NODE_BC_DATA: two different displacement values specified in the same direction')
                                  end select
                                  EXIT CV_B_LOOP2
                               end if
                               if(k==ndim) then
                                  ! Too many constraints. 
                                  call TLS_panic ('NODE_BC_DATA: too many displacement constraints on at least one node')
                                end if
                            end do NODE_BC_LOOP2
!                            EXIT NODE_DISP_LOOP2
                         end if
                      end do CV_B_LOOP2
                   end if
                end do NODE_DISP_LOOP2
             end if
          end do INTERFACE_LOOP

       case(BC_FREE_INTERFACE_OP, BC_NORMAL_CONSTRAINT_OP, BC_CONTACT_OP)

          select case (ibcop)
          case(BC_FREE_INTERFACE_OP)
             bctype = FREE_INTERFACE
          case(BC_NORMAL_CONSTRAINT_OP)
             bctype = NORMAL_CONSTRAINT
          case(BC_CONTACT_OP)
             bctype = CONTACT
          end select

          ! Gather coordinates of boundary neighbor nodes
          do j = 1,ndim
             NULLIFY(Btemp)
             call NN_GATHER_BOUNDARYDATA (SOURCE=Vertex%coord(j),BOUNDARY=Btemp)
             if (j == 1) then
                allocate(Bcoords(ndim,SIZE(Btemp)), stat = status)
                if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: Bcoords')
             end if
             Bcoords(j,:) = Btemp
             DEALLOCATE(Btemp)          
          end do
          !
          INTERFACE_LOOP2: do n = 1, SIZE(Interface_List)
             if (Interface_List(n) /= 0) then
                Narea = 0.0
                Farea = 0.0
                do i = 1,ndim
                   Fnorm = 0.0
                   Fcount = 0.0
                   Ncount = 0.0
                   ! Transfer the subface normal component to an array for sum scattering 
                   do j = 1, nvf
                      do k = 1,ibcfaces
                         if ((Interface_Id(CV_Boundary(ibcop)%Face(k),CV_Boundary(ibcop)%Cell(k)) == Interface_List(n)) .and. &
                              (Mesh(CV_Boundary(ibcop)%Cell(k))%Cell_Shape <= GAP_ELEMENT_1))  then
                            Fnorm(CV_Boundary(ibcop)%Face_Node(j,k),CV_Boundary(ibcop)%Cell(k)) = &
                                 CV_Boundary(ibcop)%Face_Normal(i,j,k)
                            ! Face area ia a scalar
                            if (i == 1) Farea(CV_Boundary(ibcop)%Face_Node(j,k),CV_Boundary(ibcop)%Cell(k)) = &
                                 CV_Boundary(ibcop)%Face_Area(j,k)
                            Fcount(CV_Boundary(ibcop)%Face_Node(j,k),CV_Boundary(ibcop)%Cell(k)) = 1.0
                         end if
                      end do
                   end do
                   call EN_SUM_SCATTER(Nnorm(:,i), Fnorm)
                   call EN_SUM_SCATTER(Ncount, Fcount)
                   if (i == 1) call EN_SUM_SCATTER(Narea, Farea)
                   where (Ncount /= 0.0) 
                      Nnorm(:,i) = Nnorm(:,i) / Ncount
                   end where
                end do
                ! Normalize the normal vectors
                do inode = 1, nnodes
                   vmag = sqrt(SUM(Nnorm(inode,:)**2))
                   if (vmag > 0.0) &
                        Nnorm(inode,:) = Nnorm(inode,:)/vmag
                end do

                ! Gather Ncount for boundary neighbor nodes so we can tell which off processor nodes
                ! are on this interface
                NULLIFY(Btemp)
                call NN_GATHER_BOUNDARYDATA (SOURCE=Ncount,BOUNDARY=Btemp)
                allocate(BNcount(SIZE(Btemp)), stat = status)
                if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: BNcount')
                BNcount(:) = Btemp
                DEALLOCATE(Btemp)          

                do j = 1,ndim
                   NULLIFY(Btemp)
                   call NN_GATHER_BOUNDARYDATA (SOURCE=Nnorm(:,j),BOUNDARY=Btemp)
                   if (j == 1) then
                      allocate(BNnorm(SIZE(Btemp),ndim), stat = status)
                      if (status /= 0) call TLS_panic ('NODE_BC_DATA: allocation error: BNnorm')
                   end if
                   BNnorm(:,j) = Btemp
                   DEALLOCATE(Btemp)          
                end do

                ! Get gap nodes for all nodes on this interface 

                do i = 1, nfc
                   htemp(i) = MAXVAL(Cell%Halfwidth(i),(Mesh%Cell_Shape >= GAP_ELEMENT_1))
                end do
                mindist = PGSLib_Global_MAXVAL(htemp)
                GAP_NODE: do inode = 1,ibcnodes
                   ! Is the node on this interface?
                   if (Ncount(CV_Boundary(ibcop)%Node(inode,1)) > 0) then
                      found_node = .false.
                      nnum = CV_Boundary(ibcop)%Node(inode,1)
                      nmax = SIZES(Vertex_Ngbr_All(nnum))
                      NN_Vec => FLATTEN(Vertex_Ngbr_All(nnum))
                      ! Now look for a node that is opposite or nearly opposite this one across the interface 
                      NGBR_LOOP: do nnbr =  1,nmax
                         if (NN_Vec(nnbr) > 0) then 
                            ncnt = Ncount(NN_Vec(nnbr))
                         else
                            ncnt = BNcount(abs(NN_Vec(nnbr)))
                         end if
                         if(ncnt > 0) then
                            dist = 0.0
                            do j = 1,ndim
                               if (NN_Vec(nnbr) > 0) then
                                  dist = dist + (Vertex(nnum)%Coord(j) - Vertex(NN_Vec(nnbr))%Coord(j))**2
                               else
                                  dist = dist + (Vertex(nnum)%Coord(j) - Bcoords(j,(-NN_Vec(nnbr))))**2
                               end if
                            end do
                            dist = sqrt(dist)
                            if (NN_Vec(nnbr) > 0) then
                               Nbrvec(:) = Nnorm(NN_Vec(nnbr),:)
                            else
                               Nbrvec(:) = BNnorm(-NN_Vec(nnbr),:)
                            end if
                            dotprod = SUM(Nbrvec(:) * Nnorm(nnum,:))
                            if ((dotprod < -0.9).and.(dist < mindist)) then 
                               if (.not. found_node) then
                                  GNode(nnum) = NN_Vec(nnbr)
                                  found_node = .true.
                               else
                                  call TLS_panic ('NODE_BC_DATA: gap node matches more than one node')
                               end if
                            endif
                         end if
                      end do NGBR_LOOP
                      if (.not. found_node) then
                         GNode(nnum) = nnum
                         found_node = .true.
                         call TLS_info ('')
                         call TLS_info (' Node found at edge of internal gap interface.')
                         call TLS_info ('')
                      end if
                   end if
                end do GAP_NODE
                !
                ! Fill in the mech bc data structures
                NODE_DISP_LOOP3: do inode = 1, ibcnodes
                   nnum = CV_Boundary(ibcop)%Node(inode,1)
                   ! If on this interface
                   if (Ncount(nnum) > 0) then
                      CV_B_LOOP3: do j = 1,ndim
                         ! Find the first empty normal vector in the CV_Boundary structure
                         if (ALL(CV_Boundary(ibcop)%Node_Normal(:,inode,j) == 0.0)) then
                            CV_Boundary(ibcop)%Node(inode,j+1) = GNode(nnum)
                            CV_Boundary(ibcop)%Node_Normal(:,inode,j) = Nnorm(nnum,:)
                            CV_Boundary(ibcop)%Node_Value(inode) = Node_Value(ibcop,nnum)
                            ! Find the first available vector in the node bc structure
                            NODE_BC_LOOP3: do k = 1,ndim
                               if (ALL(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) == 0.0)) then
                                  Node_Displacement_BC%Interface(k,NBC_index(nnum)) = n
                                  Node_Displacement_BC%Gap_Node(NBC_index(nnum),k) = GNode(nnum)
                                  Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = Nnorm(nnum,:)
                                  Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = bctype
                                  Node_Displacement_BC%Value(NBC_index(nnum),k) = Node_Value(ibcop,nnum)
                                  Node_Displacement_BC%Area(NBC_index(nnum),k) = Narea(nnum)
                                  select case(k)
                                  case(1)
                                     Node_Displacement_BC%Combination(NBC_index(nnum)) = ONE_NORM_CONST 
                                  case(2)
                                     ! Check to see if we have the same normal vector
                                     costheta = SUM(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) * &
                                          Node_Displacement_BC%Normal(NBC_index(nnum),1,:))
                                     sintheta = sqrt(1.0-costheta*costheta)
                                     ! Check to see if Nvec1 and Nvec2 are equal or nearly equal
                                     if (sintheta > small_angle) then 
                                        ! Vectors are sufficiently different
                                        select case(Node_Displacement_BC%Combination(NBC_index(nnum)))
                                        case(ONE_DISPLACEMENT)
                                           Node_Displacement_BC%Combination(NBC_index(nnum)) =  ONE_D_ONE_NC
                                        case(ONE_NORM_CONST)
                                           Node_Displacement_BC%Combination(NBC_index(nnum)) =  TWO_NORM_CONST
                                        end select
                                     else
                                        ! Vectors are essentially the same
                                        select case(Node_Displacement_BC%Combination(NBC_index(nnum)))
                                        case(ONE_DISPLACEMENT)
                                           ! Issue warning?  This seems weird, but should work.
                                           Node_Displacement_BC%Combination(NBC_index(nnum)) =  ONE_D_ONE_NC
                                        case(ONE_NORM_CONST)
                                           ! If there are two gap nodes it is OK to have duplicate vectors
                                           if(Node_Displacement_BC%Gap_Node(NBC_index(nnum),k) == &
                                                Node_Displacement_BC%Gap_Node(NBC_index(nnum),1)) then
                                              ! Average the two vectors, treating the two constraints as one
                                              Node_Displacement_BC%Normal(NBC_index(nnum),1,:) = (Nnorm(nnum,:) + &
                                                   Node_Displacement_BC%Normal(NBC_index(nnum),1,:)) / 2.0
                                              Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = 0.0
                                              Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = 0
                                              Node_Displacement_BC%Value(NBC_index(nnum),k) = 0.0
                                           else
                                              Node_Displacement_BC%Combination(NBC_index(nnum)) =  TWO_NORM_CONST
                                           end if
                                        end select
                                     end if
                                  case(3)
                                     ! Check to see if we have the same normal vector
                                     fatal = .false.
                                     ! Check for duplicate vectors
                                     costheta = SUM(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) * &
                                          Node_Displacement_BC%Normal(NBC_index(nnum),1,:))
                                     sintheta = sqrt(1.0-costheta*costheta)
                                     ! Check to see if Nvec1 and Nvec3 are equal or nearly equal
                                     dup1 = (sintheta <= small_angle) 
                                     costheta = SUM(Node_Displacement_BC%Normal(NBC_index(nnum),k,:) * &
                                          Node_Displacement_BC%Normal(NBC_index(nnum),2,:))
                                     sintheta = sqrt(1.0-costheta*costheta)
                                     ! Check to see if Nvec2 and Nvec3 are equal or nearly equal
                                     dup2 = (sintheta <= small_angle) 
                                     select case(Node_Displacement_BC%Combination(NBC_index(nnum)))
                                     case(TWO_DISPLACEMENTS)
                                        ! If there are duplicate vectors then that is a strange combination, but OK(?)
                                        Node_Displacement_BC%Combination(NBC_index(nnum)) = TWO_D_ONE_NC
                                     case(ONE_D_ONE_NC)
                                        ! If dup1 then OK(?)
                                        if ((dup2) .and. (Node_Displacement_BC%Gap_Node(NBC_index(nnum),k) == &
                                             Node_Displacement_BC%Gap_Node(NBC_index(nnum),2))) then
                                           ! If there are two gap nodes it is OK to have duplicate vectors
                                           !
                                           ! The second and third vectors are the same - ignore the third vector
                                           ! if the values are the same
                                           !
                                           ! Average the two vectors, treating the two constraints as one
                                           Node_Displacement_BC%Normal(NBC_index(nnum),2,:) = (Nnorm(nnum,:) + &
                                                Node_Displacement_BC%Normal(NBC_index(nnum),2,:)) / 2.0
                                           Node_Displacement_BC%Normal(NBC_index(nnum),k,:) = 0.0
                                           Node_Displacement_BC%BC_Type(NBC_index(nnum),k) = 0
                                           Node_Displacement_BC%Value(NBC_index(nnum),k) = 0.0
                                        else
                                           Node_Displacement_BC%Combination(NBC_index(nnum)) = ONE_D_TWO_NC
                                        end if
                                     case(TWO_NORM_CONST)
                                        ! We only allow 3 normal constraints if there is only one gap node 
                                        ! and three linearly independent normals (nested 3D corner)
                                        if (dup1.or.dup2) fatal = .true.
                                        if ((Node_Displacement_BC%Gap_Node(NBC_index(nnum),k) /= &
                                             Node_Displacement_BC%Gap_Node(NBC_index(nnum),2)) .or. & 
                                             (Node_Displacement_BC%Gap_Node(NBC_index(nnum),k) /= &
                                             Node_Displacement_BC%Gap_Node(NBC_index(nnum),1))) fatal = .true. 
                                        if (.not. fatal) then
                                           Node_Displacement_BC%Combination(NBC_index(nnum)) = THREE_NORM_CONST
                                        else
                                           if (fatal) call TLS_panic ('NODE_BC_DATA: unsupported combination of three contact or normal constraints')
                                        end if
                                     end select
                                  end select
                                  EXIT CV_B_LOOP3
                               end if
                               if(k==ndim) then
                                  ! Too many constraints.  It might be possible to examine the normals and fix
                                  ! things here, but...
                                  select case(Node_Displacement_BC%Combination(NBC_index(nnum)))
                                  case(THREE_DISPLACEMENTS)
                                     ! No problem, already fully constrained
                                  case(THREE_NORM_CONST)
                                     ! Issue warning about more than ndim normal constraints and continue(?)
                                  case(TWO_D_ONE_NC)
                                     ! Punt (hard to assume much about what is meant here)
                                  case(ONE_D_TWO_NC)
                                     ! Issue warning and continue(?)
                                  end select
                               end if
                            end do NODE_BC_LOOP3
                            EXIT NODE_DISP_LOOP3
                         end if
                      end do CV_B_LOOP3
                   end if
                end do NODE_DISP_LOOP3
                deallocate (BNcount)
             end if
          end do INTERFACE_LOOP2
       end select
       deallocate(Fnorm)
       deallocate(Nnorm)
       deallocate(Narea)
       deallocate(Fcount)
       deallocate(Ncount)
       deallocate(Farea)
    end if

  END SUBROUTINE NODE_BC_DATA

  !---------------------------------------------------------------------------------
  SUBROUTINE DISPLACEMENT_CONSTRAINT_VECTORS()
    !
    ! Purpose:
    ! Calculate vectors needed for constraint equations for nodes with two or 
    ! three normal constraints. Check for some consistency problems and some 
    ! degenerate cases.  Also get normal vectors for interface between two
    ! different gap nodes.
    !
    ! Author: David Korzekwa
    !---------------------------------------------------------------------------------
    use mech_bc_data_module
    use bc_data_types
    use legacy_mesh_api,      only: nnodes
    use lu_solve_module,      only: LU_SOLVE, factsolve
    use solid_mechanics_mesh, only: ndim

    ! Local variables
    integer :: inode, idim, solve_flag, nnum, status
    real(r8), dimension(ndim) :: Nvec1, Nvec2, Gvec1, Gvec2, Gvec3, &
                                 Tvec1, Tvec2, Tvec3, Vvec, Wvec, D, B, Avec
    integer, dimension(ndim) :: indx
    real(r8), dimension(ndim,ndim) :: M
    real(r8), allocatable, dimension(:,:,:) :: All_Norm
    real(r8) :: costheta, sintheta, vmag
    logical :: pivot = .true.

    ! Gather normals from off-processor nodes int a boundary array so we can
    ! calculate lambda3 if we need to.
    !
    ! Put normals into full nodal arrays
    Allocate(All_Norm(nnodes, ndim, ndim), stat= status)
    if (status /= 0) call TLS_panic ('DISPLACEMENT_CONSTRAINT_VECTORS: allocation error: ALL_Norm')
    All_Norm = -2.0
    do inode = 1,SIZE(Node_Displacement_BC%Node)
       nnum = Node_Displacement_BC%Node(inode)
       do idim = 1,ndim
          ALL_Norm(nnum,idim,:) = Node_Displacement_BC%Normal(inode,idim,:)
       end do
    end do

    do inode = 1,SIZE(Node_Displacement_BC%Node)
       select case(Node_Displacement_BC%Combination(inode))
       case(ONE_DISPLACEMENT)
          ! Nothing  needed

       case(TWO_DISPLACEMENTS)
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          Nvec2 = Node_Displacement_BC%Normal(inode,2,:)
          ! Compute vector Avec from the two normal vectors 
          D(1) = Node_Displacement_BC%Value(inode,1)
          D(2) = Node_Displacement_BC%Value(inode,2)
          costheta = SUM(Nvec1 * Nvec2)
          sintheta = sqrt(1.0-costheta*costheta)
          B(1) = (D(1) - costheta * D(2))/sintheta
          B(2) = (D(2) - costheta * D(1))/sintheta
          Avec(:) = B(1)*Nvec1(:) + B(2) * Nvec2(:)
          ! Put result in first vector 
          Node_Displacement_BC%Vector(inode,1,:) = Avec
          ! Compute and store the tangent vector in the second vector
          Tvec1 = cross_product(Nvec1, Nvec2)
          vmag = sqrt(SUM(Tvec1(:)**2))
          Tvec1 = Tvec1/vmag
          Node_Displacement_BC%Vector(inode,2,:) = Tvec1

       case(THREE_DISPLACEMENTS)
          !Solve the system [n1 n2 n3}^T A = D
          M(1,:) = Node_Displacement_BC%Normal(inode,1,:)
          M(2,:) = Node_Displacement_BC%Normal(inode,2,:)
          M(3,:) = Node_Displacement_BC%Normal(inode,3,:)
          B(:) = Node_Displacement_BC%Value(inode,:)
          solve_flag=factsolve
          call LU_SOLVE (M,B,indx,solve_flag,pivot)
          ! Put result in first vector 
          Node_Displacement_BC%Vector(inode,1,:) = B

       case(ONE_NORM_CONST)
          ! Nothing  needed

       case(TWO_NORM_CONST)
          ! We have normal vectors for each interface and need a tangent vector
          Gvec1 = Node_Displacement_BC%Normal(inode,1,:)
          Gvec2 = Node_Displacement_BC%Normal(inode,2,:)
          ! Check to see if we have the same normal vector
          costheta = SUM(Gvec1(:) * Gvec2(:))
          sintheta = sqrt(1.0-costheta*costheta)
          ! If Gvec1 and Gvec2 are sufficiently different
          if (sintheta > small_angle) then 
             ! Compute and store the tangent vector in the first vector
             Tvec1 = cross_product(Gvec1, Gvec2)
             vmag = sqrt(SUM(Tvec1(:)**2))
             Tvec1 = Tvec1/vmag
             Node_Displacement_BC%Vector(inode,1,:) = Tvec1
          else
             ! Average the normals and make them the same.  Currently we assume that 
             ! if the normals are the same then we do not need lambda3. 
             Node_Displacement_BC%Normal(inode,1,:) = (Gvec1(:) + Gvec2(:))/2.0
             Node_Displacement_BC%Normal(inode,2,:) = Node_Displacement_BC%Normal(inode,1,:)
             ! The tangent vector is undefined
             Node_Displacement_BC%Vector(inode,1,:) = 0.0
          end if

       case(THREE_NORM_CONST)
          ! We have normal vectors for each interface and need three 
          ! tangent vectors
          Gvec1 = Node_Displacement_BC%Normal(inode,1,:)
          Gvec2 = Node_Displacement_BC%Normal(inode,2,:)
          Gvec3 = Node_Displacement_BC%Normal(inode,3,:)
          ! Compute and store the tangent vectors
          Tvec1 = cross_product(Gvec1, Gvec2)
          vmag = sqrt(SUM(Tvec1(:)**2))
          Tvec1 = Tvec1/vmag
          Node_Displacement_BC%Vector(inode,1,:) = Tvec1
          !
          Tvec2 = cross_product(Gvec2, Gvec3)
          vmag = sqrt(SUM(Tvec2(:)**2))
          Tvec2 = Tvec2/vmag
          Node_Displacement_BC%Vector(inode,2,:) = Tvec2
          !
          Tvec3 = cross_product(Gvec3, Gvec1)
          vmag = sqrt(SUM(Tvec3(:)**2))
          Tvec3 = Tvec3/vmag
          Node_Displacement_BC%Vector(inode,3,:) = Tvec3


       case(ONE_D_ONE_NC)
          ! Surface normal
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          ! Gap normal
          Gvec1 = Node_Displacement_BC%Normal(inode,2,:)
          ! Tangent vector
          Tvec1 = cross_product(Nvec1, Gvec1)
          vmag = sqrt(SUM(Tvec1(:)**2))
          if (vmag > 0.0) Tvec1 = Tvec1/vmag
          Node_Displacement_BC%Vector(inode,1,:) = Tvec1
          ! Vector to complete an orthogonal set
          Vvec = cross_product(Tvec1, Nvec1)
          Node_Displacement_BC%Vector(inode,2,:) = Vvec

          ! If this appears to be a symmetry plane, use v for the gap normal
          if(Node_Displacement_BC%Value(inode,1) == 0.0) &
               Node_Displacement_BC%Normal(inode,2,:) = Vvec

          ! Cosine of angle between the V vector and the gap normal
          costheta = SUM(Vvec * Gvec1)
          Node_Displacement_BC%Scalar(inode,1) = costheta

       case(TWO_D_ONE_NC)
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          Nvec2 = Node_Displacement_BC%Normal(inode,2,:)
          Gvec1 = Node_Displacement_BC%Normal(inode,3,:)
          ! Compute vector Avec from the two normal vectors 
          D(1) = Node_Displacement_BC%Value(inode,1)
          D(2) = Node_Displacement_BC%Value(inode,2)
          costheta = SUM(Nvec1 * Nvec2)
          B(1) = (D(1) - costheta * D(2))/(1.0 - costheta * costheta)
          B(2) = (D(2) - costheta * D(1))/(1.0 - costheta * costheta)
          do idim = 1, ndim
             Avec(idim) = B(1)*Nvec1(idim) + B(2) * Nvec2(idim)
          end do
          ! Put result in first vector 
          Node_Displacement_BC%Vector(inode,1,:) = Avec
          ! Compute and store the tangent vector in second vector
          Tvec1 = cross_product(Nvec1, Nvec2)
          vmag = sqrt(SUM(Tvec1(:)**2))
          Tvec1 = Tvec1/vmag
          Node_Displacement_BC%Vector(inode,2,:) = Tvec1

          ! Cosine of angle between the tangent vector and the gap normal
          costheta = SUM(Tvec1 * Gvec1)
          Node_Displacement_BC%Scalar(inode,1) = costheta

          ! If these appear to be symmetry planes, use t for the gap normal
          if((Node_Displacement_BC%Value(inode,1) == 0.0) .and. (Node_Displacement_BC%Value(inode,2) == 0.0)) then
             Node_Displacement_BC%Normal(inode,3,:) = Tvec1
             if (costheta < 0.0) Node_Displacement_BC%Normal(inode,3,:) = -Tvec1
          end if
   
       Case(ONE_D_TWO_NC)
          ! We have one surface normal vector, ...
          Nvec1 = Node_Displacement_BC%Normal(inode,1,:)
          ! and two gap normal vectors
          Gvec1 = Node_Displacement_BC%Normal(inode,2,:)
          Gvec2 = Node_Displacement_BC%Normal(inode,3,:)
          ! Check to see if we have the same normal vector
          costheta = SUM(Gvec1(:) * Gvec2(:))
          sintheta = sqrt(1.0-costheta*costheta)
          ! If Gvec1 and Gvec2 are essentially the same
          if (sintheta <= small_angle) then 
             ! Average the normals and make them the same
             ! This will cause Vvec and Wvec to be equal
             Node_Displacement_BC%Normal(inode,2,:) = (Gvec1(:) + Gvec2(:))/2.0
             Node_Displacement_BC%Normal(inode,3,:) = Node_Displacement_BC%Normal(inode,2,:)
             Gvec1 = Node_Displacement_BC%Normal(inode,2,:)
             Gvec2 = Gvec1
          end if

          ! We need vectors v and w, and two cosine scalars 
          ! Compute the tangent vectors
          Tvec1 = cross_product(Nvec1, Gvec1)
          vmag = sqrt(SUM(Tvec1(:)**2))
          Tvec1 = Tvec1/vmag
          !
          Tvec2 = cross_product(Nvec1, Gvec2)
          vmag = sqrt(SUM(Tvec2(:)**2))
          Tvec2 = Tvec2/vmag
          ! Vectors v and w are orthogonal to t1,n1 and t2,n1 respectively
          Vvec =  cross_product(Tvec1, Nvec1)
          Node_Displacement_BC%Vector(inode,1,:) = Vvec
          Wvec =  cross_product(Tvec2, Nvec1)
          Node_Displacement_BC%Vector(inode,2,:) = Wvec

          ! If this appears to be a symmetry plane, use v and w for the gap normals
          if(Node_Displacement_BC%Value(inode,1) == 0.0) then 
             Node_Displacement_BC%Normal(inode,2,:) = Vvec
             Node_Displacement_BC%Normal(inode,3,:) = Wvec
          end if

          ! cosine(theta_1) and cosine(theta_2)
          Node_Displacement_BC%Scalar(inode,1) = SUM(Vvec * Gvec1)
          Node_Displacement_BC%Scalar(inode,2) = SUM(Wvec * Gvec2)

       end select
    end do
    Deallocate(All_Norm)

  end SUBROUTINE DISPLACEMENT_CONSTRAINT_VECTORS

  !---------------------------------------------------------------------------------
  !! Copied from MeshSupport
  pure function cross_product (a, b) result (axb)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: axb(3)
    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  !---------------------------------------------------------------------------------
  SUBROUTINE CELL_LOGICAL_CENTROID (Xv,cell_cen_l)
    !
    ! Calculate the logical cell centroid given vertex coordinates and the cell volume
    !
    ! Author: David Korzekwa
    !---------------------------------------------------------------------------------
    use legacy_mesh_api, only: ncells, Cell
    use solid_mechanics_mesh, only: ndim, nrot

    ! Arguments
    real(r8), Dimension(:,:,:), intent(IN) :: Xv
    real(r8), Dimension(:,:),   intent(OUT) :: cell_cen_l
    ! Local Variables
    integer :: i, i1, i2, v1, v2, v3, v4, v5, v6, v7, v8,status
    real(r8), dimension(:,:), pointer ::L, M, N, D1, D2, D3, Dv
    real(r8), dimension(:,:), pointer :: LxD3, MxD2, NxD1, D1xDv, D2xDv, D3xDv
    ! Explicitly allocate temporaries
    allocate (L(ndim,ncells), stat=status)
    if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: L')
    allocate (M(ndim,ncells),stat=status )
    if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: M')
    allocate (N(ndim,ncells),stat=status)
    if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: N')
    allocate (LxD3(nrot,ncells),stat=status)
    if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: LxD3')
    allocate (MxD2(nrot,ncells),stat=status)
    if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: MxD2')
    allocate (NxD1(nrot,ncells),stat=status)
    if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: NxD1')
    ! Only need these for 3-D.
    if (ndim == 3) then
       allocate (D1(ndim,ncells),stat=status)
       if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: D1')
       allocate (D2(ndim,ncells),stat=status)
       if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: D2')
       allocate (D3(ndim,ncells),stat=status)
       if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: D3')
       allocate (Dv(ndim,ncells),stat=status)
       if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: Dv')
       allocate (D1xDv(nrot,ncells),stat=status)
       if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: D1xDv')
       allocate (D2xDv(nrot,ncells),stat=status)
       if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: D2xDv')
       allocate (D3xDv(nrot,ncells),stat=status)
       if (status /= 0) call TLS_panic ('CELL_LOGICAL_CENTROID: allocation error: D3xDv')
    end if

    ! Set vertex numbers
    select case(ndim)
    case (2)
       v1 = 1; v2 = 2; v3 = 3; v4 = 4
    case (3)
       v1 = 1; v2 = 2; v3 = 3; v4 = 4
       v5 = 5; v6 = 6; v7 = 7; v8 = 8
    end select
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
       L = L/4
       M = M/4
       N = N/4
       cell_cen_l = 0.0_r8
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
          cell_cen_l(v1,:) = cell_cen_l(v1,:) + L(i,:)*(MxD2(i,:) - NxD1(i,:)) &
               + (N(i,:)*D2xDv(i,:) - M(i,:)*D1xDv(i,:))/12
          cell_cen_l(v2,:) = cell_cen_l(v2,:) + M(i,:)*(NxD1(i,:) - LxD3(i,:)) &
               + (L(i,:)*D1xDv(i,:) - N(i,:)*D3xDv(i,:))/12
          cell_cen_l(v3,:) = cell_cen_l(v3,:) + N(i,:)*(LxD3(i,:) - MxD2(i,:)) &
               + (M(i,:)*D3xDv(i,:) - L(i,:)*D2xDv(i,:))/12
       end do
       do i = 1,ndim
          cell_cen_l(i,:) = 0.5_r8 + cell_cen_l(i,:)/(24*Cell%Volume)
       end do
    end select
    ! Explicitly allocate temporaries
    deallocate (L)
    deallocate (M )
    deallocate (N)
    deallocate (LxD3)
    deallocate (MxD2)
    deallocate (NxD1)
    ! Only need these for 3-D.
    if (ndim == 3) then
       deallocate (D1)
       deallocate (D2)
       deallocate (D3)
       deallocate (Dv)
       deallocate (D1xDv)
       deallocate (D2xDv)
       deallocate (D3xDv)
    end if
  END SUBROUTINE CELL_LOGICAL_CENTROID

END MODULE NODE_OP_SETUP_MODULE
