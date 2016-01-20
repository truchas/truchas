!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_Displacement_Init
  !
  !-----------------------------------------------------------------------------
  ! Purpose:
  ! Setup solid mechanics BCs in conjunction with BC_INIT in init_module.
  ! First step towards using the "new" BCs independent of the old bc data
  ! structures.
  !
  ! Provides:
  ! Initialize_Displacement_BC
  ! Append_to_Displacement_BC
  ! Make_Displacement_BC_Atlases
  !
  ! Documentation - yeah, right.
  !
  ! Author: Dave Korzekwa
  !         Sriram Swaminarayan
  !-----------------------------------------------------------------------------
  use bc_data_types
  use bc_module
  use mech_bc_data_module
  use bc_initialize
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  Private
  PUBLIC :: Initialize_Displacement_BC, Append_to_Displacement_BC, &
       Make_Displacement_BC_Atlases, Node_Set_BC_Init, Displacement_BC, &
       Interface_Surface_Id

CONTAINS
  SUBROUTINE Initialize_Displacement_BC
    !-----------------------------------------------------------------------------
    ! Purpose:
    !   Initialize whatever is needed for the Displacement BC.
    !
    ! Author: Dave Korzekwa
    !-----------------------------------------------------------------------------
    use mech_bc_data_module
    use parameter_module, only: mbc_surfaces, mbc_nodes
    use legacy_mesh_api, only: ndim
    !
    !
    ! Arguments
    ! Local variables
    integer :: status
    type (BC_Operator), POINTER :: Operator
    type (BC_Region),   POINTER :: Region
    !
    ! Announce what is going on
    call TLS_info (' Initializing Displacement Boundary Conditions ')
    !
    ! Initialize the BC specifier
    call INITIALIZE(Displacement_BC, 'Displacement', BC_DISPLACEMENT_ID)
    !
    ! Initialize regions for all operators.  If they end up not being used, they
    ! will just have zero size.
    !
    ! Get the  operator for x traction
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_X_TRACTION_OP)
    !
    ! For now, all operators start life ACTIVE
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    !
    ! Get the Region, so we can initialize that
    Region => BC_OP_Get_Region(Operator)
    !
    call INITIALIZE(Region)
    ! The easiest way to build up a region is to start with nothing
    ! and keep appending data.  For that we need to start with a 0 sized region.
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    ! Do the same for the other operators
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Y_TRACTION_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Z_TRACTION_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_X_DISPLACEMENT_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Y_DISPLACEMENT_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Z_DISPLACEMENT_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORM_DISPLACEMENT_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORM_TRACTION_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_FREE_INTERFACE_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORMAL_CONSTRAINT_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)
    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_CONTACT_OP)
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)
    Region => BC_OP_Get_Region(Operator)
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)
    !
    NULLIFY(Operator, Region)

    ! Initialize nodal bc data
    allocate(Node_Disp_BC_Temp%Node(mbc_nodes*mbc_surfaces), stat=status)
    if (status /= 0) call TLS_panic ('Node_Set_BC_Init: allocation error: Node_Disp_BC_Temp%Node')
    allocate(Node_Disp_BC_Temp%Gap_Node(mbc_nodes*mbc_surfaces,ndim), stat=status)
    if (status /= 0) call TLS_panic ('Node_Set_BC_Init: allocation error: Node_Disp_BC_Temp%Gap_Node')
    allocate(Node_Disp_BC_Temp%BC_Type(mbc_nodes*mbc_surfaces,ndim), stat=status)
    if (status /= 0) call TLS_panic ('Node_Set_BC_Init: allocation error: Node_Disp_BC_Temp%BC_Type')
    allocate(Node_Disp_BC_Temp%Normal(mbc_nodes*mbc_surfaces,ndim,ndim), stat=status)
    if (status /= 0) call TLS_panic ('Node_Set_BC_Init: allocation error: Node_Disp_BC_Temp%Normal')
    allocate(Node_Disp_BC_Temp%Combination(mbc_nodes*mbc_surfaces), stat=status)
    if (status /= 0) call TLS_panic ('Node_Set_BC_Init: allocation error: Node_Disp_BC_Temp%Combination')
    allocate(Node_Disp_BC_Temp%Value(mbc_nodes*mbc_surfaces,ndim), stat=status)
    if (status /= 0) call TLS_panic ('Node_Set_BC_Init: allocation error: Node_Disp_BC_Temp%Value')
    Node_Disp_BC_Temp%Node        = 0
    Node_Disp_BC_Temp%Gap_Node    = 0
    Node_Disp_BC_Temp%BC_Type     = 0
    Node_Disp_BC_Temp%Normal      = 0.0
    Node_Disp_BC_Temp%Combination = 0
    Node_Disp_BC_Temp%Value       = 0.0
    !
  end SUBROUTINE Initialize_Displacement_BC
  !
  !
  !
  SUBROUTINE Append_to_Displacement_BC(BC_Mask, n, f)
    !-------------------------------------------------------------------
    ! Purpose:
    !  Add displacement BC from a single BC namelist to the appropriate
    !  "new" bc operator.
    !
    ! Author: Dave Korzekwa
    !---------------------------------------------------------------------
    !
    use legacy_mesh_api, only: ndim, nfc, ncells, Cell
    use bc_data_module, only: BC_Type, BC_Value
    !
    ! Arguments
    logical, dimension(nfc,ncells), intent(IN) :: BC_Mask
    integer, intent(IN) :: n, f ! The index for this bc namelist, and the face
    !
    ! Local Variables
    type (BC_Operator), POINTER :: Operator
    type (BC_Region),   POINTER :: Region
    real(r8), dimension(1,nfc,ncells) :: BC_Values
    logical,  dimension(nfc,ncells) :: BC_UseFunction
    real(r8), dimension(ndim, nfc, ncells) :: BC_Positions
    !
    ! Initialize local arrays
    !
    ! Careful - these variable are very similar in name to the BC namelist variables.
    !
    BC_Values = 0.0
    BC_UseFunction = .false.
    BC_Positions = -1.0
    !
    ! Get the operator for this namelist
    !
    select case (trim(BC_TYPE(n)))
       case('x-traction')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_X_TRACTION_OP)
       case('y-traction')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Y_TRACTION_OP)
       case('z-traction')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Z_TRACTION_OP)
       case('x-displacement')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_X_DISPLACEMENT_OP)
       case('y-displacement')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Y_DISPLACEMENT_OP)
       case('z-displacement')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Z_DISPLACEMENT_OP)
       case('normal-displacement')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORM_DISPLACEMENT_OP)
       case('normal-traction')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORM_TRACTION_OP)
       case('free-interface')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_FREE_INTERFACE_OP)
       case('normal-constraint')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORMAL_CONSTRAINT_OP)
       case('contact')
          Operator => BC_Spec_Get_Operator(Displacement_BC, BC_CONTACT_OP)
    end select
    !
    if (f == 1) then
       call TLS_info ('')
       call TLS_info ('   Adding '//trim(BC_TYPE(n))//' operator ... ', .false., TLS_VERB_NOISY)
    end if
    !
    ! Get the region
    !
    Region => BC_OP_Get_Region(Operator)
    where (BC_Mask(f,:))
       ! Only room for 1 value in the new BC data structure
       BC_Values(1,f,:) = BC_Value(1,n)
       BC_UseFunction(f,:) = .false.
       ! Position doesn't get used for displacements or tractions for now.  Set to face
       ! centroids.
       BC_Positions(1,f,:) = Cell(:)%Face_Centroid(1,f)
       BC_Positions(2,f,:) = Cell(:)%Face_Centroid(2,f)
       BC_Positions(3,f,:) = Cell(:)%Face_Centroid(3,f)
    end where
    !
    ! Add to the data for this region/operator and distribute
    !
    call APPEND(Region, BC_Mask, BC_Values, BC_UseFunction, POSITIONS=BC_Positions)
    !
    call CANONICALIZE(Region)
    !
    nullify(Operator)
    nullify(Region)
  end SUBROUTINE Append_to_Displacement_BC
  !
  !
  !
  SUBROUTINE Make_Displacement_BC_Atlases
    !---------------------------------------------------------------
    ! Purpose:
    !  Take the regions constructed from the bc input and make the 
    !  corresponding atlases'
    !
    ! Author: Dave Korzekwa
    !
    !----------------------------------------------------------------

    ! Local Variables
    type(BC_Operator), POINTER :: Operator
    type(BC_Region),   POINTER :: Region
    type(BC_Atlas),    POINTER :: Atlas

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_X_TRACTION_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made X-traction atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Y_TRACTION_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made Y-traction atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Z_TRACTION_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made Z-traction atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_X_DISPLACEMENT_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made X-displacement atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Y_DISPLACEMENT_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made Y-displacement atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_Z_DISPLACEMENT_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made Z-displacement atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORM_DISPLACEMENT_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made normal-displacement atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORM_TRACTION_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made normal-traction atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_FREE_INTERFACE_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made free interface atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_NORMAL_CONSTRAINT_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made normal constraint atlas.', TLS_VERB_NOISY)

    Operator => BC_Spec_Get_Operator(Displacement_BC, BC_CONTACT_OP)
    Region => BC_OP_Get_Region(Operator)
    Atlas => BC_OP_Get_Atlas(Operator)
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    call CANONICALIZE(Atlas)
    call TLS_info (' Made contact atlas.', TLS_VERB_NOISY)

  end SUBROUTINE Make_Displacement_BC_Atlases
  !
  !-----------------------------------------------------------------------------
  !
  SUBROUTINE Node_Set_BC_Init(m,n)
    !
    ! Purpose:
    !  Take the coordinates for a set of nodes from the BC input and find any nodes
    !  on this processor that match.
    !
    ! Author(s): Dave Korzekwa
    !
    !-----------------------------------------------------------------------------
    use legacy_mesh_api,      only: ndim, nnodes, Vertex, UnPermute_Vertex_Vector
    use bc_data_module,       only: BC_Value, BC_Type, Conic_Tolerance, &
                                    Node_Disp_Coords
    use mech_bc_data_module,  only: Node_Disp_BC_Temp, X_DISPLACEMENT, Y_DISPLACEMENT, Z_DISPLACEMENT, &
                                    NORMAL_DISPLACEMENT
    use parameter_module,     only: mbc_nodes
    use pgslib_module,        only: PGSLib_GLOBAL_ANY, PGSLib_COLLATE
    use parallel_info_module
    !
    ! Arguments
    integer, intent(IN) :: m, n
    ! Local variables
    integer :: status, i, j, idim, inode, bc_gnode
    integer, pointer, dimension(:) :: Collated_Nodes
    real(r8) :: dist
    logical :: node_found
    character(128) :: message
    !
    if (p_info%IOP) then
       ALLOCATE(Collated_Nodes(p_info%nPE))
    else
       ALLOCATE(Collated_Nodes(0))
    end if
    NODAL_BC_LOOP: do i = 1,mbc_nodes
       ! See if data has been changed from default
       NODE_COORDS: if ((Node_Disp_Coords(1,i,n) > -1.0d10)) then
          node_found = .false.
          bc_gnode = 0
          ! Loop over cells and vertices
          NODE_LOOP: do inode = 1, nnodes
             dist = 0.0
             do idim = 1,ndim
                dist = dist + (Vertex(inode)%Coord(idim) - Node_Disp_Coords(idim,i,n))**2
             end do
             dist = sqrt(dist)
             ! If the node is close enough to the specified point and is on this processor...
             if (ABS(dist) <= Conic_Tolerance(m, n)) then
                ! Only use the first index for this temporary collection of node BC info
                nbc_nodes = nbc_nodes + 1
                Node_Disp_BC_Temp%Node(nbc_nodes) = inode
                Node_Disp_BC_Temp%Value(nbc_nodes,1) = BC_Value(1,n)
                select case (BC_Type(n))
                case ('x-displacement')
                   Node_Disp_BC_Temp%BC_Type(nbc_nodes,1) = X_DISPLACEMENT
                case ('y-displacement')
                   Node_Disp_BC_Temp%BC_Type(nbc_nodes,1) = Y_DISPLACEMENT
                case ('z-displacement')
                   Node_Disp_BC_Temp%BC_Type(nbc_nodes,1) = Z_DISPLACEMENT
                case ('normal-displacement')
                   Node_Disp_BC_Temp%BC_Type(nbc_nodes,1) = NORMAL_DISPLACEMENT
                end select
                bc_gnode = UnPermute_Vertex_Vector(inode)
                node_found = .true.
!                EXIT NODE_LOOP
             end if
          end do NODE_LOOP
          ! Check to see if a node was found.
          if (PGSLib_GLOBAL_ANY(node_found)) then
             call PGSLib_COLLATE (Collated_Nodes, bc_gnode)
!             node_found = .false.
             do j = 1, SIZE(Collated_Nodes)
                if (Collated_Nodes(j) /= 0) then
!                   if (.not. node_found) then
                      write (message,101) Collated_Nodes(j), TRIM(BC_Type(n))
101                   format(9x,'Node number',i5,' matches bc coordinates for ',a,&
                           ' boundary condition')
                      call TLS_info (message)
!                      node_found = .true.
!                   else
!                      Output_String = blank_line
!                      write (Output_String,103) Node_Disp_Coords(:,i,n)
!103                   format (/,' FATAL: Nodes found on more than one processor for coordinates ',3es12.4/)
!                      call PUNT (Output_String, 'Node_Set_BC_Init')
!                   end if
                end if
             end do
          else
             write(message,'(a,3es12.4)') 'Node_Set_BC_Init: no node found for nodal BC coordinates ', Node_Disp_Coords(:,i,n)
             call TLS_fatal (message)
          end if
       end if NODE_COORDS
    end do NODAL_BC_LOOP
    deallocate(Collated_Nodes)
  end SUBROUTINE Node_Set_BC_Init
!
  !-----------------------------------------------------------------------------
  !
  SUBROUTINE Interface_Surface_Id(Mask,p,f,msrf)
    !
    ! Purpose:
    !  Collect mesh surface id info for multiple surface with interface boundary 
    !  conditions.  This is used to handle corners and edges between interfaces
    !  where constraints need to be added instead of using an average normal.
    !
    ! Author(s): Dave Korzekwa
    !
    !-----------------------------------------------------------------------------
    use mech_bc_data_module,  only: Interface_ID, Interface_List
    use legacy_mesh_api, only: ncells, nfc
    use parallel_info_module
!    use legacy_mesh_api,          only: Mesh, GAP_ELEMENT_1
    !
    ! Arguments
    integer, intent(IN) :: p, f, msrf
    logical, dimension(nfc,ncells), intent(IN) :: Mask
    character(128) :: message

    ! Local Variables
    integer :: n, icell, status

    if (.not. ASSOCIATED(Interface_ID)) then
       allocate(Interface_ID(nfc,ncells), stat = status)
       if (status /= 0) call TLS_panic ('INTERFACE_SURFACE_ID: allocation error: Interface_ID')
       Interface_ID = 0 
    end if

    ! Add to list of interface ids
    do n = 1, SIZE(Interface_List)
       if (Interface_List(n) == msrf) EXIT
       if (Interface_List(n) == 0) then
          Interface_List(n) = msrf
          EXIT
       end if
    end do

    ! Assign each face to an interface id
    do icell = 1, ncells
       if (Mask(f,icell)) then
!       if (Mask(f,icell) .and. (Mesh(icell)%Cell_Shape < GAP_ELEMENT_1)) then
          if (Interface_ID(f,icell) == 0) then
             Interface_ID(f,icell) = msrf
          else
             write(message,'(a,2(i0,a))') 'Interface_Surface_Id: ', icell, ' face ', f, ' is on more than one interface'
             call TLS_panic (message)
          end if
       end if
    end do
  END SUBROUTINE Interface_Surface_Id

end Module BC_Displacement_Init

