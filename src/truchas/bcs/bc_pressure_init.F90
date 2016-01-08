!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module BC_Pressure_Init

  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the special case of BC stuff needed to test/develop
  !   use of BCs for pressure.
  !   Lots of stuff is hardwired.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@lanl.gov)
  !         Sriram Swaminarayan
  !-----------------------------------------------------------------------------
  use bc_data_types
  use bc_initialize
  use parameter_module
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  Private
  PUBLIC :: Initialize_Pressure_BC

CONTAINS
  subroutine Initialize_Pressure_BC(Pressure_BC)
    !-----------------------------------------------------------------------------
    ! Purpose:
    !   Initialize whatever is needed for the pressure BC.
    !   Eventually this routine will be vanish and a general
    !   BC initialization will be done.
    !
    !-----------------------------------------------------------------------------

    ! Arguments
    type (BC_Specifier), intent(INOUT), target :: Pressure_BC

    ! Announce what is going on
    call TLS_info (' Initializing Pressure Boundary Conditions')

    ! Initialize the BC specifier
    call INITIALIZE(Pressure_BC, 'PRESSURE', BC_PRESSURE_ID)

    !!!!!!!!!! Dirichlet operator !!!!!!!!!!
    call TLS_info ('   Constructing Dirichlet Operator ... ', TLS_VERB_NOISY)
    call Set_Pressure_Operator(Pressure_BC, BC_DIRICHLET_Op, 'DIRICHLET')

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The REFLECTIVE operator uses a region similar to HNeumann, but
    ! points are reflected from cell centers across HNeumann faces.
    call TLS_info ('   Constructing Reflective Operator ... ', TLS_VERB_NOISY)
    call Set_Pressure_Operator(Pressure_BC, BC_REFLECTIVE_Op, 'REFLECTIVE')


    !!! EXTERIOR !!!
    call TLS_info ('   Constructing EXTERIOR Operator ... ', TLS_VERB_NOISY)
    call Set_Pressure_Operator(Pressure_BC, BC_EXTERIOR_Op, 'EXTERIOR')

    !!! NEUMANN !!!
    call TLS_info ('   Constructing Neumann Operator ... ', TLS_VERB_NOISY)
    call Set_Pressure_Operator(Pressure_BC, BC_NEUMANN_Op, 'NEUMANN')

    !!! HNEUMANN !!!
    call TLS_info ('   Constructing HNeumann Operator ... ', TLS_VERB_NOISY)
    call Set_Pressure_Operator(Pressure_BC, BC_HNEUMANN_Op, 'HNEUMANN')

    call TLS_info (' Pressure BCs initialized.')

  end subroutine Initialize_Pressure_BC

  subroutine Set_Pressure_Operator(BC_Spec, BC_Operator_ID, BC_Name_String)
    use bc_regions
    use bc_data_types
    use bc_initialize
    use parallel_util_module
    use pgslib_module, only: PGSLib_Global_SUM
    
    type(BC_Specifier), intent(INOUT), target :: BC_Spec
    integer, intent(IN) :: BC_Operator_ID
    character(*), intent(IN) :: BC_Name_String
    
    ! Local variables
    type(BC_Operator), POINTER :: Operator
    type(BC_Region),   POINTER :: Region
    type(BC_Atlas),    POINTER :: Atlas
    integer :: NumBoundaryPoints
    character(128) :: message
    
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    ! Get the  operator
    Operator => BC_Spec_Get_Operator(BC_Spec, BC_Operator_ID)
    
    ! For now, all operators start life ACTIVE
    call BC_Op_Set_State(Operator, BC_OP_ACTIVE)

    call TLS_info ('   Building the Region ... ', .false., TLS_VERB_NOISY)
    
    ! Get the Region, so we can initialize that
    Region => BC_OP_Get_Region(Operator)
    
    ! Initialize the region, and put the data into it.
    call Set_Pressure_BC_Region(REGION = Region, BC_OPERATOR_ID = BC_Operator_ID)
         
    call TLS_info ('done.', TLS_VERB_NOISY)
    
    ! Build the Atlas from the Region data
    
    call TLS_info ('   Constructing the Atlas ... ', .false., TLS_VERB_NOISY)

    ! Get the atlas and initialize it
    Atlas => BC_OP_Get_Atlas(Operator)
    
    call BC_Atlas_From_Region(ATLAS = Atlas, REGION = Region)
    
    call TLS_info ('done.', TLS_VERB_NOISY)
    
    call TLS_info ('   Canonicalizing Atlas ... ', .false., TLS_VERB_NOISY)
    call CANONICALIZE(Atlas)
    call TLS_info ('done.', TLS_VERB_NOISY)
    
    NumBoundaryPoints = PGSLib_Global_SUM(DATA_SIZE(Atlas) )
    write(message,'(3x,2a,i0,a)') trim(adjustl(bc_name_string)), ': ', &
                                  NumBoundaryPoints, ' boundary points'
    call TLS_info (message)
    
  end subroutine Set_Pressure_Operator

  subroutine Set_Pressure_BC_Region(Region, BC_Operator_ID)
    !-----------------------------------------------------------------------------
    ! Purpose:
    !   Initialize whatever is needed for the pressure Region
    !   The pressure BC determinants are not particularly clear,
    !   this code is subject to lots of changes.
    !   Right now this scans the whole mesh and for each BC type (BC_Operator_ID)
    !   and sets up one region at a time.  Eventually will probably want to set up
    !   all the regions at once, so we only walk through mesh once.
    !
    !-----------------------------------------------------------------------------
    use bc_module,            only: BC_P_REFLECTIVE, DIRICHLET, BC, Prs, BC_Prs
    use mesh_module,          only: Cell, Mesh

    ! Arguments
    type(BC_Region), intent(INOUT) :: Region
    integer, intent(IN ) :: BC_Operator_ID

    ! Local variables
    logical, dimension(nfc,ncells) :: Face_Is_External, Face_Is_Dirichlet
    logical, dimension(nfc,ncells) :: BC_Mask
    real(r8), dimension(1,nfc,ncells) :: BC_Values
    logical,  dimension(nfc,ncells) :: BC_UseFunction
    real(r8), dimension(ndim, nfc, ncells) :: BC_Positions

    integer :: f, c

    ! Initialize the Region to be ready to accumulate the data
    call INITIALIZE(Region)
    ! The easiest way to build up a region is to start with nothing
    ! and keep appending data.  For that we need to start with a 0 sized region.
    ! For now, we assert that all types of pressure BCs have only a single DOF.
    call ALLOC(Region, SIZE=0, DIMENSIONALITY = ndim, DOF = 1)

    ! For reflected BCs
    ! positions are reflected across face centroid.  Position of Bdy Pt B is
    ! B = Xface + (Xface - Xcell) = 2*Xface - Xcell
    ! For all other BCs, position is face center itself.

    ! To build the region, use the already established pressure BC
    ! data structures.  These are:
    ! Dirichlet:
    !    DIRICHLET(BC%Flag, Prs%Face_bit(f)) == .TRUE.
    ! Homogeneous Neumann
    !    All external boundaries (Mesh(c)%Ngbr_cell(f) == 0) which
    !    are NOT Dirichlet

    ! We need to know which are dirichlet and which are external faces
    do f = 1, nfc
       Face_Is_Dirichlet(f,:) = DIRICHLET(BC%Flag, Prs%Face_Bit(f))
       Face_Is_External (f,:) = (Mesh%Ngbr_Cell(f) == 0)
    end do

    BC_Values    = 0.0
    BC_Positions = 0.0
    ! For pressure, right now never use a function.
    BC_UseFunction = .FALSE.
    do c = 1, ncells
       do f = 1, nfc
          ! Assume this isn't a boundary face
          BC_Mask(f,c) = .FALSE.

          ! Is this a boundary face that we want to process?
          select case (BC_Operator_ID)
          case (BC_DIRICHLET_Op)
             if (Face_Is_Dirichlet(f,c)) then
                BC_Mask(f,c)   = .TRUE.
                BC_Values(1,f,c) = BC_Prs(f,c)
                BC_Positions(:,f,c) = Cell(c)%Face_Centroid(:,f)
             end if
             CYCLE
             
          case (BC_HNEUMANN_Op)
             if (Face_Is_External(f,c) .AND. (.NOT. Face_Is_Dirichlet(f,c))) then
                BC_Mask(f,c)   = .TRUE.
                BC_Values(1,f,c) = 0.0
                BC_Positions(:,f,c) = Cell(c)%Face_Centroid(:,f)
             end if
             CYCLE
             
          case (BC_P_REFLECTIVE)
             ! For now we are putting a reflective BC on every external face.
             ! We do this to avoid singularities in the LSLR calculations
             if (Face_Is_External(f,c) ) then
                BC_Mask(f,c)   = .TRUE.
                ! Don't actually need values for Reflective case
                BC_Values(1,f,c) = 0.0
                ! BC_Positions are reflected across face
                BC_Positions(:, f, c) = 2.0*Cell(c)%Face_Centroid(:,f) - Cell(c)%Centroid(:)
             end if
             CYCLE

          case (BC_EXTERIOR_Op)
             ! Identify every exterior face.  This is used in the pressure calculations,
             ! among other places.  The position is just the face centroid.
             ! External faces have no neighbor.
             if (Face_Is_External(f,c) ) then
                BC_Mask(f,c)   = .TRUE.
                ! Don't need values, since that is context dependent.
                BC_Values(1,f,c) = 0.0
                ! BC_Positions are face centroids
                BC_Positions(:, f, c) = Cell(c)%Face_Centroid(:,f)
             end if
             CYCLE

          case DEFAULT
            ! Didn't find the BC we were looking for
            CYCLE
          end select

       end do
    end do
    

    ! Append this to the region (which is empty, since this is only
    ! thing getting added in this example.)
    call APPEND(Region, BC_Mask, BC_Values, BC_UseFunction, POSITIONS=BC_Positions)

    call CANONICALIZE(Region)

  end subroutine Set_Pressure_BC_Region

END Module BC_Pressure_Init


    
    
    
