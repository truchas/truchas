!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MECH_BC_DATA_MODULE
!=======================================================================
  ! Purpose(s):
  !
  !   Defines data variables for solid mechanics boundary conditions. 
  !   Note: This module is declared public and the variables are saved 
  !   upon use of the module.
  !
  !   Currently a quick hack to set up test problems for the T-M operator
  !
  !
  ! Author(s): Dave Korzekwa (dak@lanl.gov)
  !=======================================================================
  use kinds, only: r8
  use bc_data_types
  implicit none
  public
  
  !:: Nodal_BC, Displacement_BC, nbc_nodes, Node_Displacement_BC, &
  !       Interface_ID, Interface_Lis
  !
  ! Variables needed in processing the namelist input
  !
  ! Define a structure for nodal BC's
  type Nodal_BC
     integer, pointer, dimension(:)     :: Node        ! Node number
     integer, pointer, dimension(:,:)   :: Gap_Node    ! Node(s) across gap(s)
     integer, pointer, dimension(:,:)   :: BC_Type     ! Type of nodal constraint
     real(r8), pointer, dimension(:,:,:)   :: Normal      ! Normal vectors for each separate constraint
     integer, pointer, dimension(:)     :: Combination ! Identifier for combination of multiple BCs
     real(r8), pointer, dimension(:,:)     :: Value       ! BC values (displacements or friction parameters)
     integer, pointer, dimension(:,:)   :: Interface   ! BC values (displacements or friction parameters)
     ! Vectors and scalars needed by the displacement constraint algorithms.  These 
     ! are different for each combination of constraints, and will change if the 
     ! mesh changes.
     real(r8), pointer, dimension(:,:,:)   :: Vector      ! Vectors constructed for the constraint algorithm
     real(r8), pointer, dimension(:,:)     :: Scalar      ! Scalars such as cosines between vectors
     real(r8), pointer, dimension(:,:)     :: Lambda      ! Contact function, between 0 (no contact) 
                                                                       ! and 1 (full contact)
     real(r8), pointer, dimension(:,:)     :: Gap_Disp    ! The gap displacement, if applicable. Positive
                                                                       ! corresponds to an open gap.
     real(r8), pointer, dimension(:,:)     :: Normal_Traction   ! The normal traction at this node.  Only applicable to gap
                                                                       ! nodes for now.
     real(r8), pointer, dimension(:,:)     :: Area        ! The surface area associated with this node.
  end type Nodal_BC
  ! Nodal_BC instance for all nodal BCs
  type(Nodal_BC), save    :: Node_Displacement_BC
  ! Nodal BC type for collecting individual BCs prior to combining with nodes specified by face sets
  type(Nodal_BC), save    :: Node_Disp_BC_Temp
  !
  ! Displacement BC data in "new" BC Atlas format
  type (BC_Specifier), save, target            :: Displacement_BC
  !
  integer, save                         :: nbc_nodes = 0  ! Number of nodal bcs
  integer, pointer, dimension(:,:)      :: Interface_ID => Null()    ! Interface surface ids
  integer, save, dimension(50)          :: Interface_List = 0 ! List of Surface ids
  ! Combination Identifiers
  integer, parameter         :: TRACTION_ONLY       = 0  
  integer, parameter         :: ONE_DISPLACEMENT    = 1  
  integer, parameter         :: TWO_DISPLACEMENTS   = 2
  integer, parameter         :: THREE_DISPLACEMENTS = 3
  integer, parameter         :: ONE_NORM_CONST      = 4  
  integer, parameter         :: TWO_NORM_CONST      = 5
  integer, parameter         :: THREE_NORM_CONST    = 6
  integer, parameter         :: ONE_D_ONE_NC        = 7
  integer, parameter         :: TWO_D_ONE_NC        = 8
  integer, parameter         :: ONE_D_TWO_NC        = 9

  ! BC types for solid mechanics
  integer, parameter         :: X_TRACTION          = 1
  integer, parameter         :: Y_TRACTION          = 2
  integer, parameter         :: Z_TRACTION          = 3
  integer, parameter         :: X_DISPLACEMENT      = 4
  integer, parameter         :: Y_DISPLACEMENT      = 5
  integer, parameter         :: Z_DISPLACEMENT      = 6
  integer, parameter         :: NORMAL_DISPLACEMENT = 7
  integer, parameter         :: NORMAL_TRACTION     = 8
  integer, parameter         :: FREE_INTERFACE      = 9
  integer, parameter         :: NORMAL_CONSTRAINT   = 10
  integer, parameter         :: CONTACT             = 11

  integer, pointer, dimension(:)  :: NBC_index

  real(r8), pointer, dimension(:)       :: Face_Gap => Null()
  integer, pointer, dimension(:,:,:) :: Face_Node => Null()

END MODULE MECH_BC_DATA_MODULE
