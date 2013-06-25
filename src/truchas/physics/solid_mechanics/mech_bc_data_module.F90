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
  use kind_module,      only: int_kind, real_kind, log_kind
  use bc_data_types
  implicit none

  !  private
  ! Public module
  public
  !:: Nodal_BC, Displacement_BC, nbc_nodes, Node_Displacement_BC, &
  !       Interface_ID, Interface_Lis
  !
  ! Variables needed in processing the namelist input
  !
  ! Define a structure for nodal BC's
  type Nodal_BC
     integer (KIND=int_kind), pointer, dimension(:)     :: Node        ! Node number
     integer (KIND=int_kind), pointer, dimension(:,:)   :: Gap_Node    ! Node(s) across gap(s)
     integer (KIND=int_kind), pointer, dimension(:,:)   :: BC_Type     ! Type of nodal constraint
     real (KIND=real_kind), pointer, dimension(:,:,:)   :: Normal      ! Normal vectors for each separate constraint
     integer (KIND=int_kind), pointer, dimension(:)     :: Combination ! Identifier for combination of multiple BCs
     real (KIND=real_kind), pointer, dimension(:,:)     :: Value       ! BC values (displacements or friction parameters)
     integer (KIND=int_kind), pointer, dimension(:,:)   :: Interface   ! BC values (displacements or friction parameters)
     ! Vectors and scalars needed by the displacement constraint algorithms.  These 
     ! are different for each combination of constraints, and will change if the 
     ! mesh changes.
     real (KIND=real_kind), pointer, dimension(:,:,:)   :: Vector      ! Vectors constructed for the constraint algorithm
     real (KIND=real_kind), pointer, dimension(:,:)     :: Scalar      ! Scalars such as cosines between vectors
     real (KIND=real_kind), pointer, dimension(:,:)     :: Lambda      ! Contact function, between 0 (no contact) 
                                                                       ! and 1 (full contact)
     real (KIND=real_kind), pointer, dimension(:,:)     :: Gap_Disp    ! The gap displacement, if applicable. Positive
                                                                       ! corresponds to an open gap.
     real (KIND=real_kind), pointer, dimension(:,:)     :: Normal_Traction   ! The normal traction at this node.  Only applicable to gap
                                                                       ! nodes for now.
     real (KIND=real_kind), pointer, dimension(:,:)     :: Area        ! The surface area associated with this node.
  end type Nodal_BC
  ! Nodal_BC instance for all nodal BCs
  type(Nodal_BC), save    :: Node_Displacement_BC
  ! Nodal BC type for collecting individual BCs prior to combining with nodes specified by face sets
  type(Nodal_BC), save    :: Node_Disp_BC_Temp
  !
  ! Displacement BC data in "new" BC Atlas format
  type (BC_Specifier), save, target            :: Displacement_BC
  !
  integer (KIND=int_kind), save                         :: nbc_nodes = 0  ! Number of nodal bcs
  integer (KIND=int_kind), pointer, dimension(:,:)      :: Interface_ID => Null()    ! Interface surface ids
  integer (KIND=int_kind), save, dimension(50)          :: Interface_List = 0 ! List of Surface ids
  ! Combination Identifiers
  integer (KIND=int_kind), parameter         :: TRACTION_ONLY       = 0  
  integer (KIND=int_kind), parameter         :: ONE_DISPLACEMENT    = 1  
  integer (KIND=int_kind), parameter         :: TWO_DISPLACEMENTS   = 2
  integer (KIND=int_kind), parameter         :: THREE_DISPLACEMENTS = 3
  integer (KIND=int_kind), parameter         :: ONE_NORM_CONST      = 4  
  integer (KIND=int_kind), parameter         :: TWO_NORM_CONST      = 5
  integer (KIND=int_kind), parameter         :: THREE_NORM_CONST    = 6
  integer (KIND=int_kind), parameter         :: ONE_D_ONE_NC        = 7
  integer (KIND=int_kind), parameter         :: TWO_D_ONE_NC        = 8
  integer (KIND=int_kind), parameter         :: ONE_D_TWO_NC        = 9

  ! BC types for solid mechanics
  integer (KIND=int_kind), parameter         :: X_TRACTION          = 1
  integer (KIND=int_kind), parameter         :: Y_TRACTION          = 2
  integer (KIND=int_kind), parameter         :: Z_TRACTION          = 3
  integer (KIND=int_kind), parameter         :: X_DISPLACEMENT      = 4
  integer (KIND=int_kind), parameter         :: Y_DISPLACEMENT      = 5
  integer (KIND=int_kind), parameter         :: Z_DISPLACEMENT      = 6
  integer (KIND=int_kind), parameter         :: NORMAL_DISPLACEMENT = 7
  integer (KIND=int_kind), parameter         :: NORMAL_TRACTION     = 8
  integer (KIND=int_kind), parameter         :: FREE_INTERFACE      = 9
  integer (KIND=int_kind), parameter         :: NORMAL_CONSTRAINT   = 10
  integer (KIND=int_kind), parameter         :: CONTACT             = 11

  integer (KIND=int_kind), pointer, dimension(:)  :: NBC_index

  real (KIND=real_kind), pointer, dimension(:)       :: Face_Gap => Null()
  integer (KIND=int_kind), pointer, dimension(:,:,:) :: Face_Node => Null()

END MODULE MECH_BC_DATA_MODULE
