MODULE NODE_OPERATOR_MODULE
  !======================================================================
  ! Purpose:
  !
  !   Define data structures and procedures associated with the node based
  !   operator for thermo-mechanics
  !=======================================================================
  ! Nxtot is for the boundary kludge
  use parameter_module, only: ndim, nvc, nvf
  use kind_module,      only: int_kind, real_kind, log_kind
  use bc_data_types,    only: BC_MAX_OPERATORS 
  implicit none
  ! Private Module
  private
  ! Public types and variables
  public :: CV_Face_Internal, &
            CV_Face_Boundary, &
            CV_Internal,      &
            CV_Boundary,      &
            nipc, nipbf,      &
            stress_reduced_integration, &
            cv_init, &
            Nodal_Volume, &
            mech_precond_init
  ! Boundary variables
  public :: nbface, nbnode, nmechbc
  ! Number of valid BC operators - currently x,y,z tractions and x,y,z displacements,
  ! normal displacement, normal traction, free (unconstrained) interface, normally
  ! constrained interface and contact.  Not all of these are fully implemented.
  integer(KIND = int_kind), parameter                 :: nmechbc = 11
  ! Boundary face and node counters, one for each BC operator
  integer(KIND = int_kind), save,  dimension(BC_MAX_OPERATORS) :: nbface, nbnode
  ! Integration points per cell
  integer(KIND = int_kind), parameter :: nipc = nvc*ndim/2
  ! Integration points per boundary face
  integer(KIND = int_kind), parameter :: nipbf = nvf
  !
  ! Flag to call setup routines only once
  Logical(KIND =  log_kind),save            :: cv_init=.false.
  ! Reduced integration flag
  Logical(KIND =  log_kind),save            :: stress_reduced_integration
  ! Flag to call precondition setup routines only when needed
  Logical(KIND =  log_kind),save            :: mech_precond_init=.false.
  ! Control volume faces internal to the mesh cells
  type CV_Face_Internal
     ! Normal vector for control volume face
     real    (KIND=real_kind), pointer, Dimension(:,:,:)   :: Face_Normal
     ! Area of control volume face
     real    (KIND=real_kind), pointer, Dimension(:,:)     :: Face_Area
     ! IP coordinates for control volume face
     real    (KIND=real_kind), pointer, Dimension(:,:,:)   :: Face_Coord
     ! Inverse jacobian for gradient
     real    (KIND=real_kind), pointer, Dimension(:,:,:,:) :: Face_Ijac
  end type CV_Face_Internal
  ! Control volume faces on boundaries
  type CV_Face_Boundary
     ! Cell that contains these integration points
     integer (KIND=int_kind),  pointer, Dimension(:)       :: Cell
     ! Nodes associated with these integration points
     integer (KIND=int_kind),  pointer, Dimension(:,:)     :: Node
     ! BC values for these nodes, only one component
     real    (KIND=real_kind), pointer, Dimension(:)       :: Node_Value
     ! Surface normal components for these nodes
     real    (KIND=real_kind), pointer, Dimension(:,:,:)   :: Node_Normal
     ! Number of interfaces (and normal vectors) for these nodes
     integer (KIND=int_kind), pointer, Dimension(:)        :: NN_Count
     ! Cell faces that contains these CV faces
     integer (KIND=int_kind),  pointer, Dimension(:)       :: Face
     ! Node numbers associated with these CV faces
     integer (KIND=int_kind),  pointer, Dimension(:,:)     :: Face_Node
     ! Normal vectors for control volume faces
     real    (KIND=real_kind), pointer, Dimension(:,:,:)   :: Face_Normal
     ! Areas of control volume faces
     real    (KIND=real_kind), pointer, Dimension(:,:)     :: Face_Area
     ! IP coordinates for control volume faces
     real    (KIND=real_kind), pointer, Dimension(:,:,:)   :: Face_Coord
     ! Inverse jacobian for gradient
     real    (KIND=real_kind), pointer, Dimension(:,:,:,:) :: Face_Ijac
  end type CV_Face_Boundary
  ! Create instances of control volume derived types
  type(CV_Face_Internal),save                             :: CV_Internal
  ! Control volume faces on boundaries
  type(CV_Face_Boundary),save, allocatable, dimension(:)  :: CV_Boundary
  ! Nodal volume array for preconditioner
  real (KIND=real_kind), pointer, Dimension(:)                     :: Nodal_Volume
END MODULE NODE_OPERATOR_MODULE

