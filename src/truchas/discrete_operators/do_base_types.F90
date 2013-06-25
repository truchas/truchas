MODULE DO_BASE_TYPES
  !=======================================================================
  ! Purpose(s):
  !   Define base data structures for discrete operations
  !   such as gradient, divergence, etc. DO_Specifier is the
  !   "container" pointer structure instantiated to "hold"
  !   all of the persistent solve vectors, matrices, etc
  !   associated with a call to the gradient operators from
  !   a given model component (such as "projection").
  !
  !   Public Data Structures: DO_Specifier
  !                           SVDU_Type
  !                           dX_Type
  !                           SField_Type
  !                           NIdx_Type
  !
  !   Public Parameters:      DO_SOLVE_DEFAULT
  !                           DO_SOLVE_ORTHO
  !                           DO_SOLVE_LU_LSLR
  !                           DO_SOLVE_SVD_LSLR
  !                           DO_NUM_ST
  !
  ! Contains:
  !
  ! Author(s): Doug Kothe (dbk@lanl.gov)
  !            Jeff Durachta (durachta@verizon.net)
  !            Robert Ferrell (ferrell@cpca.com)
  !
  !=======================================================================
  use bc_data_types
  use kind_module,       only: real_kind, int_kind, log_kind
  use var_vector_module, only: REAL_VAR_VECTOR

  implicit none

  ! Private Module
  private

  ! Public Data Structures
  public :: DO_Specifier, DO_Diag_Specifier, DO_DiagListSpec, &
            SVDU_Type, dX_Type, SField_Type, NIdx_Type, cMat_row_type

  ! Public Parameters
  public :: DO_SOLVE_DEFAULT
  public :: DO_SOLVE_ORTHO
  public :: DO_SOLVE_LU_LSLR
  public :: DO_SOLVE_SVD_LSLR
  public :: DO_NUM_ST


  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Internal data structures

  ! Structure Primitives

  ! Structure to hold set of SVD U Matrices
    type SVDU_Type
      real(KIND=real_kind),pointer,dimension(:,:) :: Mat=>NULL()
    end type SVDU_Type

  ! Structure to hold set of geometric dX vectors
    type dX_Type
      real(KIND=real_kind),pointer,dimension(:,:) :: FData=>NULL()
    end type dX_Type

  ! Structure to hold set of scalar field values
    type SField_Type
      real(KIND=real_kind),pointer,dimension(:)   :: FData=>NULL()
    end type SField_Type

  ! Structure to hold index to face neighbors
    type NIdx_Type
       integer(KIND=int_kind),pointer,dimension(:) :: ptr=>NULL()
    end type NIdx_Type 

  type DO_Specifier
! The contents of this structure have been made "public" in
! order to avoid the complication of "accessor functions".
! The "public" use of the data structure contents is to be
! limitied to within the discrete operator (do_*.F90) set
! of modules. This is done for reasons of transparency and
! performance.
! Developers must NOT directly access any of this data structure's
! parts nor should this data structure be modified in any way outside
! of the discrete operators infrastructure.
    integer(KIND=int_kind)                   :: ID=-9999
    integer(KIND=int_kind)                   :: Method=-9999

    type(dX_Type)    ,pointer,dimension(:,:) :: dX_Struct   =>NULL()
    type(SField_Type),pointer,dimension(:,:) :: TotW_Struct =>NULL()
    type(SField_Type),pointer,dimension(:,:) :: GeoW_Struct =>NULL()
    type(SField_Type),pointer,dimension(:,:) :: PHI_Struct  =>NULL()
    type(NIdx_Type)  ,pointer,dimension(:)   :: NF_Struct   =>NULL()
    type(NIdx_Type)  ,pointer,dimension(:)   :: NFN_Struct  =>NULL()
    type(NIdx_Type)  ,pointer,dimension(:,:) :: NghIdx      =>NULL()

  ! Ortho Solution
    real(KIND=real_kind),pointer,dimension(:,:,:) :: dX_Scaled     =>NULL()
    real(KIND=real_kind),pointer,dimension(:)     :: W_Ortho       =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:)   :: W_Ortho_Nghbr =>NULL()
    
  ! LSLR Data Structures
    logical(KIND=log_kind),pointer,dimension(:,:)     :: UpdateDecomp =>NULL()
    logical(KIND=log_kind),pointer,dimension(:,:,:)   :: SolveFlag    =>NULL()
    
  ! LSLR via LU w/ full pivot
    real(KIND=real_kind)  ,pointer,dimension(:,:,:,:) :: ALU     =>NULL()
    logical(KIND=log_kind),pointer,dimension(:,:,:)   :: PivFlag =>NULL()
    integer(KIND=int_kind),pointer,dimension(:,:,:)   :: row1    =>NULL()
    integer(KIND=int_kind),pointer,dimension(:,:,:)   :: row2    =>NULL()
    
  ! LSLR via SVD
    type(SVDU_Type)     ,pointer,dimension(:)       :: SVD_cM =>NULL()
  ! Indices into block arrays defined above
    integer(KIND=int_kind),pointer,dimension(:,:)   :: SVDlb  =>NULL()
    integer(KIND=int_kind),pointer,dimension(:,:)   :: SVDub  =>NULL()
    
  ! LSLR "Standard Uncertainty"
    real(KIND=real_kind),pointer,dimension(:,:,:)   :: StandUncert =>NULL()

    logical(KIND=log_kind),pointer,dimension(:,:)   :: done =>NULL()
    
  ! Dirichlet operator
    type (BC_Operator)    ,pointer                   :: DIR_Op =>NULL()
  ! Mask to determine whether worth doing dirichlet loop
  ! for given cell & face
    logical(KIND=log_kind),pointer,dimension(:,:)    :: DIR_Mask =>NULL()
    
  ! Max number of BCs to process; use pointer for consistency.
    integer(KIND=int_kind),pointer                   :: MaxBCLen =>NULL()
    
  ! Max number of face and BC "neighbors" in ncell set; use pointer for consistency.
  ! Calculated in DO_INIT_LSLR_SS, it is used to define size of RHS for LSLR_SVD
    integer(KIND=int_kind),pointer                   :: NumNghbrsMax =>NULL()
    
  ! Variable vector for Neighbor weights
    type(REAL_VAR_VECTOR),pointer,dimension(:)       :: CurNbr_Weights =>NULL()

  ! Block arrays which are the actual allocated objects; arrays defined above
  ! point into these
    integer(KIND=int_kind),pointer,dimension(:)        :: IDXs           =>NULL()
    real(KIND=real_kind),pointer,dimension(:)          :: PHIs           =>NULL()
    real(KIND=real_kind),pointer,dimension(:)          :: GeoW           =>NULL()
    real(KIND=real_kind),pointer,dimension(:)          :: W              =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:)        :: dX             =>NULL()
    integer(KIND=int_kind),pointer,dimension(:)        :: NumFaces       =>NULL()
    integer(KIND=int_kind),pointer,dimension(:)        :: NumFaceNghbrs  =>NULL()
  end type DO_Specifier

  type DO_Diag_Specifier
    integer(KIND=int_kind)                        :: FieldType=-9999

    real(KIND=real_kind),pointer,dimension(:)     :: Phi_Set    =>NULL()

    real(KIND=real_kind),pointer,dimension(:,:)   :: FPhi_Set   =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:,:) :: FGrad_Set  =>NULL()

    real(KIND=real_kind),pointer,dimension(:,:)   :: FPhi_DO    =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:,:) :: FGrad_DO   =>NULL()

    real(KIND=real_kind),pointer,dimension(:,:)   :: FPhi_ERR   =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:)   :: RFPhi_ERR  =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:)   :: FPhi_AERR  =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:,:) :: FGrad_ERR  =>NULL()

    real(KIND=real_kind),pointer,dimension(:)   :: Phi_rmsERR  =>NULL()
    real(KIND=real_kind),pointer,dimension(:)   :: RPhi_rmsERR =>NULL()
    real(KIND=real_kind),pointer,dimension(:)   :: Phi_rmsAERR =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:) :: Grad_rmsERR =>NULL()

    real(KIND=real_kind),pointer,dimension(:,:,:)   :: StandUncert =>NULL()
    real(KIND=real_kind),pointer,dimension(:)       :: SUPhi_rms  =>NULL()
    real(KIND=real_kind),pointer,dimension(:,:)     :: SUGrad_rms =>NULL()
    real(KIND=real_kind),pointer,dimension(:)       :: SUrms      =>NULL()
  end type DO_Diag_Specifier

  type DO_DiagListSpec
    character(LEN=1)  :: field_type      ='F'
    logical(log_kind) :: FacePhiErr      =.false.
    logical(log_kind) :: FaceGradientErr =.false.
    logical(log_kind) :: TaylorErr       =.false.
    logical(log_kind) :: AnalyticErr     =.false.
  end type DO_DiagListSpec

  type cMat_row_type
    real(real_kind),  pointer,dimension(:) :: coeff =>NULL()  ! compress coeff vector
    integer(int_kind),pointer,dimension(:) :: g_idx =>NULL()  ! global index
    integer(int_kind)                      :: n_phi=0         ! num nghbr cells
  end type cMat_row_type


  ! Solution Techniques:
  integer(KIND=int_kind),parameter :: DO_SOLVE_DEFAULT  =0
  integer(KIND=int_kind),parameter :: DO_SOLVE_ORTHO    =1
  integer(KIND=int_kind),parameter :: DO_SOLVE_LU_LSLR  =2
  integer(KIND=int_kind),parameter :: DO_SOLVE_SVD_LSLR =3
  integer(KIND=int_kind),parameter :: DO_NUM_ST         =4  ! **Number of Solution Techniques**


END MODULE DO_BASE_TYPES
