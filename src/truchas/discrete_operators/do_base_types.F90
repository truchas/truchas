!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use kinds, only: r8
  use bc_data_types
  use var_vector_module, only: REAL_VAR_VECTOR
  implicit none
  private

  ! Public Data Structures
  public :: DO_Specifier, &
            SVDU_Type, dX_Type, SField_Type, NIdx_Type

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
      real(r8), pointer, dimension(:,:) :: Mat=>NULL()
    end type SVDU_Type

  ! Structure to hold set of geometric dX vectors
    type dX_Type
      real(r8), pointer, dimension(:,:) :: FData=>NULL()
    end type dX_Type

  ! Structure to hold set of scalar field values
    type SField_Type
      real(r8), pointer, dimension(:) :: FData=>NULL()
    end type SField_Type

  ! Structure to hold index to face neighbors
    type NIdx_Type
      integer, pointer, dimension(:) :: ptr=>NULL()
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
    integer :: ID=-9999
    integer :: Method=-9999

    type(dX_Type)    ,pointer, dimension(:,:) :: dX_Struct   =>NULL()
    type(SField_Type),pointer, dimension(:,:) :: TotW_Struct =>NULL()
    type(SField_Type),pointer, dimension(:,:) :: GeoW_Struct =>NULL()
    type(SField_Type),pointer, dimension(:,:) :: PHI_Struct  =>NULL()
    type(NIdx_Type)  ,pointer, dimension(:)   :: NF_Struct   =>NULL()
    type(NIdx_Type)  ,pointer, dimension(:)   :: NFN_Struct  =>NULL()
    type(NIdx_Type)  ,pointer, dimension(:,:) :: NghIdx      =>NULL()

  ! Ortho Solution
    real(r8), pointer, dimension(:,:,:) :: dX_Scaled     =>NULL()
    real(r8), pointer, dimension(:)     :: W_Ortho       =>NULL()
    real(r8), pointer, dimension(:,:)   :: W_Ortho_Nghbr =>NULL()
    
  ! LSLR Data Structures
    logical, pointer, dimension(:,:)   :: UpdateDecomp =>NULL()
    logical, pointer, dimension(:,:,:) :: SolveFlag    =>NULL()
    
  ! LSLR via LU w/ full pivot
    real(r8), pointer, dimension(:,:,:,:) :: ALU     =>NULL()
    logical,  pointer, dimension(:,:,:)   :: PivFlag =>NULL()
    integer,  pointer, dimension(:,:,:)   :: row1    =>NULL()
    integer,  pointer, dimension(:,:,:)   :: row2    =>NULL()
    
  ! LSLR via SVD
    type(SVDU_Type), pointer, dimension(:) :: SVD_cM =>NULL()
  ! Indices into block arrays defined above
    integer, pointer, dimension(:,:) :: SVDlb  =>NULL()
    integer, pointer, dimension(:,:) :: SVDub  =>NULL()
    
  ! LSLR "Standard Uncertainty"
    real(r8), pointer, dimension(:,:,:) :: StandUncert =>NULL()

    logical, pointer, dimension(:,:) :: done =>NULL()
    
  ! Dirichlet operator
    type(BC_Operator), pointer :: DIR_Op =>NULL()
  ! Mask to determine whether worth doing dirichlet loop
  ! for given cell & face
    logical, pointer, dimension(:,:) :: DIR_Mask =>NULL()
    
  ! Max number of BCs to process; use pointer for consistency.
    integer, pointer :: MaxBCLen =>NULL()
    
  ! Max number of face and BC "neighbors" in ncell set; use pointer for consistency.
  ! Calculated in DO_INIT_LSLR_SS, it is used to define size of RHS for LSLR_SVD
    integer, pointer :: NumNghbrsMax =>NULL()
    
  ! Variable vector for Neighbor weights
    type(REAL_VAR_VECTOR), pointer, dimension(:) :: CurNbr_Weights =>NULL()

  ! Block arrays which are the actual allocated objects; arrays defined above
  ! point into these
    integer,  pointer, dimension(:)   :: IDXs           =>NULL()
    real(r8), pointer, dimension(:)   :: PHIs           =>NULL()
    real(r8), pointer, dimension(:)   :: GeoW           =>NULL()
    real(r8), pointer, dimension(:)   :: W              =>NULL()
    real(r8), pointer, dimension(:,:) :: dX             =>NULL()
    integer,  pointer, dimension(:)   :: NumFaces       =>NULL()
    integer,  pointer, dimension(:)   :: NumFaceNghbrs  =>NULL()
  end type DO_Specifier

  ! Solution Techniques:
  integer, parameter :: DO_SOLVE_DEFAULT  =0
  integer, parameter :: DO_SOLVE_ORTHO    =1
  integer, parameter :: DO_SOLVE_LU_LSLR  =2
  integer, parameter :: DO_SOLVE_SVD_LSLR =3
  integer, parameter :: DO_NUM_ST         =4  ! **Number of Solution Techniques**

END MODULE DO_BASE_TYPES
