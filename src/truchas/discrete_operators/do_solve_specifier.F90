MODULE DO_SOLVE_SPECIFIER
  !=======================================================================
  ! Purpose(s): Instantiate the solve specifier
  !
  !   Public Interfaces: DO_INIT_SS, DO_DESTROY_SS
  !                      DO_GET_cM_compressed, DO_DESTROY_cM_compressed
  !                      DO_GET_cM_full, DO_DESTROY_cM_full
  !
  ! Contains:
  !
  ! Author(s): Doug Kothe (dbk@lanl.gov)
  !            Jeff Durachta (durachta@verizon.net)
  !            Robert Ferrell (ferrell@cpca.com)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  ! Public Subroutines
  public :: DO_INIT_SS
  public :: DO_DESTROY_SS
  public :: DO_GET_cM_compressed, DO_DESTROY_cM_compressed
  public :: DO_GET_cM_full, DO_DESTROY_cM_full

  ! Public Data Structures
  public :: dX_Scaled

  INTERFACE DO_GET_cM_compressed
    MODULE PROCEDURE DO_GET_cM_compressed
  END INTERFACE

  INTERFACE DO_GET_cM_full
    MODULE PROCEDURE DO_GET_cM_full
  END INTERFACE

  INTERFACE DO_DESTROY_cM_compressed
    MODULE PROCEDURE DO_DESTROY_cM_compressed
  END INTERFACE

  INTERFACE DO_DESTROY_cM_full
    MODULE PROCEDURE DO_DESTROY_cM_full
  END INTERFACE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Number of DO_Specifier instantiations
  integer,save      :: NUM_FIELDS=0  ! Initialized to 0

  ! As it's defined, there is only 1 possible dX_Scaled per problem
  ! This coding will allow dX_Scaled to be "destroyed"; this is not a
  ! normal operation, but may be required in particular circumstances
  ! such as operator analysis. 
  ! This object has been made "public" in order that it's status can be 
  ! checked in do_discrete_operators:face_ortho to prevent use if it's
  ! been deallocated. Only other do_* modules are intended to have access
  ! to this variable.
  real(r8),allocatable,target,dimension(:,:,:),save :: dX_Scaled

CONTAINS

  SUBROUTINE DO_INIT_SS(NewSolveSpec,SolveTech,BC_Spec,Weights,GeoExp)
    !=======================================================================
    ! Purpose(s):
    !
    !   Create and initialize the new Solve Specifier data structure.
    !
    !   INPUT:
    !     NewSolveSpec: (Null) pointer to the "to be defined" Solve Specifier
    !     SolveTech   : Solution technique; possible values are:
    !                   DO_SOLVE_ORTHO  :  Use the ORTHO operator
    !                   DO_SOLVE_LU_LSLR:  Use LSLR and the LU solver
    !                   DO_SOLVE_SVD_LSLR: Use LSLR and the SVD solver
    !                   DO_SOLVE_DEFAULT:  Use the value of use_ortho_face_gradient
    !                                      to set either ORTHO or LU_LSLR
    !     BC_Spec     : (optional) Boundary Condition specifier
    !     Weights     : (optional) Non-geomtric weights
    !     GeoExp      : (optional) Geometric distance weighting exponent set to
    !                              1 by default implying the use of 1/(d**2) as
    !                              the geometric weight.
    !
    !   OUTPUT:
    !     NewSolveSpec: Pointer to the Solve Specifier attributes
    !
    !=======================================================================
    use bc_data_types,     only: BC_Specifier
    use do_base_types,     only: DO_Specifier,DO_NUM_ST,DO_SOLVE_ORTHO, &
                                 DO_SOLVE_LU_LSLR,DO_SOLVE_SVD_LSLR,DO_SOLVE_DEFAULT
    use discrete_ops_data, only: use_ortho_face_gradient
    use parameter_module,  only: ncells

    ! Argument list
    type(DO_Specifier), pointer :: NewSolveSpec
    integer, intent(IN) :: SolveTech
    type(BC_Specifier), optional, target, intent(IN) :: BC_Spec
    real(r8), optional, dimension(ncells), intent(IN) :: Weights
    integer, optional, intent(IN) :: GeoExp

    ! Local Variables
    integer :: ST
    integer :: istat

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(SolveTech<0 .or. SolveTech>DO_NUM_ST)then
      call TLS_panic ('DO_INIT_SS: attempt to use undefined SolveTech.')
    endif

    if(ASSOCIATED(NewSolveSpec))then
      call TLS_panic ('DO_INIT_SS: attempt to re-initialize associated DO_Specifier')
    else
    ! All solve specifier pointer structure components are set to NULL()
    ! in the DO_Specifier type definition 
      ALLOCATE(NewSolveSpec,STAT=istat)
      if (istat /= 0) call TLS_panic ('DO_INIT_SS: memory allocation failure for SolveSpec')
    endif

    ! Number of Solve Specifiers
    NUM_FIELDS = NUM_FIELDS + 1
    NewSolveSpec%ID = NUM_FIELDS

    ST = SolveTech
    if(ST==DO_SOLVE_DEFAULT)then
      if(use_ortho_face_gradient)then
        ST=DO_SOLVE_ORTHO
      else
        ST = DO_SOLVE_LU_LSLR
      endif
    endif
    NewSolveSpec%Method = ST

    METHOD: SELECT CASE(ST)

      case(DO_SOLVE_ORTHO)
        call DO_INIT_ORTHO(SS=NewSolveSpec,Weights=Weights,GeoExp=GeoExp)
      case(DO_SOLVE_LU_LSLR)
        call DO_INIT_LU_LSLR(SS=NewSolveSpec,BC_Spec=BC_Spec,Weights=Weights,GeoExp=GeoExp)
      case(DO_SOLVE_SVD_LSLR)
        call DO_INIT_SVD_LSLR(SS=NewSolveSpec,BC_Spec=BC_Spec,Weights=Weights,GeoExp=GeoExp)
      case DEFAULT
        call TLS_fatal ('DO_INIT_SS: attempt to use unknown SolveTech')
    END SELECT METHOD

  END SUBROUTINE DO_INIT_SS


  SUBROUTINE DO_INIT_ORTHO(SS,Weights,GeoExp)
    !=======================================================================
    ! Purpose(s):
    !
    !   Create and initialize a new ORTHO Solve Specifier data structure.
    !
    !   INPUT:
    !     NewSolveSpec: (Null) pointer to the "to be defined" Solve Specifier
    !     Weights     : (optional) Non-geomtric weights
    !     GeoExp      : (optional) Geometric distance weighting exponent set to
    !                              1 by default implying the use of 1/(d**2) as
    !                              the geometric weight.
    !
    !   OUTPUT:
    !     NewSolveSpec: Pointer to the Solve Specifier attributes
    !
    !=======================================================================
    use cutoffs_module,   only: alittle
    use do_base_types,    only: DO_Specifier
    use gs_module,        only: EE_GATHER
    use mesh_module,      only: Cell, Mesh
    use parameter_module, only: ncells, nfc, ndim

    ! Argument list
    type(DO_Specifier), target, intent(INOUT) :: SS
    real(r8), optional, dimension(ncells), intent(IN) :: Weights
    integer, optional, intent(IN) :: GeoExp

    ! Local Variables
    integer :: n, f, f2, j, j2, istat
    real(r8), dimension(ndim) :: dX
    real(r8) :: dX_tmp
    real(r8), dimension(ndim,nfc,ncells) :: X_e
    real(r8) :: Distance
    real(r8) :: pval=10.d0**(- precision(dX_tmp))
    integer, save :: GExp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(.not. ALLOCATED(dX_Scaled))then
      GExp = 1; if(PRESENT(GeoExp))GExp = GeoExp

    ! Gather coordinates and Phi of face neighbors.

      ALLOCATE(dX_Scaled(ndim,nfc,ncells),STAT=istat)
      if (istat /= 0) call TLS_panic ('DO_INIT_ORTHO: memory allocation failure for dX_Scaled')
      SS%dX_Scaled =>dX_Scaled  ! There is only one dX_Scaled

      do n = 1,ndim
         call EE_GATHER (DEST=X_e(n,:,:), SRC=Cell%Centroid(n))
      end do
      CELL_LOOP: do j=1,ncells
      FACE_LOOP: do f = 1,nfc
         ! Compute Deltas (Cell to Face Neighbors)
         ! Physical Coordinate Deltas
         dX       = 0.0_r8
         Distance = 0.0_r8
         do n = 1,ndim
            dX(n)  = Cell(j)%Centroid(n) - X_e(n,f,j)
            if(abs(dX(n)) < pval)dX(n)=0.0_r8
            Distance = Distance + dX(n)**2
         end do
         Distance = Distance**GExp
         if(Mesh(j)%Ngbr_Cell(f) == 0)then
           dX_scaled(:,f,j) = 0.0_r8
         else
           do n = 1,ndim
             dX_tmp = dX(n)/(Distance + alittle)
             dX_scaled(n,f,j) = dX_tmp
           end do
         end if
      end do FACE_LOOP
      end do CELL_LOOP
    else
    ! As it's defined, there is only 1 possible dX_Scaled per problem
      SS%dX_Scaled =>dX_Scaled
    endif

    ALLOCATE(SS%done(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for DONE')
    ! "Done" vector allows calculation of shared face only once so long as
    ! neighboring cell centers are "on processor".
    SS%done = .false.
    CELL_LOOP00: do j = 1,ncells
      ! Analyze face connectivity
      FNGHBR_LOOP: do f = 1, nfc
        if(Mesh(j)%Ngbr_Cell(f) > 0 .and. .not. SS%done(f,j)) then
          j2 = Mesh(j)%Ngbr_Cell(f)
          f2 = Mesh(j)%Ngbr_Face(f)
          SS%done(f2,j2) = .true.
        endif
      end do FNGHBR_LOOP
    end do CELL_LOOP00

    if(PRESENT(Weights))then           ! Weights may vary between fields
      ALLOCATE(SS%W_Ortho(ncells),STAT=istat)
      if (istat /= 0) call TLS_panic ('DO_INIT_ORTHO: memory allocation failure for W_Ortho')
      SS%W_Ortho = 0.0_r8
      where(Weights > 0.0_r8)SS%W_Ortho = 1.0_r8

      ALLOCATE(SS%W_Ortho_Nghbr(nfc,ncells),STAT=istat)
      if (istat /= 0) call TLS_panic ('DO_INIT_ORTHO: memory allocation failure for W_Ortho_Nghbr')

      call EE_GATHER (DEST=SS%W_Ortho_Nghbr, SRC=SS%W_Ortho)
    endif  ! if WEIGHTS
    ALLOCATE(SS%MaxBCLen,STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_ORTHO: memory allocation failure for MaxBCLen')
    SS%MaxBCLen = 0  ! BCs are not currently used with ORTHO solution
  END SUBROUTINE DO_INIT_ORTHO


  SUBROUTINE DO_INIT_LU_LSLR(SS,BC_Spec,Weights,GeoExp)
    !=======================================================================
    ! Purpose(s):
    !
    !   Create and initialize the new Solve Specifier data structure.
    !
    !   INPUT:
    !     NewSolveSpec: (Null) pointer to the "to be defined" Solve Specifier
    !     BC_Spec     : (optional) Boundary Condition specifier
    !     Weights     : (optional) Non-geomtric weights
    !     GeoExp      : (optional) Geometric distance weighting exponent set to
    !                              1 by default implying the use of 1/(d**2) as
    !                              the geometric weight.
    !
    !   OUTPUT:
    !     NewSolveSpec: Pointer to the Solve Specifier attributes
    !
    !=======================================================================
    use bc_data_types,    only: BC_Specifier
    use do_base_types,    only: DO_Specifier
    use parameter_module, only: ncells, nfc, ndim

    ! Argument list
    type(DO_Specifier), target, intent(INOUT) :: SS
    type(BC_Specifier), target, optional, intent(IN) :: BC_Spec
    real(r8), optional, dimension(ncells), intent(IN) :: Weights
    integer, optional, intent(IN) :: GeoExp

    ! local variables
    integer :: istat
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Create and initialize common aspects of LSLR solve specifier
    call DO_INIT_LSLR_SS(SS=SS,BC_Spec=BC_Spec,Weights=Weights,GeoExp=GeoExp)

    ! Allocate LU unique SS components
    ALLOCATE(SS%ALU(ndim+1,ndim+1,nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_LU_LSLR: memory allocation failure for ALU')

    ALLOCATE(SS%PivFlag(ndim+1,nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_LU_LSLR: memory allocation failure for PivFlag')
    ALLOCATE(SS%row1(ndim+1,nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_LU_LSLR: memory allocation failure for row1')
    ALLOCATE(SS%row2(ndim+1,nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_LU_LSLR: memory allocation failure for row2')
    ! Initialize LU unique SS components
    call LU_Init(SS)
  END SUBROUTINE DO_INIT_LU_LSLR


  SUBROUTINE DO_INIT_SVD_LSLR(SS,BC_Spec,Weights,GeoExp)
    !=======================================================================
    ! Purpose(s):
    !
    !   Create and initialize the new Solve Specifier data structure.
    !
    !   INPUT:
    !     NewSolveSpec: (Null) pointer to the "to be defined" Solve Specifier
    !     BC_Spec     : (optional) Boundary Condition specifier
    !     Weights     : (optional) Non-geomtric weights
    !     GeoExp      : (optional) Geometric distance weighting exponent set to
    !                              1 by default implying the use of 1/(d**2) as
    !                              the geometric weight.
    !
    !   OUTPUT:
    !     NewSolveSpec: Pointer to the Solve Specifier attributes
    !
    !=======================================================================
    use bc_data_types,    only: BC_Specifier
    use do_base_types,    only: DO_Specifier
    use parameter_module, only: ncells, nfc, ndim

    ! Argument list
    type(DO_Specifier), target, intent(INOUT) :: SS
    type(BC_Specifier), optional, target, intent(IN) :: BC_Spec
    real(r8), optional, dimension(ncells), intent(IN) :: Weights
    integer, optional, intent(IN) :: GeoExp

    ! local variables
    integer, pointer, dimension(:) :: NumFacesCell => NULL()
    integer :: NFC_Tot
    integer :: c, f, istat

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Create and initialize common aspects of LSLR solve specifier
    call DO_INIT_LSLR_SS(SS=SS,BC_Spec=BC_Spec,Weights=Weights,GeoExp=GeoExp)

    ! Allocate SVD unique SS components
    ALLOCATE(SS%SVD_cM(ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_SVD_LSLR: memory allocation failure for SVD_U')
    ALLOCATE(SS%SVDlb(nfc,ncells),SS%SVDub(nfc,ncells), STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_SVD_LSLR: memory allocation failure for lb or ub')
    do c = 1,ncells
    ! Find number of nghbr & bc faces affecting this face
      NumFacesCell => SS%NF_Struct(c)%ptr
    ! Find sum over NFC faces
      NFC_Tot = sum(NumFacesCell(:))
    ! Allocate total and later point into array to reduce memory splintering
    ! and possibly enhance cache loading
      ALLOCATE(SS%SVD_cM(c)%Mat(ndim+1,NFC_Tot),STAT=istat)
      if (istat /= 0) call TLS_panic ('DO_INIT_SVD_LSLR: memory allocation failure for Cmatrix')
    ! Initialize index pointers into SVD_cM(c1)%Mat array
      do f=1,nfc
        SS%SVDub(f,c) = sum(NumFacesCell(1:f))
        SS%SVDlb(f,c) = SS%SVDub(f,c)-NumFacesCell(f)+1
      end do
    end do

    ! Initialize SVD unique SS components
    call SVD_Init(SS)
  END SUBROUTINE DO_INIT_SVD_LSLR


  SUBROUTINE DO_INIT_LSLR_SS(SS,BC_Spec,Weights,GeoExp)
    !=======================================================================
    ! Purpose(s): To instantiate and initialize the persistent data
    !             structures for face neighbor and BC computations
    !
    !             This routine is "simply" a "game" of counting the number
    !             of valid face neighbors and BCs and packing the
    !             persistent data (dXs, geometric weights, and indices into
    !             the list of valid face neighbor points) into
    !             arrays.
    !
    !=======================================================================
    use bc_data_types,     only: BC_Specifier,BC_Chart_ID,BC_Spec_Get_Operator, &
                                 BC_DIRICHLET_Op,BC_Op_Start_Search,BC_Get_Length, &
                                 BC_Get_Atlas,BC_Get_Chart,BC_Chart_Length, &
                                 BC_Chart_Positions, &
                                 DIMENSIONALITY
    use cutoffs_module,    only: alittle
    use do_base_types,     only: DO_Specifier,SField_Type
    use gs_module,         only: EE_GATHER
    use mesh_module,       only: Cell,Mesh,Is_Face_Ngbr
    use parameter_module,  only: ncells, nfc, ndim
    use var_vector_module, only: REAL_VAR_VECTOR, CREATE, SIZES, FLATTEN, DESTROY

    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS
    type(BC_Specifier),  OPTIONAL, target, intent(IN) :: BC_Spec
    real(r8), OPTIONAL, dimension(ncells), intent(IN) :: Weights
    integer, OPTIONAL, intent(IN) :: GeoExp

    ! Local Variables
  ! Various index variables
    integer :: d,j,n,f,j1,f1,j2,f2
  ! Running count of the number of valid neighbor faces
    integer :: NumNghbrs
    integer, dimension(ncells) :: NumNghbrsAll
  ! The ...tmp arrays are temporary structures to allocated to hold the
  ! maximum possible extent of the required data until the true extent is
  ! determined by Is_Face_Ngbr, DEGENERATE_FACE and the BC routines
    integer, allocatable, dimension(:,:,:) :: IDXtmp

    type(REAL_VAR_VECTOR), pointer, dimension(:,:) :: X_Centers => NULL()

    type(SField_Type), dimension(ndim) :: Xc
    integer, pointer, dimension(:) :: Ngbr_Cells_Face => NULL()
    integer :: FN
    type (BC_Chart_ID), pointer :: ChartID => NULL()
    integer :: Length
    real(r8), dimension(:,:), pointer :: Positions => NULL()
    
  ! This is for all the neighbor weights
    real(r8), dimension(:), pointer :: CellsNbrWeights => NULL()
    
  ! This is for all the neighbor and BC faces
    integer, pointer, dimension(:) :: NumFacesCell => NULL()
    
  ! This is for all the neighbor faces
    integer, pointer, dimension(:) :: NumNghbrCell => NULL()

    real(r8) :: dX1
    integer :: lb,ub,idx,istat,NFC_offset,NFC_F,j_offset
    integer :: NumComputedF
    integer :: GExp
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call DO_CREATE_SS(SS)

    SS%UpdateDecomp=.false.
    SS%done = .false.

    ! Allocate a vector to store neighboring coordinate data.
    ALLOCATE (X_Centers(ndim,ncells))
    GATHER_LOOP: do d = 1,ndim
      call CREATE (ARRAY = X_Centers(d,:),                &
                   SIZES = SIZES(Mesh(1:ncells)%Ngbr_Cells_All) )
      call EE_GATHER (DEST=X_Centers(d,:), SOURCE=Cell%Centroid(d))
    end do GATHER_LOOP
                      
    if(PRESENT(Weights))then
      ALLOCATE(SS%CurNbr_Weights(ncells),STAT=istat)
      if (istat /= 0) call TLS_panic ('DO_INIT_LSLR_SS: memory allocation failure for CurNbr_Weights')
      ! Gather the weights from all neighbors, for use by all cells
      call CREATE(ARRAY = SS%CurNbr_Weights(:), &
                  SIZES = SIZES(Mesh(1:ncells)%Ngbr_Cells_All))
      call EE_Gather (DEST=SS%CurNbr_Weights, SOURCE=Weights)
    endif
      
    ! Get the DIRICHLET operators
    if(PRESENT(BC_Spec))then                             
      SS%DIR_Op => BC_Spec_Get_Operator(BC_Spec, BC_DIRICHLET_Op)
      call BC_Op_Start_Search(SS%DIR_Op)
    ! Note that BC_Get_Length returns nothing (NULL) if there are no values
    ! associated with that operator. By design, the intrinsic MAXVAL must return
    ! something... which happens to be the negative number of greatest magnitude
    ! consistent with KIND of the calling LHS. Hence the use of the array
    ! constructor to make 0 the smallest value which may be returned.
      SS%MaxBCLen = MAXVAL( (/0,BC_Get_Length(BC_Get_Atlas(SS%DIR_Op))/) )
    else
      SS%MaxBCLen = 0
    endif
      
  ! The calls to the BC_Op functions are relativly expensive. More important
  ! the results cannot be "prefetched". In principle, the compiler can arrange
  ! to "look ahead" in the mask arrays.
    SS%DIR_Mask = .false.                                   
    
    ALLOCATE(SS%NumFaces(nfc*ncells), SS%NumFaceNghbrs(nfc*ncells), STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_LSLR_SS: memory allocation failure for NumFace arrays')
    SS%NumFaces = 0                                         
    SS%NumFaceNghbrs = 0                                    
    NumComputedF = 0   ! set num computed faces neighbors = 0
    NumNghbrsAll = SIZES(X_Centers(1,:))
    SS%NumNghbrsMax = maxval(NumNghbrsAll(:)) + SS%MaxBCLen
    ALLOCATE( IDXtmp(SS%NumNghbrsMax,nfc,ncells), STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_LSLR_SS: memory allocation failure for IDXtmp')
    IDXtmp = -HUGE(IDXtmp)  ! init value almost certain to cause seg vio if used

    ! Analyze face connectivity; done created first for simplicity
    CELL_LOOP00: do j = 1,ncells
      FNGHBR_LOOP: do f = 1, nfc
        if(Mesh(j)%Ngbr_Cell(f) > 0 .and. .not. SS%done(f,j)) then
          j2 = Mesh(j)%Ngbr_Cell(f)
          f2 = Mesh(j)%Ngbr_Face(f)
          SS%done(f2,j2) = .true.
        endif
      end do FNGHBR_LOOP
    end do CELL_LOOP00

    ! Count the number of shared faces
    CELL_LOOP0: do j = 1,ncells
      j_offset = (j-1)*nfc
      SS%NF_Struct(j)%ptr  =>SS%NumFaces(1+j_offset:nfc+j_offset)
      SS%NFN_Struct(j)%ptr =>SS%NumFaceNghbrs(1+j_offset:nfc+j_offset)
      NumNghbrs = NumNghbrsAll(j)
      Ngbr_Cells_Face => FLATTEN(Mesh(j)%Ngbr_Cells_Face)
      FACE_LOOP0: do f = 1, nfc
        if(SS%done(f,j))then
        ! Enforce consistency of shared faces
        ! Find (face,cell) which actually performs calculation
          j1 = Mesh(j)%Ngbr_Cell(f)
          f1 = Mesh(j)%Ngbr_Face(f)
          SS%NFN_Struct(j)%ptr(f) = SS%NFN_Struct(j1)%ptr(f1)
          SS%NF_Struct(j)%ptr(f) = SS%NF_Struct(j1)%ptr(f1)
        else
          NEIGHBOR_LOOP0: do n=1,NumNghbrs
          ! If cell is not a boundary cell or a neighbor cell, cycle
            if(Mesh(j)%Ngbr_Cell(f)==0 .or. Is_Face_Ngbr(Ngbr_Cells_Face(n),f) )then
              SS%NFN_Struct(j)%ptr(f) = SS%NFN_Struct(j)%ptr(f) + 1  ! Increment num face neighbors
              SS%NF_Struct(j)%ptr(f) = SS%NF_Struct(j)%ptr(f) + 1    ! Increment num total faces
              FN = SS%NF_Struct(j)%ptr(f)
              IDXtmp(FN,f,j) = n  ! Track index of face neighbor
            endif  ! if(Mesh(j)%Ngbr_Cell(f)==0
          end do NEIGHBOR_LOOP0

        ! Now add the BC points
          if(SS%MaxBCLen>0)then
          ! Check to see if (f, c) is in the bdy operator
            if(BC_Get_Chart  &
                       (CHARTID=ChartID,OPERATOR=SS%DIR_Op,CELL=j,FACE=f))then
              SS%DIR_Mask(f,j) = .true.
              SS%NF_Struct(j)%ptr(f) = SS%NF_Struct(j)%ptr(f) + BC_Chart_Length(ChartID)  ! Increment num total faces
            endif  ! if(BC_Get_Chart
          endif  ! if(MaxBCLen>0)
        ! increment the number of face components actually computed
          NumComputedF = NumComputedF + SS%NF_Struct(j)%ptr(f)
        endif  ! if(SS%done(f,j))
      end do FACE_LOOP0
    end do CELL_LOOP0

  ! Now that we know how many valid items there really are, allocate the
  ! persistent storage and point to it; Note that IDXs really only needs
  ! to be the size of number of face neighbors whereas everything else is 
  ! the total number of faces in the calcuation. It turns out that it simplifies
  ! the coding and counting to make everything the same size.
    SS%NumNghbrsMax = (maxval(SS%NumFaces))  ! Set actual maximum number of FN + BC
    ALLOCATE(SS%dX(ndim+1,NumComputedF), &
                     SS%W(NumComputedF), &
                  SS%GeoW(NumComputedF), &
                  SS%PHIs(NumComputedF), &
                  SS%IDXs(NumComputedF), STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_INIT_LSLR_SS: memory allocation failure for dX,W,GeoW and PHI arrays')
    SS%IDXs = -HUGE(SS%IDXs)  ! init value almost certain to cause seg vio if used
    GExp=1; if(PRESENT(GeoExp))GExp=GeoExp
    SS%W = 0.0_r8; SS%GeoW = 0.0_r8
    idx = 0
    NFC_offset = 0
    CELL_LOOP: do j = 1,ncells
      do d = 1, ndim
        Xc(d)%FData => FLATTEN(X_Centers(d,j))
      end do
      if(PRESENT(Weights))CellsNbrWeights => FLATTEN(SS%CurNbr_Weights(j))
      NumFacesCell => SS%NF_Struct(j)%ptr
      NumNghbrCell => SS%NFN_Struct(j)%ptr
      FACE_LOOP: do f = 1, nfc
        if(SS%done(f,j))then ! Only calculate unique faces
          j1 = Mesh(j)%Ngbr_Cell(f)
          f1 = Mesh(j)%Ngbr_Face(f)
          SS%dX_Struct(f,j)%FData   =>  SS%dX_Struct(f1,j1)%FData
          SS%GeoW_Struct(f,j)%FData =>SS%GeoW_Struct(f1,j1)%FData
          SS%TotW_Struct(f,j)%FData =>SS%TotW_Struct(f1,j1)%FData
          SS%Phi_Struct(f,j)%FData  => SS%Phi_Struct(f1,j1)%FData
          SS%NghIdx(f,j)%ptr        =>     SS%NghIdx(f1,j1)%ptr
        else
          lb = idx+1; ub = idx+NumNghbrCell(f)
          if(ub>NumComputedF) &
            call TLS_panic ('DO_INIT_LSLR_SS: memory indexing error for IDXs array')
        ! IDXs vector stores index of face neighbor Phi
        ! BC Phi values are acquired from BC routines
          SS%IDXs(lb:ub) = IDXtmp(1:NumNghbrCell(f),f,j)
          SS%NghIdx(f,j)%ptr => SS%IDXs(lb:ub)  ! Point to this face's indices
          NEIGHBOR_LOOP: do n=1,NumNghbrCell(f)
             idx = idx + 1   ! increment the index into the persistent store
             FN = SS%IDXs(idx)
           ! Compute the distance vector and weight for this neighbor.
             do d = 1, ndim
                SS%dX(d+1,idx) = Xc(d)%FData(FN) - Cell(j)%Face_Centroid(d,f)
                SS%GeoW(idx)   = SS%GeoW(idx) + SS%dX(d+1,idx)*SS%dX(d+1,idx)
             end do
           ! Use distance to scale phi solution. This keeps problem matrices scale invariant
             SS%GeoW(idx) = SS%GeoW(idx)**GExp
             SS%dX(1,idx) = sqrt(SS%GeoW(idx))
           ! Geometric weight is 1/d**2, where d is the distance
           ! to the neighbor.
             if(SS%GeoW(idx) >= alittle) SS%GeoW(idx) = 1.0_r8/SS%GeoW(idx)
             SS%W(idx) = SS%GeoW(idx)
             if(PRESENT(Weights))SS%W(idx) = SS%W(idx)*CellsNbrWeights(FN)
          end do NEIGHBOR_LOOP

         ! Now add the BC points
          if(SS%MaxBCLen>0)then
            if(j==1)call BC_Op_Start_Search(SS%DIR_Op)
            ! Check to see if (f, c) is in the bdy operator
            if(BC_Get_Chart  &
                       (CHARTID=ChartID,OPERATOR=SS%DIR_Op,CELL=j,FACE=f))then
              SS%DIR_Mask(f,j) = .true.
              ! Found a non-empty chart, so grab the data from it.
              Length                = BC_Chart_Length(ChartID)
              ALLOCATE(Positions(DIMENSIONALITY(SS%DIR_Op),Length))
              Positions(:,1:Length) = BC_Chart_Positions(ChartID)
              BC_DIR_LOOP: do n = 1,Length
                idx = idx + 1   ! increment the index into the persistent store
                ! Compute the distance vector and weight for this neighbor.
                do d = 1, ndim
                  SS%dX(d+1,idx) = Positions(d,n) - Cell(j)%Face_Centroid(d,f)
                  SS%GeoW(idx)   = SS%GeoW(idx) + SS%dX(d+1,idx)*SS%dX(d+1,idx)
                end do
                SS%GeoW(idx) = SS%GeoW(idx)**GExp
                SS%dX(1,idx) = sqrt(SS%GeoW(idx))  ! Use distance to scale phi solution
                ! Geometric weight is 1/d**2, where d is the distance
                ! to the neighbor.
                if(SS%GeoW(idx) >= alittle) SS%GeoW(idx) = 1.0_r8/SS%GeoW(idx)
                SS%W(idx) = SS%GeoW(idx)  ! No Non-Geo BC weights at this time
              end do BC_DIR_LOOP
              DEALLOCATE(Positions)
            endif  ! if(BC_Get_Chart
          endif  ! if(MaxBCLen>0)
          ub = idx  ! New upper bound includes BCs; lb same as neighbor loop
          if(ub>NumComputedF)then
             call TLS_panic ('DO_INIT_LSLR_SS: memory indexing error for dX, GeoW and other arrays')
          endif
          NFC_F = NumFacesCell(f)
          if(NFC_F /=(ub-lb+1))then
            call TLS_panic ('DO_INIT_LSLR_SS: memory bounds error for dX, GeoW and other arrays')
          endif
          if (NFC_F > 1) then
            dX1 = sum(SS%dX(1,lb:ub))
            SS%dX(1,lb:ub) = dX1 / float(NFC_F)  ! Average the "scale invariant" dX1 components
          end if
          SS%dX_Struct(f,j)%FData => SS%dX(:,lb:ub)
          SS%GeoW_Struct(f,j)%FData => SS%GeoW(lb:ub)
          SS%TotW_Struct(f,j)%FData => SS%W(lb:ub)
          SS%Phi_Struct(f,j)%FData => SS%PHIs(lb:ub)
        endif  ! if(SS%done(f,j))
      end do FACE_LOOP
    end do CELL_LOOP

    DEALLOCATE(IDXtmp)
  ! Destroy and Deallocate the X_Centers working array.
    do d = 1,ndim
       call DESTROY (ARRAY = X_Centers(d,:))
    end do
    DEALLOCATE (X_Centers)
  END SUBROUTINE DO_INIT_LSLR_SS


  SUBROUTINE DO_CREATE_SS(SS)
    !=======================================================================
    ! Purpose(s): To instantiate the persistent data structures in the
    !             discrete operators container (doc) for 
    !             face neighbor and BC computations
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier
    use parameter_module, only: ncells, nfc, ndim

    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS

    ! Local Variables
    integer :: istat
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE(SS%dX_Struct(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array dX_Struct')
    ALLOCATE(SS%PHI_Struct(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array PHI_Struct')
    ALLOCATE(SS%TotW_Struct(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array W_Struct')
    ALLOCATE(SS%GeoW_Struct(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array GeoW_Struct')
    ALLOCATE(SS%NghIdx(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array NghIdx')
    ALLOCATE(SS%NF_Struct(ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array NF_Struct')
    ALLOCATE(SS%NFN_Struct(ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array NFN_Struct')
    ALLOCATE(SS%UpdateDecomp(nfc,ncells))
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array UpdateDecomp')
    ALLOCATE(SS%SolveFlag(ndim+1,nfc,ncells))
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array SolveFlag')
    ALLOCATE(SS%StandUncert(ndim+1,nfc,ncells))
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array StandUncert')
    ALLOCATE(SS%done(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array DONE')
    SS%DIR_Op => NULL()  ! Set here for completeness. Dir_Op is set to NULL()
                         ! in DO_Specifier typedef. It is created by BC routines

    ALLOCATE(SS%DIR_Mask(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array DIR_Mask')
    ALLOCATE(SS%MaxBCLen,STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array MaxBCLen')
    ALLOCATE(SS%NumNghbrsMax,STAT=istat)
    if (istat /= 0) call TLS_panic ('DO_CREATE_SS: memory allocation failure for array NumNghbrsMax')
    SS%CurNbr_Weights => NULL()  ! Set here for completeness. CurNbr_Weights is set
                                 ! to NULL() in DO_Specifier typedef. It is
                                 ! instantiated in DO_INIT_LSLR_SS if weights present

  END SUBROUTINE DO_CREATE_SS


  SUBROUTINE DO_DESTROY_SS(SS,rm_dX_Scaled)
    !=======================================================================
    ! Purpose(s): To destroy the persistent data structures in the
    !             discrete operators container (doc) for
    !             face neighbor and BC computations
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier
    use parameter_module,  only: ncells
    use var_vector_module, only: DESTROY

    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS
    logical, optional, intent(IN) :: rm_dX_Scaled

    ! Local Variables
    integer :: i
    logical :: rmdXScaled
    
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Deallocate block storage; Nullify (Face,cell) pointers

    if(ASSOCIATED(SS%dX_Struct))then
      DEALLOCATE(SS%NumFaces, SS%NumFaceNghbrs)
      DEALLOCATE(SS%dX, SS%W, SS%GeoW, SS%PHIs, SS%IDXs)

      DEALLOCATE(SS%dX_Struct);   SS%dX_Struct   =>NULL()
      DEALLOCATE(SS%TotW_Struct); SS%TotW_Struct =>NULL()
      DEALLOCATE(SS%GeoW_Struct); SS%GeoW_Struct =>NULL()
      DEALLOCATE(SS%PHI_Struct);  SS%PHI_Struct  =>NULL()
      DEALLOCATE(SS%NF_Struct);   SS%NF_Struct   =>NULL()
      DEALLOCATE(SS%NFN_Struct);  SS%NFN_Struct  =>NULL()
      DEALLOCATE(SS%NghIdx);      SS%NghIdx      =>NULL()
    endif

  ! Ortho Solution

  ! dX_Scaled is a "shared" object which may have multiple pointers to it.
  ! Still, one must allow for it to be destroyed (e.g. in the case of diagnostic
  ! analysis). Unfortunately, conflicting requirements of "permenance" vs "destroyable"
  ! create possible problems; Allocation status of dX_Scaled checked in face_ortho of
  ! do_discrete_operators.
  ! We allow deallocation ONLY if in use by this instantiation of SS in the first place

  ! Do we want to destroy it in the first place?
    rmdXScaled = .false.; if(PRESENT(rm_dX_Scaled))rmdXScaled=rm_dX_Scaled
    if(rmdXScaled .and. ASSOCIATED(SS%dX_Scaled))then  ! Was this SS using dX_Scaled at all?
    ! Will test status of this in do_discrete_operators::face_ortho to ensure that there's
    ! no attempt to use an unallocated dX_Scaled
      if(ALLOCATED(dX_Scaled))DEALLOCATE(dX_Scaled)
    endif
                                                                  SS%dX_Scaled     =>NULL()
    if(ASSOCIATED(SS%W_Ortho))DEALLOCATE(SS%W_Ortho);             SS%W_Ortho       =>NULL()
    if(ASSOCIATED(SS%W_Ortho_Nghbr))DEALLOCATE(SS%W_Ortho_Nghbr); SS%W_Ortho_Nghbr =>NULL()

  ! LSLR Data Structures
    if(ASSOCIATED(SS%UpdateDecomp))DEALLOCATE(SS%UpdateDecomp); SS%UpdateDecomp =>NULL()
    if(ASSOCIATED(SS%SolveFlag))DEALLOCATE(SS%SolveFlag);       SS%SolveFlag    =>NULL()

  ! LSLR via LU w/ full pivot
    if(ASSOCIATED(SS%ALU))DEALLOCATE(SS%ALU);         SS%ALU     =>NULL()
    if(ASSOCIATED(SS%PivFlag))DEALLOCATE(SS%PivFlag); SS%PivFlag =>NULL()
    if(ASSOCIATED(SS%row1))DEALLOCATE(SS%row1);       SS%row1    =>NULL()
    if(ASSOCIATED(SS%row2))DEALLOCATE(SS%row2);       SS%row2    =>NULL()

  ! LSLR via SVD
    if(ASSOCIATED(SS%SVD_cM))then
      do i = 1,ncells
        DEALLOCATE(SS%SVD_cM(i)%Mat)
      end do
      DEALLOCATE(SS%SVD_cM); SS%SVD_cM =>NULL()
    endif

  ! LSLR "Standard Uncertainty"
    if(ASSOCIATED(SS%StandUncert))DEALLOCATE(SS%StandUncert); SS%StandUncert =>NULL()

    if(ASSOCIATED(SS%done))DEALLOCATE(SS%done); SS%done =>NULL()

  ! Dirichlet operator
  ! SS%DIR_Op points to BC operator; SS%DIR_Op cannot be deallocated
                                                        SS%DIR_Op   =>NULL()
    if(ASSOCIATED(SS%DIR_Mask))DEALLOCATE(SS%DIR_Mask); SS%DIR_Mask =>NULL()

  ! Max number of BCs to process; use pointer for consistency.
    if(ASSOCIATED(SS%MaxBCLen))DEALLOCATE(SS%MaxBCLen); SS%MaxBCLen =>NULL()

  ! Max number of face and BC "neighbors" in ncell set; use pointer for consistency.
  ! Calculated in DO_INIT_LSLR_SS, it is used to define size of U for LSLR_SVD
    if(ASSOCIATED(SS%NumNghbrsMax))DEALLOCATE(SS%NumNghbrsMax); SS%NumNghbrsMax =>NULL()

  ! Variable vector for Neighbor weights
  ! CurNbrWeights is of type(REAL_VAR_VECTOR) and must be deconstructed
  ! by its native module
    if(ASSOCIATED(SS%CurNbr_Weights))call DESTROY(SS%CurNbr_Weights); SS%CurNbr_Weights =>NULL()

    if(ASSOCIATED(SS%SVDlb))then
      DEALLOCATE(SS%SVDlb,SS%SVDub); SS%SVDlb =>NULL(); SS%SVDub =>NULL()
    endif
  END SUBROUTINE DO_DESTROY_SS


  SUBROUTINE LU_Init(SS)
    !=======================================================================
    ! Purpose(s): Initializes LU decomposition
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier
    use do_update_module, only: UpdateFaceLU
    use parameter_module, only: ncells, nfc

    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS
    ! Local Variables
    integer :: Face,c1

    do c1 = 1,ncells
      FACE_LOOP: do Face=1,nfc
        if(SS%done(Face,c1))cycle FACE_LOOP
        call UpdateFaceLU(SS,Face,c1)
      end do FACE_LOOP
    end do
  END SUBROUTINE LU_Init


  SUBROUTINE SVD_Init(SS)
    !=======================================================================
    ! Purpose(s): Initializes SVD decomposition
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier
    use do_update_module, only: UpdateFaceSVD
    use parameter_module, only: ncells, nfc
 
    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS
    ! Local Variables
    integer :: Face,c1

    do c1 = 1,ncells
      FACE_LOOP: do Face=1,nfc
        if(SS%done(Face,c1))cycle FACE_LOOP
        call UpdateFaceSVD(SS,Face,c1)
      end do FACE_LOOP
    end do
  END SUBROUTINE SVD_Init


  FUNCTION DO_get_cM_compressed(SS,do_BCs) RESULT(cM)
    !=======================================================================
    ! Purpose(s): To instantiate a compressed form of the "coefficient
    !             matrix" described in 'Discussion of the T Diffusion
    !             Discretization', M. Hall, 11/06/02
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier,cMat_row_type,DO_SOLVE_SVD_LSLR
    use mesh_module,      only: Cell,Mesh,DEGENERATE_FACE
    use parallel_info_module, only: p_info
    use parameter_module,  only: ndim,nfc,ncells

    type(DO_Specifier), target, intent(IN) :: SS
    logical, optional, intent(IN) :: do_BCs

    type(cMat_row_type), pointer, dimension(:) :: cM
    real(r8) :: cM_tmp

    logical :: doBCs
    integer :: istat, idx, idx2

  ! SVD solution on "Design Matrix"
    real(r8), pointer, dimension(:,:) :: DO_cM => NULL()

    integer :: c,c2,f,f2,n,m,NFN,lb,ub
    integer,pointer,dimension(:) :: NumFacesCell => NULL()
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    cM => NULL()
    if(p_info%nPE /= 1) &
      call TLS_fatal ('get_cM_compressed: get_cM is not currently implemented to run in parallel')
    if(SS%Method /= DO_SOLVE_SVD_LSLR) &
      call TLS_fatal ('get_cM_compressed: coefficient matrix only available with SVD solve technique')
    doBCs=.false.; if(PRESENT(do_BCs))doBCs=do_BCs
    if(doBCs) &
      call TLS_fatal ('get_cM_compressed: get_cM is not currently implemented to handle boundary conditions')
    ALLOCATE(cM(nfc*ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('get_cM_compressed: memory allocation failure for coefficient structure')
    idx = 1
    CELL_LOOP: do c=1,ncells
      NumFacesCell => SS%NFN_Struct(c)%ptr  ! Access face nghbr count
    ! NOTE: processing of this object as well as subsequent use might
    ! be optimized by use of the fact that a given face is shared between
    ! two cells (for faces not on an exterior boundary).
      FACE_LOOP: do f=1,nfc
        NFN =  NumFacesCell(f)      
      ! cM(idx)%n_phi initialized to 0 => 0 for degenerate case
        if((Mesh(c)%Ngbr_Cell(f) /= DEGENERATE_FACE))then
          if(SS%done(f,c))then
          ! Note: DO_cM has not been filled for "done" components; must copy
            c2 = Mesh(c)%Ngbr_Cell(f)
            f2 = Mesh(c)%Ngbr_Face(f)
            idx2 = (c2-1)*nfc + f2  ! Calculate idx of already calc'd face
            ALLOCATE(cM(idx)%coeff(NFN), &
                     cM(idx)%g_idx(NFN), STAT=istat)
            if (istat /= 0) call TLS_panic ('get_cM_compressed: memory allocation failure for coefficient matrix')
            cM(idx)%n_phi = NFN
            cM(idx)%coeff(:) = cM(idx2)%coeff(:)
            cM(idx)%g_idx(:) = cM(idx2)%g_idx(:)
          else
            cM(idx)%n_phi = NFN
            ALLOCATE(cM(idx)%coeff(NFN), &
                     cM(idx)%g_idx(NFN), STAT=istat)
            if (istat /= 0) call TLS_panic ('get_cM_compressed: memory allocation failure for coefficient matrix')
          ! Find the last index in the contiguous SVD_U for a given cell
            lb = SS%SVDlb(f,c); ub = SS%SVDub(f,c)                
            DO_cM =>SS%SVD_cM(c)%Mat(:,lb:ub)
            cM(idx)%coeff(:) = 0.0_r8  ! initialize coefficient array
            NGHBR_LOOP: do n=1,NFN
            ! Calculate "coefficient vector" A_face dot BinvG (MH - pg 7)
            ! It is assumed that the "diffusion coefficient" has been captured
            ! in the calculation of the discrete_operator LSLR coeffs.
            ! NOTE: DO_cM(1,:) refers to face Phi portion of matrix excluded here;
            !       DO_cM(2:ndim+1,:) refers to face gradient portion of matrix
            ! NOTE2: DO_cM may contain BC components; 1st NFN is face neighbor
            !        only portion
              cM_tmp = 0.0_r8
              do m=1,ndim
                cM_tmp = cM_tmp + Cell(c)%Face_Normal(m,f)*DO_cM(m+1,n)
              end do
              cM(idx)%coeff(n) = Cell(c)%Face_Area(f)*cM_tmp
              cM(idx)%g_idx(n) = SS%NghIdx(f,c)%ptr(n)
            end do NGHBR_LOOP
          endif  ! if(SS%done
        endif  ! if((Mesh(c)%Ngbr_Cell(f)
        idx = idx + 1  ! Increment global matrix index
      end do FACE_LOOP
    end do CELL_LOOP
  END FUNCTION DO_get_cM_compressed


  SUBROUTINE DO_destroy_cM_compressed(cM)
    !=======================================================================
    ! Purpose(s): To destroy the data structure for the compressed form 
    !             of the "coefficient matrix".
    !
    !=======================================================================
    use do_base_types,    only: cMat_row_type
    use parameter_module,  only: nfc,ncells

    type(cMat_row_type), pointer, dimension(:) :: cM
    integer :: idx
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(.not.ASSOCIATED(cM)) &
      call TLS_panic ('DO_destroy_cM_compressed: attempt to DEALLOCATE unassociated coefficient matrix')
    do idx=1,nfc*ncells
      DEALLOCATE(cM(idx)%coeff,cM(idx)%g_idx)
    end do
    DEALLOCATE(cM)
  END SUBROUTINE DO_destroy_cM_compressed


  FUNCTION DO_get_cM_full(SS,do_BCs) RESULT(cM)
    !=======================================================================
    ! Purpose(s): To instantiate the full form of the "coefficient
    !             matrix" described in 'Discussion of the T Diffusion
    !             Discretization', M. Hall, 11/06/02
    !
    !=======================================================================
    use do_base_types, only: DO_Specifier,DO_SOLVE_SVD_LSLR
    use mesh_module, only: Cell,Mesh,DEGENERATE_FACE
    use parallel_info_module, only: p_info
    use parameter_module, only: ndim,nfc,ncells

    type(DO_Specifier), target,   intent(IN) :: SS
    logical,  optional, intent(IN) :: do_BCs

    real(r8), pointer, dimension(:,:) :: cM 
    real(r8) :: cM_tmp

    logical :: doBCs
    integer :: istat, idx, idx2

  ! SVD solution on "Design Matrix"
    real(r8),pointer,dimension(:,:) :: DO_cM => NULL()

    integer :: c,c2,f,f2,j,n,m,NFN,lb,ub
    integer, pointer, dimension(:) :: NumFacesCell => NULL()
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    cM => NULL()
    if(p_info%nPE /= 1) &
      call TLS_fatal ('get_cM_full: get_cM is not currently implemented to run in parallel')
    if(SS%Method /= DO_SOLVE_SVD_LSLR) &
      call TLS_fatal ('get_cM_full: coefficient matrix only available with SVD solve technique')
    doBCs=.false.; if(PRESENT(do_BCs))doBCs=do_BCs
    if(doBCs) &
      call TLS_fatal ('get_cM_full: get_cM is not currently implemented to handle boundary conditions')
    ALLOCATE(cM(ncells,nfc*ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('get_cM_full: memory allocation failure for coefficient structure')
    cM = 0.0_r8
    idx = 1
    CELL_LOOP: do c=1,ncells
      NumFacesCell => SS%NFN_Struct(c)%ptr  ! Access face nghbr count
    ! NOTE: processing of this object as well as subsequent use might
    ! be optimized by use of the fact that a given face is shared between
    ! two cells (for faces not on an exterior boundary).
      FACE_LOOP: do f=1,nfc
        if((Mesh(c)%Ngbr_Cell(f) /= DEGENERATE_FACE))then
          NFN =  NumFacesCell(f)
          if(SS%done(f,c))then
          ! Note: DO_cM has not been filled for "done" components; must copy
            c2 = Mesh(c)%Ngbr_Cell(f)
            f2 = Mesh(c)%Ngbr_Face(f)
            idx2 = (c2-1)*nfc + f2  ! Calculate idx of already calc'd face
            NGHBR_LOOP_copy: do n=1,NFN
            ! Copy "coefficient vector" A_face dot BinvG (MH - pg 7)
            ! It is assumed that the "diffusion coefficient" has been captured
            ! in the calculation of the discrete_operator LSLR coeffs DO_cM
            ! Find Phi vector component associated with matrix element
              j = SS%NghIdx(f2,c2)%ptr(n)
              cM(j,idx) = cM(j,idx2)
            end do NGHBR_LOOP_copy
          else
          ! Find the last index in the contiguous SVD_U for a given cell
            lb = SS%SVDlb(f,c); ub = SS%SVDub(f,c)
            DO_cM =>SS%SVD_cM(c)%Mat(:,lb:ub)
            NGHBR_LOOP: do n=1,NFN
            ! Calculate "coefficient vector" A_face dot BinvG (MH - pg 7)
            ! It is assumed that the "diffusion coefficient" has been captured
            ! in the calculation of the discrete_operator LSLR coeffs DO_cM
            ! NOTE: DO_cM(1,:) refers to face Phi portion of matrix excluded here;
            !       DO_cM(2:ndim+1,:) refers to face gradient portion of matrix
            ! NOTE2: DO_cM may contain BC components; 1st NFN is face neighbor
            !        only portion
              cM_tmp = 0.0_r8
              do m=1,ndim
                cM_tmp = cM_tmp + Cell(c)%Face_Normal(m,f)*DO_cM(m+1,n)
              end do
            ! Find Phi vector component associated with matrix element
              j = SS%NghIdx(f,c)%ptr(n)
              cM(j,idx) = Cell(c)%Face_Area(f)*cM_tmp
            end do NGHBR_LOOP
          endif  ! if(SS%done
        endif  ! if((Mesh(c)%Ngbr_Cell(f)
        idx = idx + 1  ! Increment global matrix index
      end do FACE_LOOP
    end do CELL_LOOP
  END FUNCTION DO_get_cM_full

  SUBROUTINE DO_destroy_cM_full(cM)
    !=======================================================================
    ! Purpose(s): To destroy the data structure for the full form 
    !             of the "coefficient matrix".
    !
    !=======================================================================
    real(r8), pointer, dimension(:,:) :: cM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(.not.ASSOCIATED(cM)) &
      call TLS_panic ('DO_destroy_cM_full: attempt to DEALLOCATE unassociated coefficient matrix')
    DEALLOCATE(cM)
  END SUBROUTINE DO_destroy_cM_full

END MODULE DO_SOLVE_SPECIFIER
