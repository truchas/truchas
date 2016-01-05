!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DO_UPDATE_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Update utilities for various components of the LSLR discrete
  !   operators solution
  !
  !  "Public Interfaces" (for use within discrete operators (do_*.F90) ONLY)
  !                       FGetPhiValues
  !                       UpdateLSLRWeights
  !                       UpdateFaceLU
  !                       UpdateFaceSVD
  !
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

  ! These interfaces are for other discrete operator (do_*.F90) routines ONLY
  public :: FGetPhiValues, UpdateLSLRWeights, UpdateFaceLU, UpdateFaceSVD

CONTAINS

  SUBROUTINE UpdateLSLRWeights(SS,Weights)
    !=======================================================================
    ! Purpose(s): Subroutine to update weights associated with given cell
    !
    !=======================================================================
    use do_base_types,     only: DO_Specifier
    use gs_module,         only: EE_GATHER
    use parameter_module,  only: ncells,nfc
    use mesh_module,       only: Mesh, DEGENERATE_FACE
    use var_vector_module, only: CREATE,FLATTEN,SIZES,REAL_VAR_VECTOR

    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS
    real(r8), dimension(ncells),intent(IN) :: Weights

    ! Local Variables
    integer :: icell,f,FN,NFN,n,Nidx,istat
    logical, save :: first_time=.true.
    ! This is for all the neighbor weights
    type(REAL_VAR_VECTOR),pointer,dimension(:), save :: NewNbr_Weights
    real(r8), pointer, dimension(:) :: NewCellsNbrWeights
    real(r8), pointer, dimension(:) :: CurCellsNbrWeights

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Create persistent storage for all neighbor weights; Since neighbor set depends
    ! only on geometry, not field type, can create single storage one time.
    ! This will have to be revisited if definition of "neighbor" ever becomes
    ! field or time dependent.
    if(first_time)then
      first_time=.false.
      ALLOCATE(NewNbr_Weights(ncells),STAT=istat)
      if (istat /= 0) call TLS_panic ('UpdateLSLRWeights: memory allocation failure for NewNbr_Weights')
      call CREATE(ARRAY=NewNbr_Weights(:), &
                  SIZES=SIZES(Mesh(1:ncells)%Ngbr_Cells_All))
    endif
    call EE_Gather(DEST=NewNbr_Weights(:), SOURCE=Weights(:))

    SS%UpdateDecomp=.false.
    CELL_LOOP: do icell = 1,ncells
      CurCellsNbrWeights => FLATTEN(SS%CurNbr_Weights(icell))
      NewCellsNbrWeights => FLATTEN(NewNbr_Weights(icell))

    ! Weights have already been gathered; ANY is "local"... don't need PGSLib_Global_ANY
      if(.not.ANY(CurCellsNbrWeights(:)/=NewCellsNbrWeights(:)))cycle  ! Weights unchanged

      FACE_LOOP: do f = 1, nfc
        if(SS%done(f,icell))cycle FACE_LOOP
        if (Mesh(icell)%Ngbr_Cell(f) == DEGENERATE_FACE) cycle FACE_LOOP

    ! FN is the running index for the list of weight values associated
    ! with face neighbor; Nidx is the index into the full set of neighbors
    ! for that face neighbor
        FN = 0
        NFN = SS%NFN_Struct(icell)%ptr(f)
        NEIGHBOR_LOOP: do n=1,NFN
          FN = FN + 1
          Nidx = SS%NghIdx(f,icell)%ptr(FN)
          if(CurCellsNbrWeights(Nidx) /= NewCellsNbrWeights(Nidx))then
            CurCellsNbrWeights(Nidx) = NewCellsNbrWeights(Nidx)
            SS%UpdateDecomp(f,icell)=.true.  ! At least one weight has changed; update decomp
            SS%TotW_Struct(f,icell)%FData(FN) = NewCellsNbrWeights(Nidx)*SS%GeoW_Struct(f,icell)%FData(FN)
          endif
        end do NEIGHBOR_LOOP
      end do FACE_LOOP
    end do CELL_LOOP

  ! Weights have already been gathered; ANY is "local"... don't need PGSLib_Global_ANY
    if(ANY(SS%UpdateDecomp))call UpdateDecomposition(SS)

  END SUBROUTINE UpdateLSLRWeights


  SUBROUTINE UpdateDecomposition(SS)
    !=======================================================================
    ! Purpose(s): Updates decomposition associated with given SolveSpecifier
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier,DO_SOLVE_LU_LSLR,DO_SOLVE_SVD_LSLR
    use parameter_module, only: ncells,nfc

    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS

    ! Local Variables
    integer :: Face,c1

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    METHOD: SELECT CASE(SS%Method)
      case(DO_SOLVE_LU_LSLR)
        do c1 = 1,ncells
          FACE_LOOP_LU: do Face=1,nfc
            if(SS%done(Face,c1))cycle FACE_LOOP_LU
            if(SS%UpdateDecomp(Face,c1))call UpdateFaceLU(SS,Face,c1)
          end do FACE_LOOP_LU
        end do
      case(DO_SOLVE_SVD_LSLR)
        do c1 = 1,ncells
          FACE_LOOP_SVD: do Face=1,nfc
            if(SS%done(Face,c1))cycle FACE_LOOP_SVD
            if(SS%UpdateDecomp(Face,c1))call UpdateFaceSVD(SS,Face,c1)
          end do FACE_LOOP_SVD
        end do
    END SELECT METHOD
  END SUBROUTINE UpdateDecomposition


  SUBROUTINE UpdateFaceLU(SS,Face,c1)
    !=======================================================================
    ! Purpose(s):  Updates LU decomposition
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier,dX_Type,SField_Type
    use do_solve_module,  only: do_lu_decomp,do_lu_solve
    use mesh_module,      only: Mesh, DEGENERATE_FACE
    use parameter_module, only: ndim

    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS
    integer, intent(IN) :: Face
    integer, intent(IN) :: c1

    ! Local Variables
    integer :: i,j,n,d1,d2,c2,f2,PivCnt
    real(r8) :: W, SolutionQuality
    real(r8), dimension(ndim+1,ndim+1) :: LHSi
    integer, pointer :: NumFacesCell(:)
    real(r8), pointer :: dX(:)
    type(dX_Type), pointer :: dXList(:)
    type(SField_Type), pointer :: WList(:)

    integer,  pointer :: R1(:), R2(:)
    real(r8), pointer :: LHS(:,:)
    real(r8), pointer :: StandUncert(:)
    logical,  pointer :: PF(:), SF(:)
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    dXList => SS%dX_Struct(:,c1)
    WList => SS%TotW_Struct(:,c1)
    NumFacesCell => SS%NF_Struct(c1)%ptr
    LHS =>SS%ALU(:,:,Face,c1)
    R1  =>SS%row1(:,Face,c1)
    R2  =>SS%row2(:,Face,c1)
    PF  =>SS%PivFlag(:,Face,c1)
    SF  =>SS%SolveFlag(:,Face,c1)
    StandUncert => SS%StandUncert(:,Face,c1)

    LHS = 0.0_r8
    if (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE)then
      do d1 = 1,ndim+1
        LHS(d1,d1) = 1.0_r8
      end do
    else
    ! Loop over all face components
      NGHBR_and_BC_LOOP1: do n = 1, NumFacesCell(Face)
        dX =>dXList(Face)%FData(:,n)
        W = WList(Face)%FData(n)
          
        ! Add coordinate deltas (dX) for this neighbor to the
        ! LSLR matrix components.
        MATRIX_LOOP1: do d1=1,ndim+1
          LHS(d1,d1) = LHS(d1,d1) + W*dX(d1)*dX(d1)
          MATRIX_LOOP1a: do d2 = d1+1,ndim+1
             LHS(d1,d2) = LHS(d1,d2) + W*dX(d1)*dX(d2)
             LHS(d2,d1) = LHS(d1,d2)
          end do MATRIX_LOOP1a
        end do MATRIX_LOOP1
      end do NGHBR_and_BC_LOOP1
    endif  ! if (Mesh(c1)%Ngbr_Cell(f

  ! Find LU decomposition
    call do_lu_decomp(LHS,ndim+1,R1,R2)

    PF(:) = .false.
    SF(:) = .false.
    PivCnt=0
    do i=1,ndim+1
      if(LHS(i,i)>1d200)cycle
      PivCnt = PivCnt+1
      PF(i) = .true. ! This component contributes to solution
    end do
  ! A single component is by definition singular.
  ! Solution for this face will be set to 0
    if(PivCnt==1)PF(:) = .false.

  ! PivFlag is based on permuted LU decomp. Set unpermuted
  ! SolveFlag for "external consumption" since this is the order
  ! in which the solution is returned.
    do i=1,4
      j = R2(i)
      SF(j) = PF(i)
    end do

    LHSi=0.0_r8; SolutionQuality=0.0_r8; StandUncert=0.0_r8
    do i=1,ndim+1
      LHSi(i,i) = 1.0_r8
    end do
    do i=1,ndim+1  
      if(SF(i))then  ! Only include components used in solution
      ! Calculate the inverse matrix
        call do_lu_solve(LHS,LHSi(:,i),ndim+1,ndim,R1,R2)
      ! Accumulate standard uncertainties for components used in solution
        StandUncert(i) = sqrt(LHSi(i,i))
        SolutionQuality = SolutionQuality + StandUncert(i)
      endif
    end do

  ! Standard uncertainties greater than order unity indicate poor solution quality.
  ! Thus, we (somewhat arbitrarily) reject solution for this face when standard
  ! uncertainty is greater than about 2.0 for each dimension of the solution.
  ! JWD: I think components needing removal need to be "pivoted" out of the solution,
  !      not just masked out
!   if(SolutionQuality > 2.0_r8*((ndim+1)+1))then
!     SF(:) = .false.
!     PF(:) = .false.
!   endif

! Set Pivot flag for neighbor cell sharing face if it exists
! Mesh(c1)%Ngbr_Cell(Face) > 0 implies neighbor on processor.
    if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
      c2 = Mesh(c1)%Ngbr_Cell(Face)
      f2 = Mesh(c1)%Ngbr_Face(Face)
      SS%PivFlag(:,f2,c2) = PF(:)
      SS%SolveFlag(:,f2,c2) = SF(:)
      SS%StandUncert(:,f2,c2) = StandUncert(:)
    endif
  END SUBROUTINE UpdateFaceLU


  SUBROUTINE UpdateFaceSVD(SS,Face,c1)
    !=======================================================================
    ! Purpose(s):  Updates SVD decomposition
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier,dX_Type,SField_Type
    use do_solve_module,  only: do_sv_decomp
    use mesh_module,      only: Mesh, DEGENERATE_FACE
    use parameter_module, only: Nx_tot,ndim

    ! Arguments
    type(DO_Specifier), target, intent(INOUT) :: SS
    integer, intent(IN) :: Face
    integer, intent(IN) :: c1

    ! Local Variables
    integer, pointer :: NumFacesCell(:)
    integer :: lb, ub
    integer, save :: NFC_F=0, NFC_Fmax=0
    integer :: i,j,k,n,d1,c2,f2,istat
    real(r8) :: Weight, sqrtW, cMij
    real(r8),dimension(ndim+1) :: sigmaW, invW
    type(dX_Type), pointer :: dXList(:)
    type(SField_Type), pointer :: WList(:)
    real(r8), pointer :: dX(:), cM(:,:), StandUncert(:)
    real(r8), allocatable, save :: U(:,:), W(:), V(:,:)
    logical, pointer :: SF(:)

    real(r8), parameter :: sol_scale=1.d-6
    real(r8), parameter :: stdU_cutoff=2.0_r8
    real(r8) :: Wcutoff
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    dXList =>SS%dX_Struct(:,c1)
    WList  =>SS%TotW_Struct(:,c1)
    NumFacesCell =>SS%NF_Struct(c1)%ptr

    NFC_F =  NumFacesCell(Face)
    if(NFC_Fmax==0)then
      ALLOCATE(U(NFC_F,ndim+1), &
               W(ndim+1), V(ndim+1,ndim+1),   &
               STAT=istat)
      if (istat /= 0) call TLS_panic ('UpdateFaceSVD: memory allocation failure for U')
      NFC_Fmax = NFC_F
    elseif(NFC_F > NFC_Fmax)then  ! Ensure max dim of U is large as necessary
      DEALLOCATE(U)
      ALLOCATE(U(NFC_F,ndim+1), STAT=istat)
      if (istat /= 0) call TLS_panic ('UpdateFaceSVD: subsequent memory allocation failure for U')
      NFC_Fmax = NFC_F
    endif

    ub = SS%SVDub(Face,c1)
    lb = SS%SVDlb(Face,c1)
    cM =>SS%SVD_cM(c1)%Mat(:,lb:ub)
    SF =>SS%SolveFlag(:,Face,c1)
    SF(:)=.true.
    StandUncert => SS%StandUncert(:,Face,c1)

    if(Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE)then
    ! cM is pointing to 0 size array if DEGENERATE
      return
    else
      U=0.0_r8; W=0.0_r8; V=0.0_r8
    ! Loop over all face components
      NGHBR_and_BC_LOOP1: do n = 1,NFC_F
        dX =>dXList(Face)%FData(:,n)
        Weight = WList(Face)%FData(n)
        sqrtW = sqrt(Weight)
        ! LSLR normal matrix components.
        MATRIX_LOOP1: do d1=1,ndim+1
          U(n,d1) = sqrtW*dX(d1)
        end do MATRIX_LOOP1
      end do NGHBR_and_BC_LOOP1
    endif  ! if (Mesh(c1)%Ngbr_Cell(f

    if(Nx_tot(1) == 1)U(:,2)=0.0_r8
    if(Nx_tot(2) == 1)U(:,3)=0.0_r8
    if(Nx_tot(3) == 1)U(:,4)=0.0_r8
    call do_sv_decomp(U,NFC_Fmax,NFC_F,ndim+1,W,V)

    Wcutoff = maxval(W(:))*sol_scale

    sigmaW = 0.0_r8
    do i=1,ndim+1
      if(W(i)<Wcutoff)then
        W(i)=0.0_r8  ! remove singular components
        SF(i) = .false. ! Set solve flag
      else
        do j=1,ndim+1
          sigmaW(j) = sigmaW(j) + (V(j,i)/W(i))**2
        end do
      endif
    end do
! Standard uncertainties greater than order unity indicate poor solution quality.
! Thus, we (somewhat arbitrarily) reject component for this face when 
! standard uncertainty is greater than stdU_cutoff.
    do i=1,ndim+1
      StandUncert(i) = sqrt(sigmaW(i))
! Standard Uncertainty cutoff TBD
!     if(StandUncert(i)>stdU_cutoff)then
!       W(i) = 0.0_r8
!       SF(i) = .false.
!     endif     

      if(W(i) /= 0.0_r8)then  ! Calculate 1/W for coefficient matrix cM below
        invW(i) = 1/W(i)
      else
        invW(i) = 0.0_r8
      endif
    end do

! Calculate "static" portion of solve 1 time (cM is "static"
! so long as weights are not updated).
    do i=1,ndim+1
      do j=1,NFC_F
        cMij = 0.0_r8
        do k=1,ndim+1
          cMij = cMij + (invW(k)*V(i,k)*U(j,k))
        end do
        cM(i,j) = cMij
      end do
    end do

  ! Set Solve flag for neighbor cell sharing face if it exists
  ! Mesh(c1)%Ngbr_Cell(Face) > 0 implies neighbor on processor.
    if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
      c2 = Mesh(c1)%Ngbr_Cell(Face)
      f2 = Mesh(c1)%Ngbr_Face(Face)
      SS%SolveFlag(:,f2,c2) = SF(:)
      SS%StandUncert(:,f2,c2) = StandUncert(:)
    endif
  END SUBROUTINE UpdateFaceSVD


  FUNCTION FGetPhiValues(SS,icell,Phi) RESULT(PhiLcell)
    !=======================================================================
    ! Purpose(s): Accessor function for neighbor and BC field values
    !             associated with given cell
    !=======================================================================
    use bc_data_types
    use do_base_types,  only: DO_Specifier,SField_Type
    use gs_module,      only: EE_GATHER
    use mesh_module,    only: Mesh,DEGENERATE_FACE
    use parameter_module, only: ncells,nfc
    use var_vector_module, only: REAL_VAR_VECTOR, CREATE, SIZES, FLATTEN

    type(DO_Specifier), target, intent(INOUT) :: SS
    integer,intent(IN) :: icell
    real(r8), dimension(ncells), intent(IN)  :: Phi

    type(SField_Type),pointer,dimension(:)         :: PhiLcell

    ! Local Variables
    type(REAL_VAR_VECTOR), pointer, dimension(:), save :: Centers
    type(SField_Type) :: Scalar_E
    integer :: FN,NFN
    integer :: Nidx
    logical :: lret
    integer :: f,j,n,istat
    type(BC_Chart_ID), POINTER :: ChartID
    integer :: Length
    integer, save :: vlen=0
    real(r8), dimension(:), allocatable,save :: Values
    real(r8), dimension(:,:), pointer :: ValuesMultiDOF
    logical, save :: first_time=.true.
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Have to re-gather field values at beginning of each iteration
    if(icell==1)then
      if(SS%MaxBCLen>0)then  ! BCs exist 
      ! Initialize BC search indices
        call BC_Op_Start_Search(SS%DIR_Op)
        if(SS%MaxBCLen>vlen)then  ! Need new Values length
          if(ALLOCATED(Values))DEALLOCATE(Values)
        ! Allocate room for the values to be gathered
          vlen = SS%MaxBCLen
          ALLOCATE(Values(vlen),STAT=istat)
          if (istat /= 0) call TLS_panic ('FGetPhiValues: memory allocation failure for Values')
        endif
      endif
      if(first_time)then
        first_time = .false.
        ALLOCATE (Centers(ncells))
        call CREATE(ARRAY=Centers(:),                &
                    SIZES=SIZES(Mesh(1:ncells)%Ngbr_Cells_All) )
      endif
      call EE_GATHER (DEST=Centers, SOURCE=Phi)
      CELL_LOOP: do j = 1,ncells
        Scalar_E%FData => FLATTEN(Centers(j))
        FACE_LOOP: do f = 1, nfc
          if(SS%done(f,j))cycle FACE_LOOP
          if (Mesh(j)%Ngbr_Cell(f) == DEGENERATE_FACE) cycle FACE_LOOP

      ! FN is the running index for the list of field values associated
      ! with face neighbor; Nidx is the index into the full set of neighbors
      ! for that face neighbor
          FN = 0
          NFN = SS%NFN_Struct(j)%ptr(f)
          NEIGHBOR_LOOP: do n=1,NFN
            FN = FN + 1
            Nidx =SS%NghIdx(f,j)%ptr(FN)
            SS%PHI_Struct(f,j)%FData(FN) = Scalar_E%FData(Nidx)
          end do NEIGHBOR_LOOP
          if(SS%MaxBCLen==0) cycle FACE_LOOP  ! Bypass remainder of loop if no BCs
          ! Check to see if (f, c) is in the bdy operator
          if(SS%DIR_Mask(f,j))then
            lret = BC_Get_Chart(CHARTID=ChartID,OPERATOR=SS%DIR_Op,CELL=j,FACE=f)
          ! Found a non-empty chart, so grab the data from it.
            Length                = BC_Chart_Length(ChartID)
            ! This is fine for BCs with 1 DOF, but will not work in general.
            ValuesMultiDOF        => BC_Chart_Values(ChartID)
            Values(1:Length)      =  ValuesMultiDOF(1,1:Length)
            BC_DIR_LOOP: do n = 1,Length
              FN = FN + 1
              SS%PHI_Struct(f,j)%FData(FN) = Values(n)
            end do BC_DIR_LOOP
          endif
        end do FACE_LOOP
      end do CELL_LOOP
    endif
    PhiLcell => SS%PHI_Struct(:,icell)
  END FUNCTION FGetPhiValues

END MODULE DO_UPDATE_MODULE
