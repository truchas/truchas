!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DO_DISCRETE_OPERATORS
  !=======================================================================
  ! Purpose(s):
  !   Define various procedures to perform discrete operations
  !   such as gradient, divergence, etc.
  !
  !   Public Interfaces:      DO_GRADIENT_FACE
  !                           DO_FACE_SOLVE
  !                           DO_UPDATE_WEIGHTS
  !                           DO_GoodSolution
  !                           DO_GoodPhiSolution
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
  public :: DO_GRADIENT_FACE, DO_FACE_SOLVE, DO_UPDATE_WEIGHTS

  ! Public Functions
  public :: DO_GoodSolution, DO_GoodPhiSolution

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  INTERFACE DO_GRADIENT_FACE
    MODULE PROCEDURE GRADIENT_FACE
  END INTERFACE

  INTERFACE DO_FACE_SOLVE
    MODULE PROCEDURE FACE_SOLVE
  END INTERFACE

  INTERFACE DO_UPDATE_WEIGHTS
    MODULE PROCEDURE UpdateWeights
  END INTERFACE

  INTERFACE DO_GoodSolution
    MODULE PROCEDURE GoodSolution
  END INTERFACE

  INTERFACE DO_GoodPhiSolution
    MODULE PROCEDURE GoodPhiSolution
  END INTERFACE

CONTAINS

  SUBROUTINE GRADIENT_FACE (Phi, SolveSpec, Grad, Phi_Face)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face-centered gradient (Gradient)
    !   and perhaps the face-centered value of
    !   of a cell-centered scalar quantity Phi on face f.
    !
    !=======================================================================
    use cutoffs_module,    only: alittle
    use do_base_types,    only: DO_Specifier, DO_SOLVE_LU_LSLR, DO_SOLVE_SVD_LSLR, DO_SOLVE_ORTHO
    use parameter_module,  only: ndim,nfc,ncells
    use mesh_module,      only: Mesh

    ! Arguments
    real(r8),dimension(ncells),             intent(IN)    :: Phi
    type(DO_Specifier), target,                         intent(INOUT) :: SolveSpec
    real(r8),dimension(ndim,nfc,ncells),    intent(INOUT) :: Grad
    real(r8),dimension(nfc,ncells),OPTIONAL,intent(OUT)   :: Phi_Face

    ! Local Variables
    integer                          :: n, f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(SolveSpec%ID<0) call TLS_panic ('GRADIENT_FACE: called with unassociated SolveSpec')

     SELECT_METHOD: SELECT CASE(SolveSpec%Method)
       CASE(DO_SOLVE_LU_LSLR)
        ! LSLR via LU Method
          call FACE_LSLR_LU(Phi,SS=SolveSpec,GRAD_PHI=Grad,FACE_PHI=Phi_Face)
        ! Eliminate Face Gradient Noise
          do n = 1, ndim
            Grad(n,:,:) = MERGE(0.0_r8,Grad(n,:,:),ABS(Grad(n,:,:)) <= alittle)
          end do
        ! Note that the following step zeros out gradient values on boundary
        ! faces. It is a hold over from the old LSLR implementation and should
        ! likely be removed at some point.
          do f = 1,nfc
            do n = 1,ndim
              where(Mesh%Ngbr_Cell(f) == 0)Grad(n,f,:) = 0.0_r8
            end do
          end do
       CASE(DO_SOLVE_SVD_LSLR)
        ! LSLR via SVD Method
          call FACE_LSLR_SVD(Phi,SS=SolveSpec,GRAD_PHI=Grad,FACE_PHI=Phi_Face)
        ! Eliminate Face Gradient Noise
          do n = 1, ndim
            Grad(n,:,:) = MERGE(0.0_r8,Grad(n,:,:),ABS(Grad(n,:,:)) <= alittle)
          end do
        ! Note that the following step zeros out gradient values on boundary
        ! faces. It is a hold over from the old LSLR implementation and should
        ! likely be removed at some point. It remains here for consistency with
        ! the LU solution.
          do f = 1,nfc
            do n = 1,ndim
              where(Mesh%Ngbr_Cell(f) == 0)Grad(n,f,:) = 0.0_r8
            end do
          end do
       CASE(DO_SOLVE_ORTHO)
        ! Orthogonal stencil
          call FACE_ORTHO(Phi,SS=SolveSpec,GRAD_PHI=Grad,FACE_PHI=Phi_Face)
       CASE DEFAULT
          call TLS_fatal ('GRADIENT_FACE: unknown solution method specified')
     END SELECT SELECT_METHOD
  END SUBROUTINE GRADIENT_FACE


  SUBROUTINE FACE_SOLVE (Phi, SolveSpec, Phi_Face)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face-centered value
    !   of a cell-centered scalar quantity Phi on face f.
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier, DO_SOLVE_LU_LSLR, DO_SOLVE_SVD_LSLR, DO_SOLVE_ORTHO
    use parameter_module,  only: nfc,ncells

    ! Arguments
    real(r8),dimension(ncells),             intent(IN)    :: Phi
    type(DO_Specifier),                                 intent(INOUT) :: SolveSpec
    real(r8),dimension(nfc,ncells),         intent(INOUT) :: Phi_Face

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(SolveSpec%ID<0) call TLS_panic ('FACE_SOLVE: called with unassociated SolveSpec')

     SELECT_METHOD: SELECT CASE(SolveSpec%Method)
       CASE(DO_SOLVE_LU_LSLR)
        ! LSLR via LU Method
          call FACE_LSLR_LU(Phi,SS=SolveSpec,FACE_PHI=Phi_Face)
       CASE(DO_SOLVE_SVD_LSLR)
        ! LSLR via SVD Method
          call FACE_LSLR_SVD(Phi,SS=SolveSpec,FACE_PHI=Phi_Face)
       CASE(DO_SOLVE_ORTHO)
        ! Orthogonal stencil
          call FACE_ORTHO(Phi,SS=SolveSpec,FACE_PHI=Phi_Face)
       CASE DEFAULT
          call TLS_fatal ('GRADIENT_FACE: unknown solution method specified')
     END SELECT SELECT_METHOD
  END SUBROUTINE FACE_SOLVE


  SUBROUTINE FACE_ORTHO(Phi,SS,Grad_Phi,Face_Phi)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face-centered gradient (Gradient)
    !   of a cell-centered scalar quantity Phi on face f.
    !
    !=======================================================================
    use do_base_types,  only: DO_Specifier
    use do_solve_specifier,  only: dX_Scaled
    use gs_module,      only: EE_GATHER
    use parameter_module,  only: ndim,nfc,ncells
    use mesh_module,    only: Mesh
    
    ! Arguments
    real(r8), dimension(ncells), intent(IN) :: Phi
    type(DO_Specifier), target, intent(IN) :: SS
    real(r8),dimension(ndim,nfc,ncells),OPTIONAL,intent(OUT) :: Grad_Phi
    real(r8),dimension(nfc,ncells)     ,OPTIONAL,intent(OUT) :: Face_Phi

    ! Local Variables
    integer                          :: n, f, f2, j, j2
    real(r8),dimension(nfc,ncells)      :: Phi_e
    real(r8)                            :: dPhi
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(.not.ALLOCATED(dX_Scaled)) call TLS_panic ('Face_ORTHO: called with unallocated dX_Scaled')
    
    call EE_GATHER (DEST=Phi_e, SRC=Phi)
    if(PRESENT(Grad_Phi))then
      ! Calculate orthogonal face gradients
      CELL_LOOP: do j=1,ncells
         FACE_LOOP: do f = 1,nfc
            if(Mesh(j)%Ngbr_Cell(f) == 0)then
              Grad_Phi(:,f,j) = 0.0_r8
            else
            ! Calculate value for given face only once;
            ! copy value to neighbor shared face
              if(SS%done(f,j))then
                j2 = Mesh(j)%Ngbr_Cell(f)
                f2 = Mesh(j)%Ngbr_Face(f)
                Grad_Phi(:,f,j) = Grad_Phi(:,f2,j2)
              else
                dPhi = Phi(j) - Phi_e(f,j)
                do n = 1,ndim
                   Grad_Phi(n,f,j) = dPhi * SS%dX_scaled(n,f,j)
                end do
               endif
            endif ! if(Mesh(j)%Ngbr_Cell(f) == 0)
         end do FACE_LOOP
      end do CELL_LOOP
     endif  ! if(PRESENT(Grad_Phi)

     ! Calculate orthogonal face values, if requested, and use the weights
     ! sent to the LSLR calculation; for now, this routine will regard weights
     ! as either zero or one; zero if weight = zero, one for any non-zero value

     if(PRESENT(Face_Phi))then
        Face_Phi = 0.0_r8 ! set a default value for those faces between two void cells
        if(ASSOCIATED(SS%W_Ortho))then
          do j = 1,ncells
             if (SS%W_Ortho(j) == 0.0_r8) CYCLE
             do f = 1,nfc
                if (SS%W_Ortho_Nghbr(f,j) == 0.0_r8) then
                   Face_Phi(f,j) = Phi(j)
                else
                   Face_Phi(f,j) = (Phi(j) + Phi_e(f,j)) / 2
                end if
             end do
          end do
        else ! no Weight specified; assume one everywhere
          do j = 1,ncells
             do f = 1,nfc
                if (Mesh(j)%Ngbr_Cell(f) == 0.0_r8) then
                   Face_Phi(f,j) = Phi(j)
                else
                   Face_Phi(f,j) = (Phi(j) + Phi_e(f,j)) / 2
                end if
             end do
          end do
        endif  ! if(ASSOCIATED(SS%W_Ortho)
     endif  ! if(PRESENT(Face_Phi)
  END SUBROUTINE FACE_ORTHO


  SUBROUTINE FACE_LSLR_LU(Phi,SS,Grad_Phi,Face_Phi)
    !=======================================================================
    ! Purpose(s):
    !   Compute the least-squares linear reconstruction gradient of Phi
    !   at face centers (Gradient_Phi) via LU decomposition with full pivot
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier,dX_Type,SField_Type
    use do_update_module, only: FGetPhiValues
    use do_solve_module,  only: do_lu_solve
    use parameter_module,  only: ndim,nfc,ncells
    use mesh_module,      only: Mesh,DEGENERATE_FACE

    ! Arguments
    real(r8), dimension(ncells), intent(IN)    :: Phi
    type(DO_Specifier), target, intent(INOUT) :: SS
    real(r8),dimension(ndim,nfc,ncells),OPTIONAL,intent(OUT) :: Grad_Phi
    real(r8),dimension(nfc,ncells)     ,OPTIONAL,intent(OUT) :: Face_Phi

    ! Local Variables
    integer                               :: Face,NFC_F
    integer                               :: n,c1,c2,f2
    integer,dimension(:),pointer,save     :: NumFacesCell
    type(dX_Type),  pointer,dimension(:)                 :: dXList
    type(SField_Type), pointer,dimension(:)              :: WList
    type(SField_Type), pointer,dimension(:)              :: PhiVal
    real(r8)                               :: Weight
    real(r8),dimension(4),save             :: dX,RHS
    real(r8)                               :: PhiValCell
  ! LU solution on "Normal Matrix"
    integer,pointer,dimension(:)   :: R1
    integer,pointer,dimension(:)   :: R2
    real(r8),pointer,dimension(:,:)   :: LHS


    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    dX  = 0.0_r8
    CELL_LOOP: do c1=1,ncells
      dXList => SS%dX_Struct(:,c1)
      PhiVal => FGetPhiValues(SS,c1,Phi)  ! Update neighbor and BC Phi values for the cell.
      WList  => SS%TotW_Struct(:,c1)
      NumFacesCell => SS%NF_Struct(c1)%ptr
      FACE_LOOP: do Face=1,nfc
        if(SS%done(Face,c1))cycle FACE_LOOP
        RHS = 0.0_r8
        if (.not. (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE))then
          ! Loop over all face components
          NFC_F =  NumFacesCell(Face)
          NGHBR_and_BC_LOOP: do n = 1, NFC_F
            PhiValCell = PhiVal(Face)%FData(n)
          ! Fill portion of dX according to problem physical dim
            dX(1:ndim+1) = dXList(Face)%FData(:,n)
            Weight = WList(Face)%FData(n)
            ! Add coordinate deltas (dX) for this neighbor to the
            ! LSLR matrix components.
            RHS(1) = RHS(1) + Weight*dX(1)*PhiValCell
            RHS(2) = RHS(2) + Weight*dX(2)*PhiValCell
            RHS(3) = RHS(3) + Weight*dX(3)*PhiValCell
            RHS(4) = RHS(4) + Weight*dX(4)*PhiValCell
          end do NGHBR_and_BC_LOOP
          LHS =>SS%ALU(:,:,Face,c1)
          R1  =>SS%row1(:,Face,c1)
          R2  =>SS%row2(:,Face,c1)
          call do_lu_solve(LHS,RHS,ndim+1,ndim,R1,R2)
        endif  ! if (.not. (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE)

        c2 = Mesh(c1)%Ngbr_Cell(Face)
        f2 = Mesh(c1)%Ngbr_Face(Face)
        if(PRESENT(Face_Phi))then
        ! Actually solved for Phi/d where d is dimensional parameter 
        ! to keep matrix components scale invariant. See DO_INIT_LSLR_SS
          Face_Phi(Face,c1) = dX(1)*RHS(1)
          if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
            Face_Phi(f2,c2) = dX(1)*RHS(1)
          endif
        endif
        if(PRESENT(Grad_Phi))then
          SELECT_DIMENSION: SELECT CASE(ndim)
            CASE (3)
              Grad_Phi(1,Face,c1) = RHS(2)
              Grad_Phi(2,Face,c1) = RHS(3)
              Grad_Phi(3,Face,c1) = RHS(4)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Grad_Phi(1,f2,c2) = RHS(2)
                Grad_Phi(2,f2,c2) = RHS(3)
                Grad_Phi(3,f2,c2) = RHS(4)
              endif
            CASE (2)
              Grad_Phi(1,Face,c1) = RHS(2)
              Grad_Phi(2,Face,c1) = RHS(3)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Grad_Phi(1,f2,c2) = RHS(2)
                Grad_Phi(2,f2,c2) = RHS(3)
              endif
            CASE (1)
              Grad_Phi(1,Face,c1) = RHS(2)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Grad_Phi(1,f2,c2) = RHS(2)
              endif
          END SELECT SELECT_DIMENSION
        endif  ! if(PRESENT(Grad_Phi))
      end do FACE_LOOP
    end do CELL_LOOP
  END SUBROUTINE FACE_LSLR_LU


  SUBROUTINE FACE_LSLR_SVD(Phi,SS,Grad_Phi,Face_Phi)
    !=======================================================================
    ! Purpose(s):
    !   Compute the least-squares linear reconstruction gradient of Phi
    !   at face centers (Gradient_Phi) via Singular Value Decomposition.
    !
    !   NOTE: This routine has an obvious extension. Since each component
    !         of Grad and Face_Phi are complete decoupled, could implement
    !         a routine which simply calculates Face_Phi with no reference
    !         to Grad and vice versa.
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier,dX_Type,SField_Type
    use do_update_module, only: FGetPhiValues
    use parameter_module,  only: ndim,nfc,ncells
    use mesh_module,      only: Mesh,DEGENERATE_FACE

    ! Arguments
    real(r8), dimension(ncells), intent(IN)    :: Phi
    type(DO_Specifier), target,                intent(INOUT) :: SS
    real(r8),dimension(ndim,nfc,ncells),OPTIONAL,intent(OUT) :: Grad_Phi
    real(r8),dimension(nfc,ncells)     ,OPTIONAL,intent(OUT) :: Face_Phi

    ! Local Variables
    integer :: Face,NFC_F,lb,ub
    integer :: n,m,c1,c2,f2
    integer, pointer, save :: NumFacesCell(:)
    type(dX_Type), pointer :: dXList(:)
    type(SField_Type), pointer :: WList(:), PhiVal(:)
    real(r8) :: Weight, RHS, dX1, PhiValCell
    real(r8), save :: X(ndim+1)

  ! SVD solution on "Design Matrix"
    real(r8), pointer :: cM(:,:)


    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    CELL_LOOP: do c1=1,ncells
      dXList => SS%dX_Struct(:,c1)
      PhiVal => FGetPhiValues(SS,c1,Phi) ! Update neighbor and BC Phi values for the cell.
      WList => SS%TotW_Struct(:,c1) ! Access current weight info
      NumFacesCell => SS%NF_Struct(c1)%ptr ! Access face count
    ! Find the last index in the contiguous SVD_U for a given cell
      FACE_LOOP: do Face=1,nfc
        if(SS%done(Face,c1))cycle FACE_LOOP
        X = 0.0_r8; RHS = 0.0_r8
        if (.not. (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE))then
          lb = SS%SVDlb(Face,c1); ub = SS%SVDub(Face,c1)
          cM =>SS%SVD_cM(c1)%Mat(:,lb:ub)
          ! Loop over all face components
          NFC_F =  NumFacesCell(Face)
          NGHBR_and_BC_LOOP: do n = 1, NFC_F
            PhiValCell = PhiVal(Face)%FData(n)
            Weight = WList(Face)%FData(n)
          ! Since SVD uses only dX on LHS rather than dX1*dX2
          ! as in normal matrix formulation, use sqrtW to weight
            RHS = sqrt(Weight)*PhiValCell
            if(ndim==3)then  ! Optimize 3D case
              X(4) = X(4) + cM(4,n)*RHS
              X(3) = X(3) + cM(3,n)*RHS
              X(2) = X(2) + cM(2,n)*RHS
              X(1) = X(1) + cM(1,n)*RHS
            else
              do m=1,ndim+1
                X(m) = X(m) + cM(m,n)*RHS
              end do
            endif
          end do NGHBR_and_BC_LOOP
        endif  ! if (.not. (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE)

        c2 = Mesh(c1)%Ngbr_Cell(Face)
        f2 = Mesh(c1)%Ngbr_Face(Face)
        if(PRESENT(Face_Phi))then
        ! dXList(Face)%FData(:,1) are identical. 
        ! Actually solved for Phi/d where d is dimensional parameter 
        ! to keep matrix components scale invariant. See DO_INIT_LSLR_SS
          dX1= dXList(Face)%FData(1,1)
          Face_Phi(Face,c1) = dX1*X(1)
          if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
            Face_Phi(f2,c2) = dX1*X(1)
          endif
        endif  
        if(PRESENT(Grad_Phi))then
          if(ndim==3)then  ! Optimize 3D case
            Grad_Phi(1,Face,c1) = X(2)
            Grad_Phi(2,Face,c1) = X(3)
            Grad_Phi(3,Face,c1) = X(4)
            if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
              Grad_Phi(1,f2,c2) = X(2)
              Grad_Phi(2,f2,c2) = X(3)
              Grad_Phi(3,f2,c2) = X(4)
            endif
          else
            do m=1,ndim
              Grad_Phi(m,Face,c1) = X(m+1)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Grad_Phi(m,f2,c2) = X(m+1)
              endif
            end do
          endif
        endif  ! if(PRESENT(Grad_Phi))
      end do FACE_LOOP
    end do CELL_LOOP
  END SUBROUTINE FACE_LSLR_SVD


  SUBROUTINE UpdateWeights(SolveSpec,Weights)
    !=======================================================================
    ! Purpose(s):
    !
    ! Subroutine to update weights associated with given cell
    !
    !=======================================================================
    use do_base_types,    only: DO_Specifier, DO_SOLVE_LU_LSLR, DO_SOLVE_SVD_LSLR, DO_SOLVE_ORTHO
    use do_update_module, only: UpdateLSLRWeights
    use gs_module,        only: EE_GATHER
    use parameter_module,  only: ncells
 
    ! Arguments
    type(DO_Specifier), intent(INOUT) :: SolveSpec
    real(r8), dimension(ncells), intent(IN) :: Weights

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! NOTE: This will also catch error where SS (SolveSpec) has somehow become associated
    ! with both weighted ORTHO and LSLR solution techniques... clearly an error
    if(.not.(ASSOCIATED(SolveSpec%W_Ortho)) .and. .not.(ASSOCIATED(SolveSpec%TotW_Struct)) ) &
      call TLS_panic ('UpdateWeights: called with unassociated SolveSpec weights')       

    if( ASSOCIATED(SolveSpec%W_Ortho) .and. ASSOCIATED(SolveSpec%TotW_Struct) ) &
      call TLS_panic ('UpdateWeights: called with ambiguously associated SolveSpec weights')
                                                         
    METHOD: SELECT CASE(SolveSpec%Method)                       
      case(DO_SOLVE_LU_LSLR)
        call UpdateLSLRWeights(SolveSpec,Weights)
      case(DO_SOLVE_SVD_LSLR)
        call UpdateLSLRWeights(SolveSpec,Weights)
      case(DO_SOLVE_ORTHO)
        SolveSpec%W_Ortho = 0.0_r8
        where (Weights > 0.0_r8) SolveSpec%W_Ortho = 1.0_r8
        call EE_GATHER(DEST=SolveSpec%W_Ortho_Nghbr, SRC=SolveSpec%W_Ortho)
      case DEFAULT
        call TLS_fatal ('UpdateWeights: attempt to update weights for unknown solution method.')
    END SELECT METHOD
  END SUBROUTINE UpdateWeights


  FUNCTION GoodPhiSolution(SolveSpec)
    use do_base_types,  only: DO_Specifier
    use parameter_module,  only: nfc,ncells
      type(DO_Specifier), intent(IN) :: SolveSpec
      logical, pointer, dimension(:,:) :: GoodPhiSolution
      logical, allocatable, dimension(:,:), target, save :: GPS
      if(.not. ALLOCATED(GPS))ALLOCATE(GPS(nfc,ncells))
      GPS = SolveSpec%SolveFlag(1,:,:)
      GoodPhiSolution => GPS
  END FUNCTION GoodPhiSolution


  FUNCTION GoodSolution(SolveSpec)
    use do_base_types,  only: DO_Specifier
    use parameter_module,  only: ndim,nfc,ncells
      type(DO_Specifier), intent(IN) :: SolveSpec
      logical, pointer, dimension(:,:,:) :: GoodSolution
      logical, allocatable, dimension(:,:,:), target, save :: GS
      if(.not. ALLOCATED(GS))ALLOCATE(GS(ndim+1,nfc,ncells))
      GS = SolveSpec%SolveFlag
      GoodSolution => GS
  END FUNCTION GoodSolution


END MODULE DO_DISCRETE_OPERATORS
