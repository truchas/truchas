!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DISCRETE_DERIVATIVES
  !=======================================================================
  ! Purpose(s):
  !   Define various procedures to perform discrete operations
  !   such as gradient, divergence, etc.
  !
  !   This version is a first attempt to handle dynamic weights.
  !   See the logical array updateDecomp in FUNCTION FGetWeightList
  !
  !   Public Interface(s): GRADIENT_FACE
  !
  ! Contains:
  !
  ! Author(s): Doug Kothe (dbk@lanl.gov)
  !            Jeff Durachta (durachta@bellatlantic.net)
  !            Robert Ferrell (ferrell@cpca.com)
  !
  !=======================================================================
  use kinds, only: r8
  use bc_data_types
  use cutoffs_module,   only: alittle
  use truchas_logging_services
  use parameter_module, only: ncells, nfc, ndim
  implicit none
  private

  ! Public Subroutines
  public :: GRADIENT_FACE, FACE_LSLR, LSLR_SOLVE_LU, LSLR_SOLVE_pivot
  public :: GoodPhiSolution
  public :: update_LSLR

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  real(r8),allocatable,dimension(:,:,:,:) :: ALU
  logical,allocatable,dimension(:,:)   :: updateDecomp
  logical,allocatable,dimension(:,:,:) :: PivFlag
  logical,allocatable,dimension(:,:,:) :: SolveFlag
  integer,allocatable,dimension(:,:,:) :: row1
  integer,allocatable,dimension(:,:,:) :: row2

  ! Internal data structures

  ! Structure to hold geometric dX
    type dX_Type
      real(r8), pointer, dimension(:,:) :: FData
    end type dX_Type

  ! Structure to hold geometric Weights
    type W_Type
      real(r8), pointer, dimension(:) :: FData
    end type W_Type

  ! Structure to hold field values
    type SField_Type
      real(r8), pointer, dimension(:) :: FData
    end type SField_Type

  ! Structure to hold index to face neighbors
    type N_PTR
       integer, pointer, dimension(:) :: ptr
    end type N_PTR 

    type(dX_Type), allocatable, dimension(:,:), target, save :: dX_Struct
    type(W_Type), allocatable, dimension(:,:), target, save :: W_Struct
    type(W_Type), allocatable, dimension(:,:), target, save :: GeoW_Struct
    type(SField_Type), allocatable, dimension(:,:),target, save :: PHI_Struct
    type(N_PTR), allocatable, dimension(:,:), target, save :: NghPtr

  ! Total number of faces to process
    integer, allocatable, dimension(:,:), target, save :: NumFaces
  ! Number of face neighbors to process
    integer, allocatable, dimension(:,:), save :: NumFaceNghbrs
  ! Number of BCs to process
    integer, allocatable, dimension(:,:), save :: NumFaceBCs 

  ! Max number of BCs to process
    integer :: MaxBCLen
#ifdef DO_DIRICHLET
  ! Mask to determine whether worth doing dirichlet loop 
  ! for given cell & face 
    logical, allocatable, dimension(:,:), save :: DIR_Mask
  ! Dirichlet operator
    type(BC_Operator), POINTER, save :: DIR_Op
#endif

    logical, allocatable, dimension(:,:), save :: done

  INTERFACE Update_LSLR
    MODULE PROCEDURE Update_LSLR_All
    MODULE PROCEDURE Update_LSLR_SelectCell
    MODULE PROCEDURE Update_LSLR_SelectFace
  END INTERFACE

CONTAINS

  SUBROUTINE GRADIENT_FACE (Phi, Grad, Phi_Face, Weight, BC_Spec, use_ortho)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face-centered gradient (Gradient)
    !   of a cell-centered scalar quantity Phi on face f.
    !
    !=======================================================================
    use mesh_module,       only: Mesh
    use gs_module,         only: EE_GATHER

    ! Arguments
    real(r8), dimension(ncells),           intent(IN)    :: Phi
    real(r8), dimension(ncells), optional, intent(IN)    :: Weight
    real(r8), dimension(ndim,nfc,ncells),  intent(INOUT) :: Grad
    real(r8), optional, dimension(nfc,ncells), intent(OUT) :: Phi_Face
    type(BC_Specifier), OPTIONAL, intent(INOUT) :: BC_Spec
    logical, intent(IN) :: use_ortho

    ! Local Variables
    integer :: n, f, j, istat
    real(r8), allocatable, dimension(:,:,:), save :: dX_scaled
    real(r8), dimension(nfc,ncells) :: Phi_e
    real(r8) :: dPhi
    real(r8), allocatable, dimension(:)   :: Weight_Ortho
    real(r8), allocatable, dimension(:,:) :: Weight_Ortho_Ngbr

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (.not. use_ortho) then

        ! LSLR Method
        call FACE_LSLR (Phi,GRAD_PHI=Grad,FACE_PHI=Phi_Face,BC_Spec=BC_Spec,WEIGHT=Weight)

        ! Eliminate Face Gradient Noise
        do n = 1, ndim
           Grad(n,:,:) = MERGE(0.0_r8,Grad(n,:,:),ABS(Grad(n,:,:)) <= alittle)
        end do

        do f = 1,nfc
          do n = 1,ndim
            where(Mesh%Ngbr_Cell(f) == 0)Grad(n,f,:) = 0.0_r8
          end do
        end do

    else

       ! Take short cuts if this is an orthogonal hex mesh...
       ! Gather coordinates and Phi of face neighbors.

       if(.not.ALLOCATED(dX_scaled))then
         ALLOCATE(dX_scaled(ndim,nfc,ncells),STAT=istat)
         if (istat /= 0) call TLS_panic ('GRADIENT_FACE: memory allocation failure for array dX')
         call INIT_ORTHOG_OP(dX_scaled)
       endif

       call EE_GATHER (Phi_e, Phi)

       ! Calculate orthogonal face gradients
       CELL_LOOP: do j=1,ncells
          FACE_LOOP: do f = 1,nfc
             if (Mesh(j)%Ngbr_Cell(f) == 0) then
                Grad(:,f,j) = 0.0_r8
             else
                dPhi = Phi(j) - Phi_e(f,j)
                do n = 1,ndim
                   Grad(n,f,j) = dPhi * dX_scaled(n,f,j)
                end do
             end if
          end do FACE_LOOP
       end do CELL_LOOP

       ! Calculate orthogonal face values, if requested, and use the weights
       ! sent to the LSLR calculation; for now, this routine will regard weights
       ! as either zero or one; zero if weight = zero, one for any non-zero value

       if (PRESENT(Phi_Face)) then

          Phi_Face = 0.0_r8 ! set a default value for those faces between two void cells

          WEIGHTS: if (PRESENT(Weight)) then

             ALLOCATE (Weight_Ortho(ncells),STAT=istat)
             if (istat /= 0) call TLS_panic ('GRADIENT_FACE: memory allocation failure for array Weight_Ortho')
             ALLOCATE (Weight_Ortho_Ngbr(nfc,ncells),STAT=istat)
             if (istat /= 0) call TLS_panic ('GRADIENT_FACE: memory allocation failure for array Weight_Ortho_Ngbr')

             Weight_Ortho = 0.0_r8
             where (Weight > 0.0_r8) Weight_Ortho = 1.0_r8
             call EE_GATHER (Weight_Ortho_Ngbr, Weight_Ortho)
        
             do j = 1,ncells
                if (Weight_Ortho(j) == 0.0_r8) CYCLE
                do f = 1,nfc
                   if (Weight_Ortho_Ngbr(f,j) == 0.0_r8) then
                      Phi_Face(f,j) = Phi(j)
                   else
                      Phi_Face(f,j) = 0.5_r8 * (Phi(j) + Phi_e(f,j))
                   end if
                end do
             end do

             DEALLOCATE (Weight_Ortho)
             DEALLOCATE (Weight_Ortho_Ngbr)

          else ! no Weight specified; assume one everywhere

             do j = 1,ncells
                do f = 1,nfc
                   if (Mesh(j)%Ngbr_Cell(f) == 0) then
                      Phi_Face(f,j) = Phi(j)
                   else
                      Phi_Face(f,j) = 0.5_r8 * (Phi(j) + Phi_e(f,j))
                   end if
                end do
             end do

          end if WEIGHTS

       end if

    end if

  END SUBROUTINE GRADIENT_FACE

  SUBROUTINE INIT_ORTHOG_OP (dX_scaled)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face-centered gradient (Gradient)
    !   of a cell-centered scalar quantity Phi on face f.
    !
    !=======================================================================
    use mesh_module, only: Cell, Mesh
    use gs_module, only: EE_GATHER

    ! Arguments
    real(r8), dimension(:,:,:), intent(OUT) :: dX_scaled

    ! Local Variables
    integer :: n, f, j
    real(r8), dimension(ndim) :: dX
    real(r8) :: dX_tmp
    real(r8), dimension(ndim,nfc,ncells) :: X_e
    real(r8) :: Distance

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do n = 1,ndim
       call EE_GATHER (X_e(n,:,:), Cell%Centroid(n))
    end do
    CELL_LOOP: do j=1,ncells
    FACE_LOOP: do f = 1,nfc
       ! Compute Deltas (Cell to Face Neighbors)
       ! Physical Coordinate Deltas
       dX       = 0.0_r8
       Distance = 0.0_r8
       do n = 1,ndim
          dX(n)  = Cell(j)%Centroid(n) - X_e(n,f,j)
          Distance = Distance + dX(n)**2
       end do
       if (Mesh(j)%Ngbr_Cell(f) == 0)then
          dX_scaled(:,f,j) = 0.0
       else
          do n = 1,ndim
            dX_tmp = dX(n)/(Distance + alittle)
            dX_scaled(n,f,j) = dX_tmp
          end do
       end if
    end do FACE_LOOP
    end do CELL_LOOP
  END SUBROUTINE INIT_ORTHOG_OP


  SUBROUTINE FACE_LSLR (Phi,Grad_Phi,Face_Phi,BC_Spec,Weight)
    !=======================================================================
    ! Purpose(s):
    !   Compute the least-squares linear reconstruction gradient of Phi
    !   at face centers (Gradient_Phi).
    !
    !=======================================================================
    use mesh_module, only: Mesh, DEGENERATE_FACE

    ! Arguments
    real(r8), dimension(ncells), intent(IN) :: Phi
    type(BC_Specifier), OPTIONAL, intent(INOUT) :: BC_Spec
    real(r8), dimension(ndim,nfc,ncells), OPTIONAL, intent(OUT) :: Grad_Phi
    real(r8), dimension(nfc,ncells), OPTIONAL, intent(OUT) :: Face_Phi
    real(r8), dimension(ncells), OPTIONAL, intent(IN) :: Weight

    ! Local Variables
    integer :: ncells_pb, Face, n,c1,c2,f2
    integer, dimension(:), pointer, save :: NumFaces
    type(dX_Type), pointer, dimension(:) :: dXList
    type(W_Type),  pointer, dimension(:) :: WList
    type(SField_Type), pointer, dimension(:) :: PhiVal
    real(r8) :: W
    real(r8), dimension(ndim+1) :: dX
    real(r8) :: PhiValCell

    logical, save :: first_time=.true.

    real(r8), allocatable, dimension(:), save :: RHS

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(.not.PRESENT(Grad_Phi) .and. .not.PRESENT(Face_Phi)) &
        call TLS_fatal ('Face_LSLR must be called with Grad_Phi and/or Face_Phi')

    if (first_time) then
      first_time = .false.
      call INIT_DD_SUPPORT_Face(BC_Spec,Phi,Weight)
      ALLOCATE(RHS(ndim+1))
      ALLOCATE(ALU(ndim+1,ndim+1,nfc,ncells))
      ALLOCATE(PivFlag(ndim+1,nfc,ncells))
      ALLOCATE(SolveFlag(ndim+1,nfc,ncells))
      ALLOCATE(row1(ndim+1,nfc,ncells))
      ALLOCATE(row2(ndim+1,nfc,ncells))
    endif

    ncells_pb = ncells  ! This will need to be changed

    SELECT_DIMENSION: SELECT CASE(ndim)
      CASE (3)
        CELL_Loop3D: do c1=1,ncells_pb
          dXList => FGetdXList(c1)
    ! Update neighbor and BC Phi values for the cell.
    ! This must be done before accessing the weights as a time dependent
    ! weight will cause the LU decomposition for that cell to be recalculated.
    ! This requires the current value of Phi.
          PhiVal => FGetPhiValues(c1,Phi)
          WList => FGetWeightList(c1, Weight)
          NumFaces => FGetListSize(c1)
          FACE_LOOP3D: do Face=1,nfc
            if(done(Face,c1))cycle FACE_LOOP3D
            RHS = 0.0
            if (.not. (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE))then
              ! Loop over all face components
              NGHBR_and_BC_LOOP3D: do n = 1, NumFaces(Face)
                PhiValCell = PhiVal(Face)%FData(n)
                dX = dXList(Face)%FData(:,n)
                W = WList(Face)%FData(n)
                ! Add coordinate deltas (dX) for this neighbor to the
                ! LSLR matrix components.
                RHS(1) = RHS(1) + W*dX(1)*PhiValCell
                RHS(2) = RHS(2) + W*dX(2)*PhiValCell
                RHS(3) = RHS(3) + W*dX(3)*PhiValCell
                RHS(4) = RHS(4) + W*dX(4)*PhiValCell
              end do NGHBR_and_BC_LOOP3D
            endif
            call LSLR_SOLVE_pivot(ALU(1,1,Face,c1),RHS,ndim+1,ndim, &
                                       row1(1,Face,c1),row2(1,Face,c1))
            c2 = Mesh(c1)%Ngbr_Cell(Face)
            f2 = Mesh(c1)%Ngbr_Face(Face)
            if(PRESENT(Face_Phi))then
              Face_Phi(Face,c1) = dX(1)*RHS(1)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Face_Phi(f2,c2) = dX(1)*RHS(1)
              endif
            endif
            if(PRESENT(Grad_Phi))then
              Grad_Phi(1,Face,c1) = RHS(2)
              Grad_Phi(2,Face,c1) = RHS(3)
              Grad_Phi(3,Face,c1) = RHS(4)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Grad_Phi(1,f2,c2) = RHS(2)
                Grad_Phi(2,f2,c2) = RHS(3)
                Grad_Phi(3,f2,c2) = RHS(4)
              endif
            endif
          end do FACE_LOOP3D
        end do CELL_Loop3D

      CASE (2)
        CELL_Loop2D: do c1=1,ncells_pb
          dXList => FGetdXList(c1)
    ! Update neighbor and BC Phi values for the cell.
    ! This must be done before accessing the weights as a time dependent
    ! weight will cause the LU decomposition for that cell to be recalculated.
    ! This requires the current value of Phi.
          PhiVal => FGetPhiValues(c1,Phi)
          WList => FGetWeightList(c1)
          NumFaces => FGetListSize(c1)
          FACE_LOOP2D: do Face=1,nfc
            if(done(Face,c1))cycle FACE_LOOP2D
            RHS = 0.0
            if (.not. (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE))then
              ! Loop over all face components
              NGHBR_and_BC_LOOP2D: do n = 1, NumFaces(Face)
                PhiValCell = PhiVal(Face)%FData(n)
                dX = dXList(Face)%FData(:,n)
                W = WList(Face)%FData(n)
                ! Add coordinate deltas (dX) for this neighbor to the
                ! LSLR matrix components.
                RHS(1) = RHS(1) + W*dX(1)*PhiValCell
                RHS(2) = RHS(2) + W*dX(2)*PhiValCell
                RHS(3) = RHS(3) + W*dX(3)*PhiValCell
              end do NGHBR_and_BC_LOOP2D
            endif
            call LSLR_SOLVE_pivot(ALU(1,1,Face,c1),RHS,ndim+1,ndim, &
                                       row1(1,Face,c1),row2(1,Face,c1))

            c2 = Mesh(c1)%Ngbr_Cell(Face)
            f2 = Mesh(c1)%Ngbr_Face(Face)
            if(PRESENT(Face_Phi))then
              Face_Phi(Face,c1) = dX(1)*RHS(1)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Face_Phi(f2,c2) = dX(1)*RHS(1)
              endif          
            endif 
            if(PRESENT(Grad_Phi))then
              Grad_Phi(1,Face,c1) = RHS(2)
              Grad_Phi(2,Face,c1) = RHS(3)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Grad_Phi(1,f2,c2) = RHS(2)
                Grad_Phi(2,f2,c2) = RHS(3)
              endif          
            endif            
          end do FACE_LOOP2D
        end do CELL_Loop2D

      CASE (1)
        CELL_Loop1D: do c1=1,ncells_pb
          dXList => FGetdXList(c1)
    ! Update neighbor and BC Phi values for the cell.
    ! This must be done before accessing the weights as a time dependent
    ! weight will cause the LU decomposition for that cell to be recalculated.
    ! This requires the current value of Phi.
          PhiVal => FGetPhiValues(c1,Phi)
          WList => FGetWeightList(c1)
          NumFaces => FGetListSize(c1)
          FACE_LOOP1D: do Face=1,nfc
            if(done(Face,c1))cycle FACE_LOOP1D
            RHS = 0.0
            if (.not. (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE))then
              ! Loop over all face components
              NGHBR_and_BC_LOOP1D: do n = 1, NumFaces(Face)
                PhiValCell = PhiVal(Face)%FData(n)
                dX = dXList(Face)%FData(:,n)
                W = WList(Face)%FData(n)
                ! Add coordinate deltas (dX) for this neighbor to the
                ! LSLR matrix components.
                RHS(1) = RHS(1) + W*dX(1)*PhiValCell
                RHS(2) = RHS(2) + W*dX(2)*PhiValCell
              end do NGHBR_and_BC_LOOP1D
            endif
            call LSLR_SOLVE_pivot(ALU(1,1,Face,c1),RHS,ndim+1,ndim, &
                                       row1(1,Face,c1),row2(1,Face,c1))
            c2 = Mesh(c1)%Ngbr_Cell(Face)
            f2 = Mesh(c1)%Ngbr_Face(Face)
            if(PRESENT(Face_Phi))then
              Face_Phi(Face,c1) = dX(1)*RHS(1)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Face_Phi(f2,c2) = dX(1)*RHS(1)
              endif          
            endif 
            if(PRESENT(Grad_Phi))then
              Grad_Phi(1,Face,c1) = RHS(2)
              if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
                Grad_Phi(1,f2,c2) = RHS(2)
              endif          
            endif            
          end do FACE_LOOP1D
        end do CELL_Loop1D
    END SELECT SELECT_DIMENSION
  END SUBROUTINE FACE_LSLR

! The "F" in the FGet... functions below stands for "Face"
! Functions for cells can be generalized from these

  ! Accessor function for dXs associated with given cell
  FUNCTION FGetdXList(icell) RESULT(dXLcell)
  integer, intent(IN) :: icell
  type(dX_Type), pointer, dimension(:) :: dXLcell
  dXLcell => dX_Struct(:,icell)
  END FUNCTION FGetdXList

  FUNCTION FGetWeightList(icell, Weight) RESULT(WLcell)
  ! Accessor function for weights associated with given cell

  use gs_module,      only: EE_GATHER
  use mesh_module,    only: Mesh, DEGENERATE_FACE
  use var_vector_module, only: REAL_VAR_VECTOR, CREATE, SIZES, FLATTEN

  integer, intent(IN) :: icell
  type(W_Type), pointer, dimension(:) :: WLcell
  real(r8), dimension(ncells), optional, intent(IN) :: Weight

  ! Local Variables
  integer :: f,FN,n,Nptr
  logical, save :: first_time=.true.
  logical, save :: first_update=.true.
  logical, save :: update_weights=.true.
  ! This is for all the neighbor weights
  type(REAL_VAR_VECTOR), pointer, SAVE, dimension(:) :: Nbr_Weight
  real(r8), pointer, dimension(:) :: CellsNbrWeights

  if(first_time)then
    first_time=.false.
  ! Change updateDecomp(face,cell) to .true. to get weights and 
  ! the LU decomposition updated for a given iteration.
    ALLOCATE(updateDecomp(nfc,ncells))
    updateDecomp=.false.
    call LU_Init()
  elseif(any(updateDecomp(:,icell)))then
    if(.not.PRESENT(Weight)) call TLS_fatal ('FGetWeightList: weights must be present for update')   
    if(first_update)then
      first_update = .false.
      ALLOCATE(Nbr_Weight(ncells))
      ! Create persistent storage for all neighbor weights
      call CREATE(ARRAY = Nbr_Weight(:), &               
                  SIZES = SIZES(Mesh(1:ncells)%Ngbr_Cells_All))
    endif
    if(update_weights)then
      update_weights = .false.
      call EE_Gather(DEST=Nbr_Weight(:), &
                     SOURCE=Weight(:), &
                     RANGE=(/1,ncells/))
    endif
    FACE_LOOP: do f = 1, nfc
      updateDecomp(f,icell)=.false.
      if(done(f,icell))cycle FACE_LOOP
      if (Mesh(icell)%Ngbr_Cell(f) == DEGENERATE_FACE) cycle FACE_LOOP

      CellsNbrWeights => FLATTEN(Nbr_Weight(icell))

    ! FN is the running index for the list of weight values associated
    ! with face neighbor; Nptr is the index into the full set of neighbors
    ! for that face neighbor
      FN = 0
      NEIGHBOR_LOOP: do n=1,NumFaceNghbrs(f,icell)
        FN = FN + 1
        Nptr = NghPtr(f,icell)%ptr(FN)
        W_Struct(f,icell)%FData(FN) = CellsNbrWeights(Nptr)*GeoW_Struct(f,icell)%FData(FN)
      end do NEIGHBOR_LOOP
      call UpdateFaceLU(f,icell)
    end do FACE_LOOP
  endif

  if(icell==ncells)update_weights=.true. !Set flag to gather weights for next iteration
  WLcell => W_Struct(:,icell)
  END FUNCTION FGetWeightList


  SUBROUTINE LU_Init()
  ! Initializes LU decomposition

  ! Local Variables
  integer :: Face,  c1

  do c1 = 1,ncells
    FACE_LOOP: do Face=1,nfc
      if(done(Face,c1))cycle FACE_LOOP
      call UpdateFaceLU(Face,c1)
    end do FACE_LOOP
  end do
  END SUBROUTINE LU_Init


  SUBROUTINE UpdateFaceLU(Face,c1)
  ! Updates LU decomposition

  use mesh_module, only: Mesh, DEGENERATE_FACE

  integer, intent(IN) :: Face
  integer, intent(IN) :: c1

  ! Local Variables
  integer :: istat,i,j,n,d1,d2,c2,f2
  integer, dimension(:), pointer :: NumFacesCell
  integer :: PivCnt
  type(dX_Type), pointer, dimension(:) :: dXList
  type(W_Type),  pointer, dimension(:) :: WList
  type(SField_Type), pointer, dimension(:) :: PhiVal
  real(r8) :: W
  real(r8), dimension(ndim+1) :: dX
  real(r8) :: PhiValCell
  real(r8), allocatable, dimension(:,:), save :: LHS
  real(r8), allocatable, dimension(:,:), save :: LHSi

  if(.not. ALLOCATED(LHS))then
    ALLOCATE(LHS(ndim+1,ndim+1),STAT=istat)
    if (istat /= 0) call TLS_panic ('UpdateFaceLU: memory allocation failure for array LHS')
    ALLOCATE(LHSi(ndim+1,ndim+1),STAT=istat) 
  endif
  dXList => dX_Struct(:,c1)
  PhiVal => PHI_Struct(:,c1)
  WList => W_Struct(:,c1)
  NumFacesCell => NumFaces(:,c1)

  LHS = 0.0
  if (Mesh(c1)%Ngbr_Cell(Face) == DEGENERATE_FACE)then
    do d1 = 1,ndim+1
      LHS(d1,d1) = 1.0_r8
    end do
  else
  ! Loop over all face components
  NGHBR_and_BC_LOOP1: do n = 1, NumFacesCell(Face)
    PhiValCell = PhiVal(Face)%FData(n)
    dX = dXList(Face)%FData(:,n)
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
  call LSLR_SOLVE_LU(LHS,ndim+1,row1(1,Face,c1),row2(1,Face,c1))

  ALU(:,:,Face,c1)=LHS
  PivFlag(:,Face,c1)=.false.
  SolveFlag(:,Face,c1)=.false.
  PivCnt=0
  do i=1,ndim+1
    if(LHS(i,i)>1d200)cycle
    PivCnt = PivCnt+1
    PivFlag(i,Face,c1)=.true. ! This component contributes to solution
  end do
! A single component is by definition singular.
! Solution for this face will be set to 0
  if(PivCnt==1)PivFlag(:,Face,c1)=.false.

! PivFlag is based on permuted LU decomp. Set unpermuted
! SolveFlag for "external consumption" since this is the order
! in which the solution is returned.
  do i=1,4
    j=row2(i,Face,c1)
    SolveFlag(j,Face,c1)=PivFlag(i,Face,c1)
  end do

! Set Pivot flag for neighbor cell sharing face if it exists
! Mesh(c1)%Ngbr_Cell(Face) > 0 implies neighbor on processor.
  if(Mesh(c1)%Ngbr_Cell(Face) > 0) then
    c2 = Mesh(c1)%Ngbr_Cell(Face)
    f2 = Mesh(c1)%Ngbr_Face(Face)
    PivFlag(:,f2,c2) = PivFlag(:,Face,c1)
    SolveFlag(:,f2,c2) = SolveFlag(:,Face,c1)
  endif

  LHSi = 0.d0
  do i=1,ndim+1
    LHSi(i,i) = 1.d0
  end do
  do i=1,ndim+1  
    call LSLR_SOLVE_pivot(LHS,LHSi(:,i),ndim+1,ndim,row1(1,Face,c1),row2(1,Face,c1))
  end do         
! write(39,*) 'SU(Face,Cell)=(',Face,',',c1,'): ',(sqrt(LHSi(i,i)),i=1,ndim+1)
  END SUBROUTINE UpdateFaceLU


  FUNCTION FGetPhiValues(icell,Phi) RESULT(PhiLcell)
  ! Accessor function for neighbor and BC field values 
  ! associated with given cell
  integer, intent(IN) :: icell
  real(r8), dimension(ncells), intent(IN) :: Phi

  type(SField_Type),pointer,dimension(:) :: PhiLcell

  ! Have to re-gather field values at beginning of each iteration
  if(icell == 1)call FGetScalarField(Phi,PHI_Struct)

  PhiLcell => PHI_Struct(:,icell)
  END FUNCTION FGetPhiValues


  SUBROUTINE FGetScalarField(SField,S_Struct)
  ! Accessor function for neighbor and BC scalar field values
  ! associated with given cell
  use gs_module, only: EE_GATHER
  use mesh_module, only: Mesh, DEGENERATE_FACE
  use var_vector_module, only: REAL_VAR_VECTOR, CREATE, SIZES, FLATTEN, DESTROY

  real(r8), dimension(:), intent(IN) :: SField
  type(SField_Type), dimension(:,:), intent(INOUT) :: S_Struct

  type(REAL_VAR_VECTOR), pointer, dimension(:) :: Centers
!  This seems to be unnecessary, and was causing a memory leak as noted below
!  type(BOUNDARY) :: Centers_Bndy

  ! This V_HACK type is needed to conform with the REAL_VAR_VECTOR type.
  type V_HACK
     real(r8), pointer, dimension(:) :: Data
  end type V_HACK
  type(V_HACK) :: Scalar_E
  integer :: FN,Nptr,f,j,n
#ifdef DO_DIRICHLET
  integer :: Length
  logical :: lret
  type(BC_Chart_ID), POINTER :: ChartID
  real(r8), dimension(:), pointer :: Values
  real(r8), dimension(:,:), pointer :: ValuesMultiDOF
#endif

#ifdef DO_DIRICHLET
! Have to re-gather field values at beginning of each iteration
! Initialize BC search indices
  call BC_Op_Start_Search(DIR_Op)
! Allocate room for the values to be gathered
  ALLOCATE(Values(MaxBCLen))
#endif
! This was causing a memory leak 
!  NULLIFY(Centers_Bndy%Data)
  ALLOCATE (Centers(ncells))
  call CREATE(ARRAY = Centers(:),                &
              SIZES = SIZES(Mesh(1:ncells)%Ngbr_Cells_All) )
  call EE_GATHER (DEST     = Centers,         &
                  SOURCE   = SField)
! This was causing a memory leak 
!,      &
!                  BOUNDARY = Centers_Bndy%Data, &
!                  RANGE    = (/1,ncells/))
  CELL_LOOP: do j = 1,ncells
    Scalar_E%Data => FLATTEN(Centers(j))
    FACE_LOOP: do f = 1, nfc

    if(done(f,j))cycle FACE_LOOP
    if (Mesh(j)%Ngbr_Cell(f) == DEGENERATE_FACE) cycle FACE_LOOP

    ! FN is the running index for the list of field values associated
    ! with face neighbor; Nptr is the index into the full set of neighbors
    ! for that face neighbor

    FN = 0
      NEIGHBOR_LOOP: do n=1,NumFaceNghbrs(f,j)
        FN = FN + 1
        Nptr = NghPtr(f,j)%ptr(FN)
        S_Struct(f,j)%FData(FN) = Scalar_E%data(Nptr)
      end do NEIGHBOR_LOOP
#ifdef DO_DIRICHLET
        ! Check to see if (f, c) is in the bdy operator
      if(DIR_Mask(f,j))then
        lret = BC_Get_Chart(CHARTID=ChartID,OPERATOR=DIR_Op,CELL=j,FACE=f)
        ! Found a non-empty chart, so grab the data from it.
        Length                = BC_Chart_Length(ChartID)
        ! Changed this to conform with multi-DOF BCs.
        ! This is fine for DIRICHLET with 1 DOF, but will not work in general.
        ValuesMultiDOF        => BC_Chart_Values(ChartID)
        Values(1:Length)      = ValuesMultiDOF(1,1:Length)
        BC_DIR_LOOP: do n = 1,Length
          FN = FN + 1
          S_Struct(f,j)%FData(FN) = Values(n)
        end do BC_DIR_LOOP
      endif
#endif
    end do FACE_LOOP
  end do CELL_LOOP
#ifdef DO_DIRICHLET
  DEALLOCATE(Values)
#endif
  call DESTROY (ARRAY = Centers)
  DEALLOCATE (Centers)
  END SUBROUTINE FGetScalarField


  FUNCTION FGetListSize(icell) RESULT(NumFacesCell)
  ! Accessor function for the number of neighbor and BC faces to processed
  ! with a given cell
  integer, intent(IN) :: icell
  integer, dimension(:), pointer :: NumFacesCell
  NumFacesCell => NumFaces(:,icell)
  END FUNCTION FGetListSize

  SUBROUTINE INIT_DD_SUPPORT_Face(BC_Spec,Phi,Weight)
    !=======================================================================
    ! Purpose(s): To instantiate and initialize the persistent data 
    !             structures for face neighbor and BC computations
    !
    !             This routine is "simply" a "game" of counting the number
    !             of valid face neighbors and BCs and packing the 
    !             persistent data (dXs, geometric weights, and indices into
    !             the list of valid face neighbor points) into
    !             arrays for the FGet... accessor functions.
    !
    !=======================================================================
    use gs_module, only: EE_GATHER
    use mesh_module, only: Cell,Mesh,BOUNDARY,Is_Face_Ngbr
    use var_vector_module, only: REAL_VAR_VECTOR, CREATE, SIZES, FLATTEN, DESTROY

    ! Arguments
    type(BC_Specifier), OPTIONAL, intent(INOUT) :: BC_Spec
    real(r8), dimension(ncells), intent(IN) :: Phi
    real(r8), dimension(ncells), OPTIONAL, intent(IN) :: Weight

    ! Local Variables
  ! Various index variables
    integer :: d,j,n,f,j2,f2
  ! Running count of the number of valid neighbor faces
    integer :: NumNghbrs
  ! The ...tmp arrays are temporary structures to allocated to hold the
  ! maximum possible extent of the required data until the true extent is
  ! determined by Is_Face_Ngbr, DEGENERATE_FACE and the BC routines
    integer,  allocatable, dimension(:,:)   :: PTRtmp
    real(r8), allocatable, dimension(:,:)   :: CellWtmp
    real(r8), allocatable, dimension(:,:)   :: GeoWtmp
    real(r8), allocatable, dimension(:,:,:) :: dXtmp
  ! The persistent storage for valid indices, field values (which are updated
  ! each interation), geometric weights and dXs
    integer,  pointer, dimension(:),   save :: PTRs
    real(r8), pointer, dimension(:),   save :: PHIs
    real(r8), pointer, dimension(:),   save :: GeoW
    real(r8), pointer, dimension(:),   save :: W
    real(r8), pointer, dimension(:,:), save :: dX
    type(REAL_VAR_VECTOR), pointer, dimension(:,:),save :: X_Centers

    type(BOUNDARY), dimension(ndim) :: X_Centers_Bndy

    ! This V_HACK type is needed to conform with the REAL_VAR_VECTOR type.
    type V_HACK
       real(r8), pointer, dimension(:) :: Data
    end type V_HACK
    type(V_HACK), dimension(ndim) :: Xc
    integer, pointer, dimension(:) :: Ngbr_Cells_Face
    integer :: FN,NFN
    integer :: istat
    logical, save :: first_time=.true.

  ! This is for all the neighbor weights
    type(REAL_VAR_VECTOR), pointer, SAVE, dimension(:) :: Nbr_Weight
    real(r8), dimension(:), pointer :: CellsNbrWeights

    real(r8) :: dX1
    
#ifdef DO_DIRICHLET
    integer :: Length
    type(BC_Chart_ID), POINTER :: ChartID
    real(r8), dimension(:,:), pointer :: Positions
    
#endif
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if(first_time)then
      first_time = .false.
    else
      call TLS_fatal ('INIT_DD_SUPPORT: cannot be called more than once') 
    endif

    ALLOCATE(done(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array DONE')
    done = .false.

  ! All of the allocating in here really needs to be error checked
    ALLOCATE(dX_Struct(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array dX_Struct')
    ALLOCATE(PHI_Struct(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array PHI_Struct')
    ALLOCATE(W_Struct(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array W_Struct')
    ALLOCATE(GeoW_Struct(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array GeoW_Struct')
    ALLOCATE(NghPtr(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array NghPtr')
    ALLOCATE(NumFaces(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array NumFaces')
    ALLOCATE(NumFaceNghbrs(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array NumFaceNghbrs')
    ALLOCATE(NumFaceBCs(nfc,ncells),STAT=istat)
    if (istat /= 0) call TLS_panic ('INIT_DD_SUPPORT: memory allocation failure for array NumFaceBCs')

    ! Initialize boundary arrays.
    do d = 1,ndim
       NULLIFY(X_Centers_Bndy(d)%Data)
    end do
    ! Allocate a vector to store neighboring coordinate data.
    ALLOCATE (X_Centers(ndim,ncells))
    GATHER_LOOP: do d = 1,ndim
      call CREATE (ARRAY = X_Centers(d,:),                &
                   SIZES = SIZES(Mesh(1:ncells)%Ngbr_Cells_All) )
      call EE_GATHER (DEST     = X_Centers(d,:),         &
                      SOURCE   = Cell%Centroid(d),      &
                      BOUNDARY = X_Centers_Bndy(d)%Data, &
                      RANGE    = (/1,ncells/))
    end do GATHER_LOOP

    if(PRESENT(Weight))then
      ALLOCATE(Nbr_Weight(ncells))
      ! Gather the weights from all neighbors, for use by all cells
      call CREATE(ARRAY = Nbr_Weight(:), &
                  SIZES = SIZES(Mesh(1:ncells)%Ngbr_Cells_All))
      call EE_Gather (DEST = Nbr_Weight, SOURCE = Weight)
    endif


    ! Get the DIRICHLET operators
#ifdef DO_DIRICHLET
    DIR_Op => BC_Spec_Get_Operator(BC_Spec, BC_DIRICHLET_Op)
    call BC_Op_Start_Search(DIR_Op)
  ! Note that BC_Get_Length returns nothing (NULL) if there are no values
  ! associated with that operator. By design, the intrinsic MAXVAL must return
  ! something... which happens to be the negative number of greatest magnitude
  ! consistent with KIND of the calling LHS. Hence the use of the array
  ! constructor to make 0 the smallest value which may be returned.
    MaxBCLen = MAXVAL( (/0,BC_Get_Length(BC_Get_Atlas(DIR_Op))/) )
#else
    MaxBCLen = 0
#endif

! The calls to the BC_Op functions are relativly expensive. More important
! the results cannot be "prefetched". In principle, the compiler can arrange
! to "look ahead" in the mask arrays.
# ifdef DO_DIRICHLET
    ALLOCATE(DIR_Mask(nfc,ncells))
    DIR_Mask = .false.
# endif

    NumFaces = 0
    NumFaceNghbrs = 0
    NumFaceBCs = 0
    CELL_LOOP: do j = 1,ncells
      do d = 1, ndim
        Xc(d)%Data => FLATTEN(X_Centers(d,j))
      end do
      Ngbr_Cells_Face => FLATTEN(Mesh(j)%Ngbr_Cells_Face)
      if(PRESENT(Weight))CellsNbrWeights => FLATTEN(Nbr_Weight(j))
      NumNghbrs = SIZES(X_Centers(1,j))
      ALLOCATE(dXtmp(ndim+1,NumNghbrs+MaxBCLen,nfc))
      if(PRESENT(Weight))ALLOCATE(CellWtmp(NumNghbrs+MaxBCLen,nfc))
      ALLOCATE(GeoWtmp(NumNghbrs+MaxBCLen,nfc))
      ALLOCATE(PTRtmp(NumNghbrs+MaxBCLen,nfc))
      GeoWtmp = 0.0
      NEIGHBOR_LOOP: do n=1,NumNghbrs
        FACE_LOOP: do f = 1, nfc
         ! If cell is not a boundary cell or a neighbor cell, cycle
           if(.not.( Mesh(j)%Ngbr_Cell(f)==0 .or. Is_Face_Ngbr(Ngbr_Cells_Face(n),f) ))cycle FACE_LOOP
           NumFaces(f,j) = NumFaces(f,j) + 1
           NumFaceNghbrs(f,j) = NumFaceNghbrs(f,j) + 1
           FN = NumFaces(f,j)
           PTRtmp(FN,f) = n
           ! Compute the distance vector and weight for this neighbor.
           do d = 1, ndim
              dXtmp(d+1,FN,f) = Xc(d)%Data(n) - Cell(j)%Face_Centroid(d,f)
              GeoWtmp(FN,f)   = GeoWtmp(FN,f) + dXtmp(d+1,FN,f)*dXtmp(d+1,FN,f)
           end do
           dXtmp(1,FN,f) = sqrt(GeoWtmp(FN,f))  ! Use distance to scale phi solution
           ! Geometric weight is 1/d**2, where d is the distance 
           ! to the neighbor.
           if (GeoWtmp(FN,f) >= alittle) GeoWtmp(FN,f) = 1.0/GeoWtmp(FN,f)
           if(PRESENT(Weight))CellWtmp(FN,f) = CellsNbrWeights(n)
        end do FACE_LOOP
      end do NEIGHBOR_LOOP

       ! Now add the BC points
# ifdef DO_DIRICHLET
      DIR_FACE_LOOP: do f = 1, nfc
          ! Check to see if (f, c) is in the bdy operator
        if(.not.(BC_Get_Chart  &
                   (CHARTID=ChartID,OPERATOR=DIR_Op,CELL=j,FACE=f)))cycle
        DIR_Mask(f,j) = .true.
        ! Found a non-empty chart, so grab the data from it.
        Length                = BC_Chart_Length(ChartID)
        ALLOCATE(Positions(DIMENSIONALITY(DIR_Op),Length))
        Positions(:,1:Length) = BC_Chart_Positions(ChartID)
        BC_DIR_LOOP: do n = 1,Length
          NumFaces(f,j) = NumFaces(f,j) + 1
          NumFaceBCs(f,j) = NumFaceBCs(f,j) + 1
          FN = NumFaces(f,j)
          PTRtmp(FN,f) = j
          ! Compute the distance vector and weight for this neighbor.
          do d = 1, ndim
             dXtmp(d+1,FN,f) = Positions(d,n) - Cell(j)%Face_Centroid(d,f)
             GeoWtmp(FN,f)   = GeoWtmp(FN,f) + dXtmp(d+1,FN,f)*dXtmp(d+1,FN,f)
          end do
          dXtmp(1,FN,f) = sqrt(GeoWtmp(FN,f))  ! Use distance to scale phi solution
          ! Geometric weight is 1/d**2, where d is the distance
          ! to the neighbor.
          if (GeoWtmp(FN,f) >= alittle) GeoWtmp(FN,f) = 1.0/GeoWtmp(FN,f)
        end do BC_DIR_LOOP
        DEALLOCATE(Positions)
      end do DIR_FACE_LOOP
# endif ! DO_DIRICHLET
  ! Now that we know how many valid items there really are, allocate the
  ! persistent storage and point to it.
  ! Again, this allocation really needs to be error checked
      FACE_LOOP2: do f = 1, nfc
        ALLOCATE( dX(ndim+1,NumFaces(f,j)) )
        ALLOCATE(    W(NumFaces(f,j)) )
        ALLOCATE( GeoW(NumFaces(f,j)) )
        ALLOCATE( PHIs(NumFaces(f,j)) )
        ALLOCATE( PTRs(NumFaces(f,j)) )
        dX1 = sum( dXtmp(1,1:NumFaces(f,j),f) )
        dX(1,1:NumFaces(f,j)) = dX1 / float(NumFaces(f,j))
        dX(2:ndim+1,:) = dXtmp(2:ndim+1,1:NumFaces(f,j),f)
      GeoW   =  GeoWtmp(  1:NumFaces(f,j),f)
         W   =  GeoWtmp(  1:NumFaces(f,j),f)
         NFN = NumFaceNghbrs(f,j)
  ! If non-geometric weights are present, convolve with geometric weight to produce
  ! complete weight value.
        if(PRESENT(Weight))W(1:NFN) = W(1:NFN)*CellWtmp(1:NFN,f)
        PTRs =   PTRtmp(  1:NumFaces(f,j),f)
         dX_Struct(f,j)%FData => dX
       GeoW_Struct(f,j)%FData => GeoW
          W_Struct(f,j)%FData => W
        Phi_Struct(f,j)%FData => PHIs
            NghPtr(f,j)%ptr   => PTRs
        if(Mesh(j)%Ngbr_Cell(f) > 0 .and. .not. done(f,j)) then
          j2 = Mesh(j)%Ngbr_Cell(f)
          f2 = Mesh(j)%Ngbr_Face(f)
          done(f2,j2) = .true.
        endif
      end do FACE_LOOP2
      DEALLOCATE(dXtmp)
      if(PRESENT(Weight))DEALLOCATE(CellWtmp)
      DEALLOCATE(GeoWtmp)
      DEALLOCATE(PTRtmp)
    end do CELL_LOOP
  ! Destroy the X_Centers working array.
    do d = 1,ndim
       call DESTROY (ARRAY = X_Centers(d,:))
    end do
  ! Deallocate a vector to store neighboring coordinate data.
    DEALLOCATE (X_Centers)
    do d = 1,ndim
       DEALLOCATE(X_Centers_Bndy(d)%Data)
    end do
    if(PRESENT(Weight))then
      call DESTROY(ARRAY = Nbr_Weight(:))
      DEALLOCATE(Nbr_Weight)
    endif
  END SUBROUTINE INIT_DD_SUPPORT_Face


  SUBROUTINE LSLR_SOLVE_LU(A,n,row,col)
    !=======================================================================
    ! Purpose:
    !
    !   Compute LU decomposition of A, with row and column pivoting.
    !   Row and column permutation matricies are stored in arrays row and
    !   row2.  For example, the solution to :
    !
    !      A x = b
    !
    !   is also given by the solution to:
    !
    !      R A C C^-1 x = R b
    !
    !   Where R is a matrix which permutes two rows of A, and C is
    !   a matrix which permutes columns of A.
    !
    !
    !   The action of R on a vector,       (Rx)(i) = x(row(i))
    !   The action of C^-1 on a vector, (C^-1x)(i) = x(col(i))
    !
    !   Using this routine in conjuction with LSLR_SOLVE_pivot(), we
    !   first solve:
    !
    !   R A C y = R b
    !
    !   and then solve  y = C^-1 x, so that x(col(i)) = y(i)
    !
    !   Notes:
    !   1. A must be rank 3 or 4.  Singular rows are solved by
    !   setting te solution on that row to 0
    !   2. A will be overwritten with its LU decomposition.
    !   3. Singular is defined by the normalized relation:
    !   abs(pivot_value(row)) <= alittle
    !
    !
    !=======================================================================
    integer, intent(IN) :: n

    ! arguments
    integer, intent(OUT) :: row(n)
    integer, intent(OUT) :: col(n)
    real(r8), dimension(n,n), intent(INOUT) :: A

    ! local variables
    integer :: piv,i,j,rowp,colp
    real(r8) :: tmp, dont_care
  ! dont_care_magnitude sets the relative magnitude below which values
  ! are "zero" wrt a cell face's LU solve component. It use implies that
  ! components with magnitudes less than 'dont_care' are "unphysical" wrt
  ! the solution being sought and can be summarily eliminated. 
  ! The absolute value
  ! for a cell face is the product of the largest magnitude component of
  ! the (ndim+1)x(ndim+1) matrix with the "dont_care_magnitude"
    real(r8), parameter :: dont_care_magnitude=1.d-15
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do i=1,n
       row(i)=i
       col(i)=i
    enddo
    dont_care = maxval(abs(A))*dont_care_magnitude
    do piv=1,n
       ! find largest pivot, at location before piviting (rowp,colp)
       tmp=-1
       do i=piv,n
          do j=piv,n
             if ( tmp < abs(A(i,j)) ) then
                tmp = abs(A(i,j))
                rowp=i           
                colp=j
             endif    
          enddo   
       enddo   
       ! check for zero pivots.  
       if (tmp <= dont_care)A(rowp,colp)=2d200  ! component will be 0'd out in solve.
       if (piv .ne. rowp) then
          ! swap rows:  piv & row
          do j=1,n
             tmp=A(piv,j)
             A(piv,j)=A(rowp,j)
             A(rowp,j)=tmp
          enddo
          j = row(piv)
          row(piv)=row(rowp)
          row(rowp)=j
       endif
       if (piv .ne. colp) then
          ! swap columns piv and colp
          do i=1,n
             tmp=A(i,piv)
             A(i,piv)=A(i,colp)
             A(i,colp)=tmp
          enddo
          j = col(piv)
          col(piv)=col(colp)
          col(colp)=j
       endif
       ! row reduce the rest of the matrix:
       do i=piv+1,n
          tmp = A(i,piv)/A(piv,piv)
          do j=piv+1,n
             A(i,j) = A(i,j) - tmp*A(piv,j)
          enddo
          ! store information to act on B:
          A(i,piv)=tmp
       enddo
    enddo
  END SUBROUTINE LSLR_SOLVE_LU


  SUBROUTINE LSLR_SOLVE_pivot(A,B,n,nm,row,col)
    !=======================================================================
    ! Purpose:
    !
    !   Solve a n x n system of equations, where A is the output
    !   from a previous call to LSLR_SOLVE_LU and B is the RHS vector,
    !   Return the result in ans.  (data in B is trashed)
    !
    !=======================================================================
    integer, intent(IN)  :: n,nm

    ! arguments
    integer, intent(IN) :: row(n)
    integer, intent(IN) :: col(n)
    real(r8), dimension(n,n), intent(IN) :: A
    real(r8), dimension(n),   intent(INOUT) :: B

    ! function return
    real(r8) :: tmp
    real(r8), dimension(n) :: B2
    ! local variables
    integer :: piv,i,j
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! solve the system Ax=B
    ! requres LU decompostion in A, as well as row() and col()
    ! (row and column permutation arrays)
    ! apply row permuations to B
    do i=1,n 
       B2(i)=B(row(i))
    enddo    
    ! apply L to B
    do piv=1,n-1
       do i=piv+1,n
          B2(i) = B2(i) - A(i,piv)*B2(piv)
       enddo 
    enddo    
    ! apply U inverse to B.
    do i=n,1,-1
       if (A(i,i)>1d200) then  ! flag indicating 0 pivot
          B2(i)=0
       else  
          tmp = B2(i)
          do j=i+1,n 
             tmp = tmp  - A(i,j)*B2(j)
          enddo
          B2(i)=tmp/A(i,i)
       endif
    enddo
    ! apply column permutation
    do i=1,n
       B(col(i))=B2(i)
    enddo
  END SUBROUTINE LSLR_SOLVE_pivot


  FUNCTION GoodPhiSolution()
    logical, pointer, dimension(:,:) :: GoodPhiSolution
    logical, allocatable, dimension(:,:), target, save :: GPS
    if(.not. ALLOCATED(GPS))ALLOCATE(GPS(nfc,ncells))
    GPS = SolveFlag(1,:,:)
    GoodPhiSolution => GPS
  END FUNCTION GoodPhiSolution

  SUBROUTINE update_LSLR_All(UFlag)
    logical, intent(IN) :: UFlag
    UpdateDecomp = UFlag
  END SUBROUTINE update_LSLR_All

  SUBROUTINE update_LSLR_SelectCell(UFlagCell)
    logical, dimension(ncells), intent(IN) :: UFlagCell
    integer :: i,j
    do j = 1,ncells
    do i = 1,nfc
      UpdateDecomp(i,j) = UFlagCell(j)
    end do
    end do
  END SUBROUTINE update_LSLR_SelectCell


  SUBROUTINE update_LSLR_SelectFace(UFlagFace)
    logical, dimension(nfc,ncells), intent(IN) :: UFlagFace
    UpdateDecomp = UFlagFace
  END SUBROUTINE update_LSLR_SelectFace


END MODULE DISCRETE_DERIVATIVES
