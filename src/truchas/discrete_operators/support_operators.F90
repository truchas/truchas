MODULE SUPPORT_OPERATORS
  !=======================================================================
  ! Purpose(s):
  !   Define an interface to the Augustus Package, which is used to 
  !   calculate the cell-based support operator diffusion discretization 
  !   matrices.
  !
  !   Public Interface(s): CALCULATE_FLUX
  !
  ! Contains:
  !
  ! Author(s): M. Berndt (berndt@lanl.gov), adopted from code written by
  !            Travis Austin (taustin@lanl.gov).
  !            Interface to Augustus developed by Michael L. Hall.
  !
  !=======================================================================

  USE bc_data_types
  USE kind_module,        ONLY: real_kind, int_kind, log_kind
  USE parameter_module,   ONLY: ncells, nfc, ndim, nvc
  USE gs_module,          ONLY: EE_GATHER, EN_GATHER
  USE fluid_data_module,  ONLY: Cell_isnt_void
  USE UbikSolve_module
  use Augustus
  
  IMPLICIT NONE

  ! Private Module
  PRIVATE 

  ! Public Subroutines and Types
  PUBLIC :: CALCULATE_FLUX
  PUBLIC :: SO_Control, SO_Ubik

  TYPE SO_Control
     REAL(KIND=real_kind), DIMENSION(:,:,:), POINTER  :: S => Null()
     REAL(KIND=real_kind), DIMENSION(:), POINTER      :: D => Null()
     INTEGER(KIND=int_kind), DIMENSION(:,:), POINTER  :: BdyInfo => Null()
     LOGICAL(KIND=log_kind)                           :: ReInitialize_SO_Matrix = .TRUE.
     LOGICAL(KIND=log_kind), DIMENSION(:), POINTER    :: Cell_isnt_treated_as_void => Null()
     LOGICAL(KIND=log_kind), DIMENSION(:,:), POINTER  :: Ngbr_isnt_treated_as_void => Null()
  END TYPE SO_Control

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  TYPE(SO_Control), POINTER, SAVE :: Work_SO_Control_Data => Null()
  type(Ubik_control_type), save   :: SO_Ubik 

  ! The following global module variables were added by MLH 
  ! to match the Augustus interface.
  !
  ! Coordinates of the local nodes of each cell.
  real(real_kind), dimension(:,:,:), pointer, save :: &
    Coordinates_Nodes_of_Cells => null()
  ! Scalar area of each local face for all cells.
  real(real_kind), dimension(:,:), pointer, save :: &
    Area_Faces_of_Cells => null()
  ! Unit normal vector of each local face for all cells.
  real(real_kind), dimension(:,:,:), pointer, save :: &
    Unit_Normal_Faces_of_Cells => null()


CONTAINS

  SUBROUTINE CALCULATE_FLUX (Phi, Flux, Cond_Face, DirBdy, SO_Control_Data)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the face-centered flux (Gradient) from a cell-centered
    !   scalar quantity Phi on face f. 
    !
    !=======================================================================
    USE bc_operations
    USE UbikSolve_module
    use ArrayAllocate_Module,  only: ArrayCreate
    USE mesh_module,           only: Cell, Mesh, Vertex, Vrtx_Face, Vrtx_Bdy
    use pgslib_module,         only: PGSLib_Sum_Prefix, PGSLib_Global_Sum
 
    IMPLICIT NONE

    ! Arguments
    REAL(KIND = real_kind), DIMENSION(ncells)    ,  INTENT(IN)         :: Phi 
    REAL(KIND = real_kind), DIMENSION(nfc,ncells),  INTENT(IN)         :: Cond_Face
    REAL(KIND=real_kind)  , DIMENSION(nfc,ncells),  INTENT(IN)         :: DirBdy

    REAL(KIND = real_kind), DIMENSION(nfc,ncells),         INTENT(INOUT) :: Flux

    TYPE(SO_Control), POINTER:: SO_Control_Data

    ! Local Variables
    INTEGER(KIND=int_kind)                       :: d, j, f, fj
    REAL(KIND=real_kind)                         :: check
    REAL(KIND=real_kind),  DIMENSION(nfc,ncells) :: Phi_e
    REAL(KIND=real_kind),  DIMENSION(nfc*ncells) :: rhs, soln
    INTEGER(KIND=int_kind),DIMENSION(ncells)     :: Global_Cell_Number
    INTEGER(KIND=int_kind),DIMENSION(nfc,ncells) :: N
    REAL(KIND=real_kind)                         :: cond_cell
    LOGICAL(KIND=log_kind),DIMENSION(:,:),POINTER:: Ngbr_isnt_treated_as_void_tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    IF( .not. associated(Coordinates_Nodes_of_Cells) ) THEN

       ! Set up Augustus interface.
       call ARRAYCREATE (Coordinates_Nodes_of_Cells, &
                         1, ndim, 1, nvc, 1, ncells, &
                         'Array Coordinates_Nodes_of_Cells(ndim,nvc,ncells)')
       call EN_GATHER (Coordinates_Nodes_of_Cells(1,:,:), Vertex%Coord(1), &
                       BOUNDARY=Vrtx_Bdy(1)%Data)
       call EN_GATHER (Coordinates_Nodes_of_Cells(2,:,:), Vertex%Coord(2), &
                       BOUNDARY=Vrtx_Bdy(2)%Data)
       call EN_GATHER (Coordinates_Nodes_of_Cells(3,:,:), Vertex%Coord(3), &
                       BOUNDARY=Vrtx_Bdy(3)%Data)

       call ARRAYCREATE (Area_Faces_of_Cells, &
                         1, nfc, 1, ncells, &
                         'Array Area_Faces_of_Cells(nfc,ncells)')
       call ARRAYCREATE (Unit_Normal_Faces_of_Cells, &
                         1, ndim, 1, nfc, 1, ncells, &
                         'Array Unit_Normal_Faces_of_Cells(ndim,nfc,ncells)')
       do f = 1, nfc
         Area_Faces_of_Cells(f,:) = Cell(:)%Face_Area(f)
         do d = 1, ndim
           Unit_Normal_Faces_of_Cells(d,f,:) = Cell(:)%Face_Normal(d,f)
         end do
       end do

    END IF
    
    Work_SO_Control_Data => SO_Control_Data

    ! 
    ! prepare special treatment for cells with zero conductivity
    ! 
    IF ( .not. ASSOCIATED(Work_SO_Control_Data%Cell_isnt_treated_as_void) ) THEN
       ALLOCATE(Work_SO_Control_Data%Cell_isnt_treated_as_void(ncells))
    ENDIF

    Work_SO_Control_Data%Cell_isnt_treated_as_void = Cell_isnt_void
    ! declare cells with zero conductivity as void
    do j=1,ncells
       cond_cell = 0d0
       do f=1,nfc
          cond_cell = cond_cell + Cond_Face(f,j)
       enddo
       if ( cond_cell == 0d0 ) then
          Work_SO_Control_Data%Cell_isnt_treated_as_void(j) = .FALSE.
       endif
    enddo

    
    IF ( .not. ASSOCIATED(Work_SO_Control_Data%Ngbr_isnt_treated_as_void) ) THEN
       ALLOCATE(Work_SO_Control_Data%Ngbr_isnt_treated_as_void(nfc,ncells))
    ENDIF
    
    ALLOCATE(Ngbr_isnt_treated_as_void_tmp(nfc,ncells))
    DO f=1,nfc
       Ngbr_isnt_treated_as_void_tmp(f,:) = Work_SO_Control_Data%Cell_isnt_treated_as_void
    ENDDO
    call EE_GATHER(Work_SO_Control_Data%Ngbr_isnt_treated_as_void,Ngbr_isnt_treated_as_void_tmp)
    DEALLOCATE(Ngbr_isnt_treated_as_void_tmp)


    IF ( .not. ASSOCIATED(Work_SO_Control_Data%S) ) THEN
       ALLOCATE(Work_SO_Control_Data%S(nfc,nfc,ncells))
       call Augustus_Set_Support_Op_Matrix ( &
         Work_SO_Control_Data%S, Cond_Face, Area_Faces_of_Cells, &
         Unit_Normal_Faces_of_Cells, Vrtx_Face, Coordinates_Nodes_of_Cells, &
         Cell%Volume, ncells, ndim)
    ELSE
       IF (Work_SO_Control_Data%ReInitialize_SO_Matrix) THEN
         call Augustus_Set_Support_Op_Matrix ( &
           Work_SO_Control_Data%S, Cond_Face, Area_Faces_of_Cells, &
           Unit_Normal_Faces_of_Cells, Vrtx_Face, Coordinates_Nodes_of_Cells, &
           Cell%Volume, ncells, ndim)
       ENDIF
    ENDIF

    IF ( .not. ASSOCIATED(Work_SO_Control_Data%D) ) THEN
       ALLOCATE(Work_SO_Control_Data%D(nfc*ncells))
       CALL SO_Lump()
    ELSE
       IF (Work_SO_Control_Data%ReInitialize_SO_Matrix) THEN
          CALL SO_Lump()
       ENDIF
    ENDIF


    Global_Cell_Number = PGSLib_Sum_Prefix( (/ (1, j=1,ncells) /) )

    ! Initialize all of the variables
    N    = 0
    soln = 0

    WHERE(Work_SO_Control_Data%BdyInfo == 1) Flux = 0

    CALL EE_GATHER (Phi_e, Phi)


    ! Calculate difference of Phi on faces (RHS for system used to get flux).
    !   COMPUTE -Div^* C Phi
    CELL_LOOP_1: DO j=1,ncells
       FACE_LOOP_1: DO f = 1,nfc

          fj = (j-1)*nfc + f

          IF( Work_SO_Control_Data%BdyInfo(f,j) == 1 ) THEN 
             
             IF (Mesh(j)%Ngbr_Cell(f) == 0) THEN 
                
                N(f,j)  = +1 ! always outward unit normal on boundary
                rhs(fj) = N(f,j) * ( DirBdy(f,j) - Phi(j) ) 

             ELSE 

                IF (Mesh(j)%Ngbr_Cell_Orig(f) > Global_Cell_Number(j)) THEN    
                   
                   N(f,j)  = +1 ! outward unit normal
                   rhs(fj) = N(f,j) * ( Phi_e(f,j) - Phi(j) )

                ELSE

                   N(f,j)  = -1 ! inward unit normal
                   rhs(fj) =  N(f,j) * ( Phi_e(f,j) - Phi(j) )
                   
                END IF
                   
             END IF

          ELSE ! if Flux boundary condition or internal face of void

             rhs(fj)  = Flux(f,j) ! Flux is specified by boundary conditions


             IF (Mesh(j)%Ngbr_Cell(f) == 0) THEN 
                N(f,j) = +1 
             ELSE
                IF (Mesh(j)%Ngbr_Cell_Orig(f) > Global_Cell_Number(j)) THEN             
                   N(f,j) = +1 ! outward unit normal
                ELSE
                   N(f,j) = -1 ! inward unit normal
                END IF
             END IF


             ! adjust right hand side to have the correct sign, note that
             ! we need plus one on the diagonal (see below)
             If( Work_SO_Control_Data%Ngbr_isnt_treated_as_void(f,j) & 
                  .and. Mesh(j)%Ngbr_Cell(f) .ne.0 ) then
                rhs(fj) = N(f,j)*rhs(fj)
             ELSE IF ( .not. Cell_isnt_void(j) ) THEN
                rhs(fj) = N(f,j)*rhs(fj)
             ENDIF
             
          END IF

       END DO FACE_LOOP_1
    END DO CELL_LOOP_1

    ! Make sure that all values in Flux vector on same face have the same value.
    where(N == -1) Flux = -Flux

    soln = RESHAPE(Flux,(/ nfc*ncells /))
 
    check =  PGSLib_Global_Sum(rhs*rhs)
    IF( check /= 0.0 ) THEN
    
       CALL Ubik_create(SO_Ubik)

       CALL Ubik_set_residual_update(SO_Ubik)

       CALL Ubik_CG(soln, rhs, SO_Ubik, Ubik_SO_Matvec, Ubik_SO_Precond)

       Flux = RESHAPE(soln,(/ nfc, ncells /))
       
       WHERE(N == -1) Flux = -Flux

       CALL Ubik_destroy(SO_Ubik)

    ELSE
   
       Flux = 0

    END IF

    RETURN

  END SUBROUTINE CALCULATE_FLUX
  
  SUBROUTINE Ubik_SO_MatVec(X_vec, Y, status)
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutine calculates S*x = y and here it is used to pass 
    !   into the UbikSolve subroutine called Ubik_CG.F which solves for
    !   the Flux Vector in CALCULATE_FLUX
    !
    !=======================================================================
    USE UbikSolve_module

    ! Arguments
    TYPE (Ubik_vector_type),                INTENT(INOUT) :: X_vec
    REAL (real_kind), DIMENSION(:), TARGET, INTENT(INOUT) :: Y
    INTEGER (int_kind),                     INTENT(OUT)   :: status
    REAL(real_kind), DIMENSION(:), POINTER :: X => null()

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    X => Ubik_values_ptr(X_vec)
    Y = 0

    Call SO_MatVec(X, Y, status)

    RETURN

  END SUBROUTINE Ubik_SO_MatVec


  SUBROUTINE SO_MatVec(X, Y, status)
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutine calculates S*x = y and here it is used to pass 
    !   into the UbikSolve subroutine called Ubik_CG.F which solves for
    !   the Flux Vector in CALCULATE_FLUX
    !
    !=======================================================================
    USE mesh_module,              ONLY: Mesh
    USE gs_module,                ONLY: EE_GATHER
    USE UbikSolve_module
    use pgslib_module,            only: PGSLib_Sum_Prefix

    ! Arguments
    REAL (real_kind), DIMENSION(:), INTENT(IN)            :: X
    REAL (real_kind), DIMENSION(:), INTENT(OUT), TARGET   :: Y
    
    INTEGER (int_kind),                     INTENT(OUT)   :: status

    ! Local Variables
    INTEGER(KIND=int_kind)                       :: c,f,d1,d2,i,j
    REAL(KIND=real_kind)                         :: rowsum
    REAL(KIND=real_kind),DIMENSION(nfc,ncells)   :: a,b
    REAL(KIND=real_kind)                         :: fact1, fact2   
    INTEGER(KIND=int_kind),DIMENSION(ncells)     :: Global_Cell_Number
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Y = 0
    Global_Cell_Number = PGSLib_Sum_Prefix( (/ (1, j=1,ncells) /) )

    CELL_LOOP: DO c=1,ncells
       FACE_ONE_LOOP: DO d1 = 1,nfc

          i = (c-1)*nfc + d1

          IF(Mesh(c)%Ngbr_Cell_Orig(d1) > Global_Cell_Number(c) .OR. Mesh(c)%Ngbr_Cell(d1) == 0 ) THEN
             fact1 =  1.0
          ELSE
             fact1 = -1.0
          ENDIF

          IF ( Work_SO_Control_Data%Cell_isnt_treated_as_void(c) ) THEN

             rowsum = 0.0

             FACE_TWO_LOOP: DO d2 = 1,nfc

                j = (c-1)*nfc + d2
                   
                IF(Mesh(c)%Ngbr_Cell_Orig(d2) > Global_Cell_Number(c) .OR. Mesh(c)%Ngbr_Cell(d2) == 0) THEN
                   fact2 = 1.0
                ELSE
                   fact2 = -1.0
                ENDIF
                rowsum = rowsum + fact1*fact2 * Work_SO_Control_Data%S(d1,d2,c) * X(j)

             END DO FACE_TWO_LOOP
             
             IF( Work_SO_Control_Data%BdyInfo(d1,c) == 1 ) THEN
                Y(i) = rowsum
             ELSEIF( Work_SO_Control_Data%BdyInfo(d1,c) == 0 ) THEN ! Enforce boundary condition
                ! we always want a positive value on the diagonal of the matrix
                Y(i) = X(i)                
             ENDIF
          
          ELSE  ! Cell is void so set equation to identity

             Y(i) = X(i)
             
          ENDIF

       END DO FACE_ONE_LOOP
    END DO CELL_LOOP

    a = RESHAPE(Y, (/ nfc, ncells /))

    CALL EE_GATHER(b, a)

    DO c=1,ncells
       DO f=1,nfc
          IF ( ( Work_SO_Control_Data%Ngbr_isnt_treated_as_void(f,c)       &
               .OR. Work_SO_Control_Data%Cell_isnt_treated_as_void(c) )    &
               .AND. Work_SO_Control_Data%BdyInfo(f,c) == 1 )              &
               THEN 
             a(f,c) = a(f,c) + b(f,c) 
          ENDIF
       ENDDO
    ENDDO
   
    Y = RESHAPE(a,(/ nfc*ncells /))

    status = 0

    RETURN

  END SUBROUTINE SO_MatVec

  SUBROUTINE Ubik_SO_Precond(b,X_vec,status)
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutine calculates M^-1 b = x where M is a preconditioner and
    !   is used to be passed into the UbikSolve subroutine called Ubik_CG.F 
    !   which solves for the Flux Vector in CALCULATE_FLUX
    !
    !=======================================================================
    USE UbikSolve_module

    ! Arguments
    TYPE (Ubik_vector_type), INTENT(INOUT) :: X_vec
    REAL (real_kind), dimension(:), target,   INTENT(IN)    :: b
    INTEGER (int_kind),                 INTENT(OUT)   :: status
    REAL(real_kind), DIMENSION(:), POINTER            :: X => null()
   

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    X => Ubik_values_ptr(X_vec)

    X = b / Work_SO_Control_Data%D

    status = 0

    RETURN

  END SUBROUTINE Ubik_SO_Precond

  SUBROUTINE SO_Lump()
    !=======================================================================
    ! Purpose(s):
    ! 
    !   Subroutine calculates S*x = y and here it is used to pass 
    !   into the UbikSolve subroutine called Ubik_CG.F which solves for
    !   the Flux Vector in CALCULATE_FLUX
    !
    !=======================================================================
    USE gs_module,                ONLY: EE_GATHER
    USE UbikSolve_module
    use pgslib_module,            only: PGSLib_Sum_Prefix

    ! Arguments

    ! Local Variables
    REAL(KIND=real_kind),DIMENSION(nfc,ncells)   :: a,b
    INTEGER(KIND=int_kind)                       :: c,f,d1,i,j
    REAL(KIND=real_kind),DIMENSION(ncells*nfc)   :: Y
    INTEGER(KIND=int_kind),DIMENSION(ncells)     :: Global_Cell_Number
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Y = 0
    Global_Cell_Number = PGSLib_Sum_Prefix( (/ (1, j=1,ncells) /) )

    CELL_LOOP: DO c=1,ncells
       FACE_ONE_LOOP: DO d1 = 1,nfc

          i = (c-1)*nfc + d1

          IF ( Work_SO_Control_Data%Cell_isnt_treated_as_void(c) ) THEN

             IF( Work_SO_Control_Data%BdyInfo(d1,c) == 1 ) THEN
                Y(i) = Work_SO_Control_Data%S(d1,d1,c)
             ELSEIF( Work_SO_Control_Data%BdyInfo(d1,c) == 0 ) THEN ! Enforce boundary condition
                ! we always want a positive value on the diagonal of the matrix
                Y(i) = 1.0   
             ENDIF
          
          ELSE  ! Cell is void so set equation to identity

             Y(i) = 1.0
             
          ENDIF

       END DO FACE_ONE_LOOP
    END DO CELL_LOOP

    a = RESHAPE(Y, (/ nfc, ncells /))

    CALL EE_GATHER(b, a)

    DO c=1,ncells
       DO f=1,nfc
          IF ( ( Work_SO_Control_Data%Ngbr_isnt_treated_as_void(f,c)     &
               .OR. Work_SO_Control_Data%Cell_isnt_treated_as_void(c) )  &
               .AND. Work_SO_Control_Data%BdyInfo(f,c) == 1 )            &
               THEN
             a(f,c) = a(f,c) + b(f,c) 
          ENDIF
       ENDDO
    ENDDO
   
    Work_SO_Control_Data%D = RESHAPE(a,(/ nfc*ncells /))

    RETURN

  END SUBROUTINE SO_Lump

END MODULE SUPPORT_OPERATORS
