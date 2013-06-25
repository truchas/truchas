MODULE Y_EQ_AX_PRS
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures necessary to perform a matrix-vector multiply
  !   (Ax), returning it in vector y
  !
  !   Public Interface:
  !
  !     * call Y_EQ_AX_PRESSURE (X, Y, status)
  !
  !         Returns in y the vector Ax, where A represents a discrete
  !         operator for a pressure projection equation del*Q = S, where
  !         Q is a "pressure flux" equal to (1/rho)*Grad(P), with rho
  !         the density and P a pressure.
  !
  ! Contains: Y_EQ_AX_PRESSURE
  !
  !=======================================================================

  Implicit None
  Private

  ! public procedures
  Public :: Y_EQ_AX_PRESSURE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE Y_EQ_AX_PRESSURE (X_vec, Y, status)
    !=======================================================================
    ! Purpose:
    !
    !   Compute y = Ax where A represents a discrete operator for a
    !   pressure projection equation del*Q = S, where Q is a "pressure
    !   flux", Q = 1/rho Grad(P), rho is the density and P is a pressure.
    !   S is a source term equation to D/dt for a regular pressure
    !   projection and D for a MAC projection, where D is a velocity
    !   divergence.
    !=======================================================================
    use bc_module,              only: BC_Prs
    use bc_operations
    use constants_module,       only: zero
    use do_interface,           only: DO_Specifier, do_init_ss,         &
                                      do_gradient_face, DO_SOLVE_ORTHO, &
                                      DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: fluidRho, Solid_Face, isPureImmobile, &
                                      MinFaceFraction, MinFluidRho
    use kind_module,            only: int_kind, real_kind
    use mesh_module,            only: Cell
    use parameter_module,       only: ncells, ndim, nfc
    use projection_data_module, only: Coeff, Vol_over_RhoCsqDt
    use time_step_module,       only: dt
    use timing_tree
    use UbikSolve_module
    use support_operators ,     only: CALCULATE_FLUX
    !use fluid_data_module,      only: fluidRho
    use ff_discrete_ops_data,   only: use_ff_support_operators, FF_SO_Control_Data

    implicit none

    ! Arguments
    type (Ubik_vector_type),                intent(INOUT) :: X_vec
    real (real_kind), dimension(:), target, intent(INOUT) :: Y
    integer (int_kind),                     intent(OUT)   :: status
    
    
    ! Local Variables
    type(DO_Specifier),pointer,save                  :: YPressure_SS =>NULL()
    integer(KIND=int_kind)                           :: YPSolveTech
    integer(KIND = int_kind)                         :: f, c, d
    real(KIND = real_kind), dimension(nfc, ncells)   :: N_Dot_Grad_X
    integer(int_kind)                                :: BdyCell, BdyFace, BdyPt
    type (BC_Operator),                   POINTER    :: EXTERIOR_Operator
    type (BC_Atlas),                      POINTER    :: EXTERIOR_Atlas
    type (BC_Operator),                   POINTER    :: DIRICHLET_Operator
    type (BC_Atlas),                      POINTER    :: DIRICHLET_Atlas
    integer( int_kind),     dimension(:), POINTER    :: BdyFaceList
    integer( int_kind),     dimension(:), POINTER    :: BdyCellList
    real(KIND = real_kind), dimension(ndim,nfc,ncells) :: fgradX
    real(real_kind), dimension(:), pointer :: X
    
    real(real_kind), dimension(nfc,ncells)               :: Flux
    real(real_kind), dimension(nfc,ncells)               :: SO_FluidRho_Face
    real(real_kind), dimension(nfc,ncells)               :: BC_Prs_zero
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    X => Ubik_values_ptr(X_vec)
    
    if ( .not. use_ff_support_operators ) then

       if(.not. ASSOCIATED(YPressure_SS))then
          YPSolveTech=DO_SOLVE_LU_LSLR; if(use_ortho_face_gradient)YPSolveTech=DO_SOLVE_ORTHO
          call do_init_ss(YPressure_SS,SOLVETECH=YPSolveTech,BC_SPEC=Pressure_BC)
       endif

    endif

    ! Start the timer.
    call start_timer ("Timer_Solver_TMP2")

    Y = zero
    
    if ( use_ff_support_operators ) then

       Flux = zero
       BC_Prs_zero = zero
              
       SO_FluidRho_Face = zero
       DO f=1,nfc
          where (FluidRho(:) /= zero) SO_FluidRho_Face(f,:) = 1.0/FluidRho(:)
       ENDDO
          
          
       CALL CALCULATE_FLUX(PHI = X,                &
            FLUX = Flux,                           &
            COND_FACE = SO_FluidRho_Face,          &
            DIRBDY = BC_Prs,                       &
            SO_CONTROL_DATA = FF_SO_Control_Data )
          
       
    else
       
       fgradx = zero
       
       call DO_GRADIENT_FACE(PHI=X, SOLVESPEC=YPressure_SS, GRAD=fgradx)

       ! Now put in the Exterior BCs.  At the moment, Exterior Faces get zero gradient.
       
       EXTERIOR_Operator => BC_Spec_Get_Operator(Pressure_BC, BC_EXTERIOR_Op)
       EXTERIOR_Atlas    => BC_OP_Get_Atlas(EXTERIOR_Operator)
       
       BdyCellList => BC_Get_Cell(EXTERIOR_Atlas)
       BdyFaceList => BC_Get_Face(EXTERIOR_Atlas)
       
       EXT_BDY_LOOP: do BdyPt = 1, DATA_SIZE(EXTERIOR_Atlas)
          BdyCell = BdyCellList(BdyPt)
          BdyFace = BdyFaceList(BdyPt)
          FGradX(:,BdyFace,BdyCell) = 0.0
       end do EXT_BDY_LOOP
       
       ! For all cells, dot the solution gradient with the face unit normal.
       do c = 1, ncells
          do f = 1, nfc
             N_Dot_Grad_X(f,c) = zero
             do d = 1,ndim
                N_Dot_Grad_X(f,c) = N_Dot_Grad_X(f,c) + fgradX(d,f,c)*Cell(c)%Face_Normal(d,f)
             end do
          end do
       end do
       
       ! If Dirichlet BCs are on this face, then 
       ! use the input solution value rather than the gradient.
       DIRICHLET_Operator => BC_Spec_Get_Operator(Pressure_BC, BC_DIRICHLET_Op)
       DIRICHLET_Atlas    => BC_OP_Get_Atlas(DIRICHLET_Operator)
       
       BdyCellList => BC_Get_Cell(DIRICHLET_Atlas)
       BdyFaceList => BC_Get_Face(DIRICHLET_Atlas)
       
       DIR_BDY_LOOP: do BdyPt = 1, DATA_SIZE(DIRICHLET_Atlas)
          BdyCell = BdyCellList(BdyPt)
          BdyFace = BdyFaceList(BdyPt)
          N_Dot_Grad_X(BdyFace, BdyCell) = X(BdyCell)
       end do DIR_BDY_LOOP

    endif


   ! Accumulate the contribution for the faces.
   do c = 1, ncells
      ! Void and Solid cells return X
      if (fluidRho(c) == zero .or. isPureImmobile(c)) then
         Y(c) = X(c) / (MinFaceFraction * MinFluidRho * Cell(c)%Volume**0.66666)
         CYCLE
      end if

      if ( use_ff_support_operators ) then

         ! If the cell is not pure immobile or void, then loop over faces
         do f = 1, nfc
            ! where the face is in pure solid, assume N_Dot_Grad_X is zero
            if (.NOT. Solid_Face(f,c)) then
               Y(c) = Y(c) + Flux(f,c)
            end if
         end do

      else 
   
         ! If the cell is not pure immobile or void, then loop over faces
         do f = 1, nfc
            ! where the face is in pure solid, assume N_Dot_Grad_X is zero
            if (.NOT. Solid_Face(f,c)) then
               Y(c) = Y(c) + Coeff(f,c) * N_Dot_Grad_X(f,c)
            end if
         end do

      endif

      Y(c) = (Y(c) - Vol_over_RhoCsqDt(c)*X(c)/dt) / Cell(c)%Volume

   end do

   call stop_timer ("Timer_Solver_TMP2")
   
   status = 0
   return

   END SUBROUTINE Y_EQ_AX_PRESSURE


END MODULE Y_EQ_AX_PRS
