!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    use bc_operations
    use do_interface,           only: DO_Specifier, do_init_ss,         &
                                      do_gradient_face, DO_SOLVE_ORTHO, &
                                      DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: fluidRho, Solid_Face, isPureImmobile, &
                                      MinFaceFraction, MinFluidRho
    use kinds, only: r8
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell
    use projection_data_module, only: Coeff, Vol_over_RhoCsqDt
    use time_step_module,       only: dt
    use truchas_timers
    use UbikSolve_module

    ! Arguments
    type (Ubik_vector_type), intent(INOUT) :: X_vec
    real(r8), dimension(:), target, intent(INOUT) :: Y
    integer, intent(OUT) :: status
    
    ! Local Variables
    type(DO_Specifier), pointer, save :: YPressure_SS =>NULL()
    integer :: YPSolveTech
    integer :: f, c, d
    real(r8), dimension(nfc, ncells) :: N_Dot_Grad_X
    integer :: BdyCell, BdyFace, BdyPt
    type(BC_Operator), POINTER :: EXTERIOR_Operator
    type(BC_Atlas),    POINTER :: EXTERIOR_Atlas
    type(BC_Operator), POINTER :: DIRICHLET_Operator
    type(BC_Atlas),    POINTER :: DIRICHLET_Atlas
    integer, dimension(:), POINTER :: BdyFaceList
    integer, dimension(:), POINTER :: BdyCellList
    real(r8), dimension(ndim,nfc,ncells) :: fgradX
    real(r8), dimension(:), pointer :: X
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    call start_timer ('ppe-matvec')

    X => Ubik_values_ptr(X_vec)
    
    if(.not. ASSOCIATED(YPressure_SS))then
       YPSolveTech=DO_SOLVE_LU_LSLR; if(use_ortho_face_gradient)YPSolveTech=DO_SOLVE_ORTHO
       call do_init_ss(YPressure_SS,SOLVETECH=YPSolveTech,BC_SPEC=Pressure_BC)
    endif

    ! Start the timer.
    !call start_timer ("Timer_Solver_TMP2")

    Y = 0
    
    call start_timer ('ppe-matvec-face-grad')
    fgradx = 0

    call DO_GRADIENT_FACE(PHI=X, SOLVESPEC=YPressure_SS, GRAD=fgradx)
    call stop_timer ('ppe-matvec-face-grad')

    ! Now put in the Exterior BCs.  At the moment, Exterior Faces get zero gradient.

    call start_timer ('ppe-matvec-bc')
    EXTERIOR_Operator => BC_Spec_Get_Operator(Pressure_BC, BC_EXTERIOR_Op)
    EXTERIOR_Atlas    => BC_OP_Get_Atlas(EXTERIOR_Operator)

    BdyCellList => BC_Get_Cell(EXTERIOR_Atlas)
    BdyFaceList => BC_Get_Face(EXTERIOR_Atlas)

    EXT_BDY_LOOP: do BdyPt = 1, DATA_SIZE(EXTERIOR_Atlas)
       BdyCell = BdyCellList(BdyPt)
       BdyFace = BdyFaceList(BdyPt)
       FGradX(:,BdyFace,BdyCell) = 0.0
    end do EXT_BDY_LOOP
    call stop_timer ('ppe-matvec-bc')

    ! For all cells, dot the solution gradient with the face unit normal.
    call start_timer ('ppe-matvec-flux')
    do c = 1, ncells
       do f = 1, nfc
          N_Dot_Grad_X(f,c) = 0
          do d = 1,ndim
             N_Dot_Grad_X(f,c) = N_Dot_Grad_X(f,c) + fgradX(d,f,c)*Cell(c)%Face_Normal(d,f)
          end do
       end do
    end do
    call stop_timer ('ppe-matvec-flux')

    ! If Dirichlet BCs are on this face, then 
    ! use the input solution value rather than the gradient.
    call start_timer ('ppe-matvec-bc')
    DIRICHLET_Operator => BC_Spec_Get_Operator(Pressure_BC, BC_DIRICHLET_Op)
    DIRICHLET_Atlas    => BC_OP_Get_Atlas(DIRICHLET_Operator)

    BdyCellList => BC_Get_Cell(DIRICHLET_Atlas)
    BdyFaceList => BC_Get_Face(DIRICHLET_Atlas)

    DIR_BDY_LOOP: do BdyPt = 1, DATA_SIZE(DIRICHLET_Atlas)
       BdyCell = BdyCellList(BdyPt)
       BdyFace = BdyFaceList(BdyPt)
       N_Dot_Grad_X(BdyFace, BdyCell) = X(BdyCell)
    end do DIR_BDY_LOOP
    call stop_timer ('ppe-matvec-bc')

   ! Accumulate the contribution for the faces.
   call start_timer ('ppe-matvec-div')
   do c = 1, ncells
      ! Void and Solid cells return X
      if (fluidRho(c) == 0 .or. isPureImmobile(c)) then
         Y(c) = X(c) / (MinFaceFraction * MinFluidRho * Cell(c)%Volume**0.66666)
         CYCLE
      end if

      ! If the cell is not pure immobile or void, then loop over faces
      do f = 1, nfc
         ! where the face is in pure solid, assume N_Dot_Grad_X is zero
         if (.NOT. Solid_Face(f,c)) then
            Y(c) = Y(c) + Coeff(f,c) * N_Dot_Grad_X(f,c)
         end if
      end do

      Y(c) = (Y(c) - Vol_over_RhoCsqDt(c)*X(c)/dt) / Cell(c)%Volume

   end do
   call stop_timer ('ppe-matvec-div')

   !call stop_timer ("Timer_Solver_TMP2")
   call stop_timer ('ppe-matvec')
   
   status = 0

   END SUBROUTINE Y_EQ_AX_PRESSURE

END MODULE Y_EQ_AX_PRS
