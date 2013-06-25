MODULE PROJECTION_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures cell-face and cell-centered velocity 
  !   projection operations.
  !
  !   Public Interface:
  !
  !     * call PROJECTION ()
  !
  !         Projection phase of the incompressible Navier-Stokes
  !         solution algorithm.
  !
  ! Contains: PROJECTION
  !
  !           MAC_PROJECTION
  !           VELOCITY_TO_FACES
  !           INTERPOLATE_VELOCITY_TO_FACES
  !           MAC_CORRECTOR
  !           MAC_RHS
  !           NONORTHO_PROJECTION
  !           ORTHO_STENCIL
  !           PROJECTION_CORRECTOR
  !           VC_DIVERGENCE_NORMS
  !           VF_DIVERGENCE_NORMS
  !
  !=======================================================================

  use kind_module, only: real_kind, int_kind

  implicit none

  ! Private Module
  private

  ! Public Subroutines
  public :: PROJECTION, ORTHO_STENCIL 

  ! Power for averaging of sound speed by material
  integer(KIND = int_kind), parameter            :: navg = 1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE PROJECTION ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Projection phase of the incompressible Navier-Stokes solution
    !   algorithm.  Begins by extrapolating time '*' cell-centered 
    !   velocities to faces, applies a MAC_PROJECTION to yield a 
    !   solenoidal face-velocity field via a cell-centered pressure, and
    !   concludes by applying the pressure gradient to the cell-centered
    !   velocities.
    !=======================================================================
    use bc_operations
    use body_data_module,       only: body_force_face_method
    use body_force_module,      only: add_face_body_force,    &
                                      compute_gravityhead
    use constants_module,       only: zero
    use cutoffs_module,         only: alittle
    use discrete_ops_data,      only: use_ortho_face_gradient
    use ff_discrete_ops_data,   only: use_ff_support_operators
    use fluid_data_module,      only: fluidRho, Solid_Face, IsPureImmobile, &
                                      Rho_Face, Rho_Face_n, Fluxing_Velocity
    use fluid_utilities_module, only: FLUIDDENSITYFACE
    use gs_module,              only: EE_GATHER
    use kind_module,            only: int_kind, real_kind
    use mesh_module,            only: Mesh, Cell
    use parameter_module,       only: ncells, ndim, nfc
    use projection_data_module, only: dt_gradP_over_Rho, Vol_over_RhoCsqDt, &
                                      ghc, ghn, dtRhoG_over_Rho, &
                                      Fcsf_new,dtCsf_over_Rho
    use surface_tension_module, only: CSF_FACE, csf_normal
    use time_step_module,       only: dt
    use timing_tree
    use truchas_logging_services

    use fischer_module

    implicit none
    
    ! Argument List

    ! Local Variables
    integer                                              :: status
    integer(int_kind)                                    :: n, f
    real(real_kind),   dimension(:,:), allocatable, save :: Scalar_Cell_Face
    real(real_kind),   dimension(:,:), allocatable, save :: Scalar_Cell_Ngbr
    logical, save                                        :: first_time = .true.
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the projection timer.
    call start_timer("timer_projection")

    ALLOCATE (Scalar_Cell_Face(nfc,ncells),       &
              Scalar_Cell_Ngbr(nfc,ncells),       &
              dt_gradP_over_Rho(ndim,nfc,ncells), &
              Vol_over_RhoCsqDt(ncells),          &
              Rho_Face(nfc,ncells),               STAT = status)
    if (status /= 0) call TLS_panic ('PROJECTION: allocation failed')

    if (.not. use_ortho_face_gradient  .or. use_ff_support_operators ) then
       ALLOCATE (dtRhoG_over_Rho(ndim,nfc,ncells), STAT = status)
       if (status /= 0) call TLS_panic ('PROJECTION: dtRhoG_over_Rho allocation failed')
    endif

    if (use_ortho_face_gradient) then
       ALLOCATE(ghc(nfc,ncells), &
                ghn(nfc,ncells), STAT = status)
       if (status /= 0) call TLS_panic ('PROJECTION: ghc or ghn allocation failed')
    end if

    if (csf_normal) then
       ALLOCATE(dtCsf_over_Rho(ndim,nfc,ncells), &
                Fcsf_new(ndim,nfc,ncells), STAT = status)
       if (status /= 0) call TLS_panic ('PROJECTION: dtCsf_over_Rho or Fcsf_new allocation failed')
    endif

    if (first_time) then
       call FISCHER_INITIALIZE()
       first_time = .false.
    endif

    ! calculate gravity head...
    if (use_ortho_face_gradient) then
       if (body_force_face_method) then
          call COMPUTE_GRAVITYHEAD()
       else
          ghc = 0.0
          ghn = 0.0
       endif
    end if

    ! Calculate Face Value of the Fluid Density
    call FluidDensityFace
        
    if (csf_normal) then
       call CSF_FACE(Fcsf_new)
    endif

    ! Transfer Cell Centered Velocity to Faces
    !  (use the Rhie-Chow algorithm for the pressure and buoyant terms)
    call VELOCITY_TO_FACES()
      
    if (csf_normal) then 
       do f=1,nfc
          do n=1,ndim 
             ! particular case when void cells. density = zero
             where (Rho_Face(f,:) <= alittle) 
                dtCsf_over_Rho(n,f,:)=zero
             elsewhere 
                dtCsf_over_Rho(n,f,:)=dt*Fcsf_new(n,f,:)/Rho_Face(f,:)
             endwhere
          enddo
       enddo
    endif

    if (.not. use_ortho_face_gradient  .or. use_ff_support_operators ) then
       ! Add gravitational acceleration to face velocities and compute Fluxing_Velocity
       call ADD_FACE_BODY_FORCE(dt, dtRhoG_over_Rho)
       do f = 1, nfc
          do n = 1, ndim
             Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + &
                                 dtRhoG_over_Rho(n,f,:)*Cell(:)%Face_Normal(n,f)
          end do
       end do
    end if

    ! Correct the Fluxing_Velocity's with a MAC projection to make them solenoidal.
    call MAC_PROJECTION(Fluxing_Velocity)

    ! Void cells don't have physical (nor solenoidal) Fluxing_Velocity's
    ! Set Fluxing_Velocity on faces that adjoin 'real' fluid cells to the 
    ! corresponding velocity on the other side of the face (corrected for direction)
    call EE_GATHER(Scalar_Cell_Face, -Fluxing_Velocity)
    do f = 1,nfc
       where (FluidRho == zero .AND. &
              Mesh%Ngbr_cell(f) /= 0 .AND. .not.IsPureImmobile) &
            Fluxing_Velocity(f,:) = Scalar_Cell_Face(f,:)
    end do

    ! Zero out Fluxing_Velocity's at solid faces.
    where (Solid_Face) Fluxing_Velocity = zero

    ! Apply the solenoidal correction to the cell-centered predictor velocities.
    call PROJECTION_CORRECTOR()
    
    ! Save old time face densities for the next step Velocity_to_Faces
    Rho_Face_n = Rho_Face

    ! Deallocate temporary storage
    DEALLOCATE (Scalar_Cell_Face)
    DEALLOCATE (Scalar_Cell_Ngbr)
    DEALLOCATE (dt_gradP_over_Rho)
    DEALLOCATE (Vol_over_RhoCsqDt)
    DEALLOCATE (Rho_Face) 

    if (.not. use_ortho_face_gradient .or. use_ff_support_operators) then
       DEALLOCATE (dtRhoG_over_Rho)
    endif
    
    if (use_ortho_face_gradient) then
       DEALLOCATE (ghc)
       DEALLOCATE (ghn)
    end if

    if (csf_normal) then
       DEALLOCATE (Fcsf_new)
       DEALLOCATE (dtCsf_over_Rho)
    endif
    
    ! Stop the projection timer
    call stop_timer("timer_projection")

    return

  END SUBROUTINE PROJECTION

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE MAC_PROJECTION (Fluxing_Velocity)
    !=======================================================================
    ! Purpose(s):
    !
    !   Project out the non-solenoidal part of the Fluxing_Velocity's
    !   using a MAC projection.  The source term for this Pressure-Poisson
    !   solve is the divergence of the Fluxing_Velocity's and the solution
    !   is a cell-centered pressure.  This is called a MAC projection
    !   because the pressure-velocity staggering is reminiscent of the
    !   original MAC scheme. Upon return, Fluxing_Velocity's have been
    !   corrected so that they are discretely solenoidal.
    !=======================================================================
    use bc_module,              only: BC_Vel
    use constants_module,       only: zero
    use fluid_data_module,      only: fluidRho, Solid_Face, &
                                      void_pressure, IsPureImmobile,  &
                                      Rho_Face, IsImmobile
    use kind_module,            only: int_kind, log_kind, real_kind
    use matl_module,            only: Matl
    use linear_solution,        only: Ubik_user
    use mesh_module,            only: Cell
    use parameter_module,       only: ncells, ndim, nfc, ncells_tot, nmat, mat_slot
    use pgslib_module,          only: PGSLIB_GLOBAL_SUM, PGSLIB_GLOBAL_ANY
    use projection_data_module, only: mac_projection_iterations,    &
                                      mac_projection_precond_iter,  &
                                      Boundary_Flag, UBIK_PRESSURE,   &
                                      Face_Density,                 &
                                      dirichlet_pressure,           &
                                      Vol_over_RhoCsqDt
    use property_data_module,   only: Sound_Speed
    use time_step_module,       only: dt
    use timing_tree
    use zone_module,            only: Zone
    use UbikSolve_module
    use truchas_logging_services
    implicit none

    ! Argument List
    real(real_kind), dimension(nfc,ncells), intent(INOUT) :: Fluxing_Velocity

    ! Local Variables
!    integer :: status
    logical(log_kind) :: Void_Cell_Found
    integer(int_kind) :: f, i, m, n, s, status
    real(real_kind),   dimension(:),   allocatable :: RHS, Solution, boundary_fv

#if ( defined(DARWIN_NAG_COMPILER_WORKAROUND) )
    integer :: nc
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the MAC projection timer.
    call start_timer("timer_projection_mac")
    ALLOCATE (RHS(ncells),      &
              Solution(ncells), &
              boundary_fv(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_PROJECTION: allocation failed')

    ! Zero out Fluxing_Velocity's at solid faces.
    where (Solid_Face) Fluxing_Velocity = zero

    ! Compute cell compressibility expression
    Vol_over_RhoCsqDt = zero
    ! Temporary use of RHS
    RHS = zero
    ! Temporarily store averaged sound speed in Vol_over_RhoCsqDt
    do m = 1,nmat
        if(.not.isImmobile(m)) then
           do s = 1, mat_slot
              where (Matl(s)%Cell%Id == m) &
                    RHS = RHS +  Matl(s)%Cell%Vof**navg
           end do
        endif
        ! calculate the contribution of material m to the reciprocal sound speed square
        if(Sound_Speed(m) > zero) then
            do s = 1, mat_slot
#if ( defined(DARWIN_NAG_COMPILER_WORKAROUND) )
               do nc = 1, ncells
                  if ( Matl(s)%Cell(nc)%Id == m) then
                    Vol_over_RhoCsqDt(nc) = Vol_over_RhoCsqDt(nc) + Matl(s)%Cell(nc)%Vof**navg/sound_speed(m)**2
                 end if
              end do
#else                     
               where (Matl(s)%Cell%Id == m) &
                    Vol_over_RhoCsqDt = Vol_over_RhoCsqDt + Matl(s)%Cell%Vof**navg/sound_speed(m)**2
#endif
            end do
        endif
    enddo
    ! Normalize by the total fluid fraction, Multiply by the Cell Volume and divide by dt
    ! (This puts the array in the most efficient form for Y_EQ_AX_PRS)
    where( RHS > zero )  Vol_over_RhoCsqDt = Cell%Volume*Vol_over_RhoCsqDt / (dt*RHS)
    RHS = zero
    ! Now divide by the fluid density
    where(FluidRho > zero ) Vol_over_RhoCsqDt = Vol_over_RhoCsqDt / FluidRho

    Face_Density = Rho_Face

    ! Get the RHS for this system
    call MAC_RHS (Fluxing_Velocity, RHS)

    do f = 1,nfc
       
       where (Boundary_Flag(f,:) == 0) &
          RHS = RHS - Fluxing_velocity(f,:)*Cell%Face_Area(f)

       boundary_fv = zero
       do n = 1,ndim
          where (Boundary_Flag(f,:) == 2) &
             boundary_fv = boundary_fv + BC_Vel(n,f,:)*Cell%Face_Normal(n,f)
       end do
       where (Boundary_Flag(f,:) == 2) &
          RHS = RHS + (boundary_fv - Fluxing_velocity(f,:))*Cell%Face_Area(f)
    end do

    ! set RHS to zero in void cells because we are solving for the change in pressure
    ! (RHS divided by dt in NONORTHO_PROJECTION)
    ! Scaling of the RHS  requires multiplying by dt here for 
    ! global scaling of RHS by dt/Volume.
    where (FluidRho(:) == zero .or. isPureImmobile(:)) RHS(:) = (void_pressure-Zone(:)%P)

    call NONORTHO_PROJECTION (RHS, Face_Density, Solution, "timer_projection_solver")

    ! Store the number of iterations.
    mac_projection_iterations   = Ubik_iter(Ubik_user(UBIK_PRESSURE)%control)
    mac_projection_precond_iter = Ubik_user(UBIK_PRESSURE)%precond_iter

    ! Apply the solenoidal correction to the face-centered velocities.
    call MAC_CORRECTOR (Solution, Face_Density, Fluxing_Velocity)
        
    ! Computed divergence norms based on the correct fluxing velocity.
    call VF_DIVERGENCE_NORMS (Fluxing_Velocity)

    ! If no dirichlet pressures specified, normalize the solution by
    ! subtracting off the global mean.
    Void_Cell_Found = .false.
    do i= 1,ncells
        if(FluidRho(i) == zero .and. .not.IsPureImmobile(i)) then
            Void_Cell_Found = .true.
            exit
        endif
    enddo

    if (.not.dirichlet_pressure .AND. .not.PGSLIB_GLOBAL_ANY(Void_Cell_Found) ) then
       Solution = (Solution - PGSLib_Global_SUM(Solution)/ncells_tot)
    end if

    ! Update the solution to the total pressure
    Zone%P = Zone%P + Solution

    DEALLOCATE (RHS)
    DEALLOCATE (Solution)
    DEALLOCATE (boundary_fv)

    ! Stop the MAC projection timer.
    call stop_timer("timer_projection_mac")

    return

  END SUBROUTINE MAC_PROJECTION

  SUBROUTINE VELOCITY_TO_FACES()
    !=======================================================================
    ! Purpose(s):
    !
    !   Transfer the cell centered velocity to cell faces
    !      Use the Rhie-Chow algorithm for the pressure gradient and buoyant forces
    !    Jim Sicilian,  11/2002.
    !=======================================================================
    use bc_module,              only: BC_Prs
    use bc_data_module,         only: BC_Pressure, BC_Zero
    use bc_operations,          only: Pressure_BC
    use constants_module,       only: zero
    use do_interface,           only: DO_Specifier, do_init_ss, &
                                      do_gradient_face, DO_SOLVE_ORTHO, &
                                      DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: fluidRho, Rho_Face, Solid_Face, Fluxing_Velocity
    use kind_module,            only: int_kind, real_kind, log_kind
    use mesh_module,            only: Cell
    use parameter_module,       only: ncells, ndim, nfc
    use projection_data_module, only: dirichlet_pressure, Boundary_Flag, ghc, ghn, &
                                      Fcsf_new
!-mf jan04
    use time_step_module,       only: dt
    use zone_module,            only: Zone
    use discrete_op_module,     only: DYNAMIC_PRESSURE_FACE_GRADIENT
    use surface_tension_module, only: csf_normal
    use support_operators ,     only: CALCULATE_FLUX
    use ff_discrete_ops_data,   only: use_ff_support_operators, FF_SO_Control_Data
    use truchas_logging_services
 
    ! local variables
    integer(int_kind)                                    :: n, f, status
    logical(log_kind), dimension(:),   allocatable       :: Mask
    real(real_kind),   dimension(:),   allocatable       :: Scalar_Cell_Center
    real(real_kind),   dimension(:),   allocatable       :: dt_over_Rho, Grad_Dot_N
    real(real_kind),   dimension(:,:), allocatable       :: Grad
    real(real_kind),   dimension(:,:,:), allocatable     :: Gradient
    integer(int_kind)                                    :: PSolveTech
    type(DO_Specifier),pointer,save                      :: Projection_SS2 =>NULL()
    real(real_kind), dimension(nfc,ncells)               :: Flux
    real(real_kind), dimension(nfc,ncells)               :: SO_FluidRho_Face

    
   ! First, allocate temporary array space
   ALLOCATE (Mask(ncells),                 &
             Scalar_Cell_Center(ncells),   &
             Gradient(ndim,nfc,ncells),    &
             Grad(ndim,ncells),            &
             dt_over_Rho(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VELOCITY_TO_FACES: allocation failed')

   if ( .not. use_ff_support_operators ) then
      ALLOCATE (Grad_Dot_N(ncells), STAT = status)
      if (status /= 0) call TLS_panic ('VELOCITY_TO_FACES: Grad_Dot_N(ncells) allocation failed')
   else
      ! This should be the first place FF_SO_Control_Data is encountered, so allocate it
      IF ( .not. ASSOCIATED(FF_SO_Control_Data) ) THEN
         ALLOCATE(FF_SO_Control_Data)
      ENDIF
      
      IF ( .not. ASSOCIATED(FF_SO_Control_Data%BdyInfo) ) THEN
         ALLOCATE(FF_SO_Control_Data%BdyInfo(nfc,ncells))
      ENDIF
   endif
   
    ! Begin by extrapolating the time '*' cell-centered velocity to the faces

    call INTERPOLATE_VELOCITY_TO_FACES

    ! Use the real Pressures at Dirichlet boundaries in this pressure gradient calculation
    BC_Prs => BC_Pressure
    
    if ( .not. use_ff_support_operators ) then

       if(.not.ASSOCIATED(Projection_SS2))then
          PSolveTech=DO_SOLVE_LU_LSLR; if(use_ortho_face_gradient)PSolveTech=DO_SOLVE_ORTHO
          call do_init_ss(Projection_SS2,SOLVETECH=PSolveTech,BC_SPEC=Pressure_BC)
       endif
       
       Scalar_Cell_Center = Zone%P
       
       if (.not. use_ortho_face_gradient) then
          ! calculate the gradient of the total pressure.
          ! Now add in the explicit pressure gradient and buoyant term at this face
          call DO_GRADIENT_FACE(PHI=Scalar_Cell_Center, SOLVESPEC=Projection_SS2, GRAD=Gradient)
       else
          ! calculate the gradient of the dynamic pressure
          ! ie; the total pressure less the gravity  head. 
          call DYNAMIC_PRESSURE_FACE_GRADIENT(Gradient, Scalar_Cell_Center, ghc, ghn)
       end if
       
    else
       
       Scalar_Cell_Center = Zone%P

       SO_FluidRho_Face = zero
       DO f=1,nfc
          where (FluidRho(:) /= zero) &
               SO_FluidRho_Face(f,:) = 1.0/FluidRho(:)
       ENDDO
       
       Flux = zero
       
       ! where there is no boundary condition or dirichlet pressure bc
       FF_SO_Control_Data%BdyInfo = 1 

       ! where there is solid face
       where (Boundary_Flag == 0) FF_SO_Control_Data%BdyInfo = 0 
       
       ! or a dirichlet velocity bc
       where (Boundary_Flag == 2) FF_SO_Control_Data%BdyInfo = 0 

       ! reinitialize the mass matrix
       FF_SO_Control_Data%ReInitialize_SO_Matrix = .TRUE. 
       
       CALL CALCULATE_FLUX(PHI = Scalar_Cell_Center,      &
            FLUX = Flux,                                  &
            COND_FACE = SO_FluidRho_Face,                 &
            DIRBDY = BC_Prs,                              &
            SO_CONTROL_DATA = FF_SO_Control_Data)
       
       ! do not reinitialize until we get back here
       FF_SO_Control_Data%ReInitialize_SO_Matrix = .FALSE.


       DO f=1,nfc
          WHERE (Cell%Face_Area(f) .ne. 0d0) 
             Flux(f,:) = Flux(f,:) / Cell%Face_Area(f)
          ENDWHERE
       END DO
       
    endif


    FACE_LOOP : do f = 1,nfc   

       ! Store the inverse face density. (needed for the surface tension model)
       dt_over_Rho = zero
       where (Rho_Face(f,:) /= zero) &
            dt_over_Rho = dt/Rho_Face(f,:)
  

       if ( .not. use_ff_support_operators ) then

          do n = 1,ndim
             Grad(n,:) = Gradient(n,f,:)
          enddo

          ! Recompute the face gradient at Dirichlet faces
          if (dirichlet_pressure) then
             
             if (use_ortho_face_gradient) then
                Grad_Dot_N = (BC_Prs(f,:)+ghn(f,:) - (Scalar_Cell_Center+ghc(f,:)))/Cell%Halfwidth(f)**2
             else
                Grad_Dot_N = (BC_Prs(f,:) - Scalar_Cell_Center)/Cell%Halfwidth(f)**2
             end if
             
             do n = 1,ndim
                where (Boundary_Flag(f,:) == 1) &
                     Grad(n,:) = Grad_Dot_N*(Cell(:)%face_centroid(n,f) - Cell(:)%Centroid(n))
             end do
          end if
          
          ! Correct the face velocities.
          do n = 1,ndim

             ! Apply the gradient correction.
             where (.not. Boundary_Flag(f,:) == 0 .and. &
                    .not. Boundary_Flag(f,:) == 2 .and. &
                    .not. Solid_Face(f,:)) &
                    Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) - &
                                dt_over_Rho*Grad(n,:)*Cell(:)%Face_Normal(n,f)

             ! Add the surface tension force 
             if (csf_normal) then
               Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + dt_over_Rho*Fcsf_new(n,f,:)*Cell(:)%Face_Normal(n,f)
             endif  

          end do

       else

          Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) - dt*Flux(f,:)
          if (csf_normal) then
             do n = 1,ndim
                Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + dt_over_Rho*Fcsf_new(n,f,:)*Cell%Face_Normal(n,f)
             enddo
          endif

       endif
          
    enddo FACE_LOOP

    ! Reset the pointer for pressures to BC_Zero for the pressure change solution
    BC_Prs => BC_Zero

    ! Clean up dynamic memory allocation
    DEALLOCATE(Scalar_Cell_Center)
    DEALLOCATE(Gradient)
    DEALLOCATE(Grad)
    if (.not. use_ff_support_operators) then
       DEALLOCATE(Grad_Dot_N)
    endif
    DEALLOCATE(dt_over_Rho)
    DEALLOCATE(Mask)

    return

  END SUBROUTINE VELOCITY_TO_FACES

  SUBROUTINE INTERPOLATE_VELOCITY_TO_FACES
    !=======================================================================
    ! Purpose(s):
    !
    !   Interpolate cell centered velocity vectors to normal velocity 
    !   components on faces
    !       Jim Sicilian,   May 2006
    !======================================================================= 
    use kind_module,            only: int_kind, real_kind
    use bc_module,              only: BC_Vel
    use constants_module,       only: one, zero
    use mesh_module,            only: Cell, Mesh
    use fluid_data_module,      only: Fluxing_Velocity, Face_Interpolation_Factor,  &
                                      fluidRho, IsPureImmobile, Solid_Face,         &
                                      Centered_GradP_Dynamic
    use gs_module,              only: EE_GATHER
    use parameter_module,       only: ncells, ndim, nfc
    use projection_data_module, only: Boundary_Flag
    use time_step_module,       only: dt
    use zone_module,            only: Zone
    use truchas_logging_services
    
    implicit none

    ! Argument List

    ! Local Variables
    real(real_kind),   dimension(:,:), allocatable, target :: Fluxing_Velocity_Ngbr
    real(real_kind),   dimension(:,:), pointer             :: fluidRho_Ngbr, Velocity_Component_Ngbr
    real(real_kind),   dimension(:),   allocatable         :: Velocity_Component, Int_factor 
    integer(int_kind)                                      :: f, n, c, status

    ! First, allocate temporary array space
    ALLOCATE (Fluxing_Velocity_Ngbr(nfc,ncells),      &
              Int_factor(ncells),                     &
              Velocity_Component(ncells),    STAT = status)
    if (status /= 0) call TLS_panic ('INTERPOLATE_VELOCITY_TO_FACES: allocation failed')

    Velocity_Component_Ngbr => Fluxing_Velocity_Ngbr
    Fluxing_Velocity = zero
    do n = 1, ndim
       Velocity_Component(:) = Zone(:)%Vc(n)
       where(fluidRho(:) > zero)
          Velocity_Component(:) = Velocity_Component(:) + dt*Centered_GradP_Dynamic(n,:)
       endwhere
       call EE_GATHER(Velocity_Component_Ngbr,Velocity_Component)
       do f = 1, nfc
          where(Boundary_Flag(f,:)<0)
             ! Interpolate at all interior faces
             Int_factor(:) = Face_Interpolation_Factor(f,:)
             Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + Cell(:)%Face_Normal(n,f) *   &
                                    ( (one-Int_factor(:))*Velocity_Component(:)  +       &
                                         Int_factor(:)   *Velocity_Component_Ngbr(f,:)      )
          elsewhere(Boundary_Flag(f,:)==0)
             ! Solid or Symmetry Boundaries
             Fluxing_Velocity(f,:) = zero
          elsewhere(Boundary_Flag(f,:)==1)
             ! Dirichlet Pressure Boundaries
             Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + Cell(:)%Face_Normal(n,f)*Velocity_Component(:)
          elsewhere(Boundary_Flag(f,:)==2)
             ! Dirichlet Velocity Boundaries
             Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + Cell(:)%Face_Normal(n,f)*BC_Vel(n,f,:)
          endwhere
       enddo
    enddo
    NULLIFY(Velocity_Component_Ngbr)

    ! Correct Special Cases
    fluidRho_Ngbr => Fluxing_Velocity_Ngbr
    call EE_Gather(fluidRho_Ngbr, fluidRho)
    do f = 1,nfc
       do c = 1, ncells
          if(Solid_Face(f,c)) then
             !  Zero Velocity on Solid Faces
             Fluxing_Velocity(f,c) = zero
          elseif(fluidRho_Ngbr(f,c)==zero .AND. Mesh(c)%Ngbr_cell(f) /= 0 ) then
             ! Use cell center velocity if adjacent cell is void
             Fluxing_Velocity(f,c) = zero
             do n = 1, ndim
                Fluxing_Velocity(f,c) = Fluxing_Velocity(f,c)+Cell(c)%Face_Normal(n,f)*  &
                                         (Zone(c)%Vc(n)+dt*Centered_GradP_Dynamic(n,c))
             enddo
          endif
       enddo
    enddo
    NULLIFY(fluidRho_Ngbr)

    ! make sure that void cell face velocities match the adjacent non-void cell
    call EE_GATHER (Fluxing_Velocity_Ngbr, -Fluxing_Velocity)
    do f = 1,nfc
       where (FluidRho == zero .AND. Mesh%Ngbr_cell(f) /= 0 .AND. .not.IsPureImmobile) &
          Fluxing_Velocity(f,:) = Fluxing_Velocity_Ngbr(f,:)
    end do

    DEALLOCATE(Fluxing_Velocity_Ngbr)
    DEALLOCATE(Velocity_Component)
    DEALLOCATE(Int_factor)
    return
  END SUBROUTINE INTERPOLATE_VELOCITY_TO_FACES

  SUBROUTINE MAC_CORRECTOR (Solution, Rho_Face, Fluxing_Velocity)
    !=======================================================================
    ! Purpose(s):
    !
    !   Apply a solenoidal correction, which is the gradient of the
    !   MAC solution field, to the face fluxing velocities.
    !=======================================================================
    use bc_module,              only: BC_Prs, BC_Vel
    use constants_module,       only: zero
    use do_interface,           only: DO_Specifier, do_init_ss, &
                                      do_gradient_face, DO_SOLVE_ORTHO, &
                                      DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: Solid_Face, FluidRho
    use kind_module,            only: int_kind, real_kind
    use mesh_module,            only: Cell
    use parameter_module,       only: ncells, ndim, nfc
    use projection_data_module, only: dirichlet_pressure, Boundary_Flag, &
                                      dt_gradP_over_Rho
    use time_step_module,       only: dt
    use support_operators ,     only: CALCULATE_FLUX
    use ff_discrete_ops_data,   only: use_ff_support_operators, FF_SO_Control_Data
    use truchas_logging_services
    
    implicit none

    ! Argument List
    real(real_kind), dimension(ncells),     intent(IN)    :: Solution
    real(real_kind), dimension(nfc,ncells), intent(IN)    :: Rho_Face
    real(real_kind), dimension(nfc,ncells), intent(INOUT) :: Fluxing_Velocity

    ! Local Variables
    type(DO_Specifier),pointer,save                :: MAC_Cor_SS =>NULL()
    integer                                        :: status
    integer(int_kind)                              :: f, n, MCSolveTech
    real(real_kind), dimension(:),     allocatable, save :: dt_over_Rho, Grad_Dot_N
    real(real_kind), dimension(:,:),   allocatable, save :: Grad
    real(real_kind), dimension(:,:,:), allocatable, save :: fgradx
 
    real(real_kind),   dimension(ncells)                 :: Scalar_Cell_Center
    real(real_kind), dimension(nfc,ncells)               :: SO_FluidRho_Face
    real(real_kind), dimension(nfc,ncells)               :: Flux

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (dt_over_Rho(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_CORRECTOR: dt_over_Rho(ncells) allocation failed')
    ALLOCATE (Grad_Dot_N(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_CORRECTOR: Grad_Dot_N(ncells) allocation failed')
    ALLOCATE (Grad(ndim,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_CORRECTOR: Grad(ndim,ncells) allocation failed')
    ALLOCATE (fgradx(ndim,nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_CORRECTOR: fgradx(ndim,nfc,ncells) allocation failed')

    fgradx = zero

    if ( use_ff_support_operators ) then
       
       ! set boundary conditions
       Scalar_Cell_Center = Solution
       
       SO_FluidRho_Face = zero
       DO f=1,nfc
          where (FluidRho(:) /= zero) SO_FluidRho_Face(f,:) = 1.0/FluidRho(:)
       ENDDO
          
       Flux = zero
       
       CALL CALCULATE_FLUX(PHI = Scalar_Cell_Center,      &
            FLUX = Flux,                   &
            COND_FACE = SO_FluidRho_Face,  &
            DIRBDY = BC_Prs,               &
            SO_CONTROL_DATA = FF_SO_Control_Data )
       
       DO f=1,nfc
          WHERE ( Cell%Face_Area(f) .ne. 0d0 ) 
             Flux(f,:) = Flux(f,:) / Cell%Face_Area(f)
          ENDWHERE
       END DO
       
    else

       if(.not.ASSOCIATED(MAC_Cor_SS))then
          MCSolveTech=DO_SOLVE_LU_LSLR; if(use_ortho_face_gradient)MCSolveTech=DO_SOLVE_ORTHO
          call do_init_ss(MAC_Cor_SS,SOLVETECH=MCSolveTech)
       endif
       
       call DO_GRADIENT_FACE(PHI=Solution, SOLVESPEC=MAC_Cor_SS, GRAD=fgradx)

    endif

    ! Loop over faces; compute the face potential gradient (dP_dx,dP_dy,dP_dz),
    ! then correct the face velocities with this gradient.

    FACE_LOOP: do f = 1,nfc

       if ( use_ff_support_operators ) then
          
          ! Correct the fluxing velocity with the normal component of the gradient.
          Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) &
               - dt*Flux(f,:)
             
       else

          do n = 1,ndim
             Grad(n,:) = fgradX(n,f,:)
          enddo
          
          ! Store the inverse face density.
          dt_over_Rho = zero
          where (Rho_Face(f,:) /= zero)  & 
               dt_over_Rho = dt/Rho_Face(f,:)
          
          ! Recompute the face gradient at Dirichlet faces
          if (dirichlet_pressure) then
             Grad_Dot_N = (BC_Prs(f,:) - Solution)/Cell%Halfwidth(f)**2
             do n = 1,ndim
                where (Boundary_Flag(f,:) == 1) &
                     Grad(n,:) = Grad_Dot_N*(cell(:)%face_centroid(n,f) - Cell(:)%Centroid(n))
             end do
          end if
          
          ! Kill the gradient at solid faces.
          do n = 1,ndim
             where (Solid_Face(f,:)) Grad(n,:) = zero
          end do
          
          ! Get the normal component of the gradient.
          Grad_Dot_N = zero
          do n = 1,ndim
             Grad_Dot_N = Grad_Dot_N + Grad(n,:)*Cell%Face_Normal(n,f)
          end do
          
          ! Correct the face velocities.
          do n = 1,ndim
             
             ! Kill the normal component of the gradient where appropriate.
             where (Boundary_Flag(f,:) == 0 .or. Boundary_Flag(f,:) == 2) &
                  Grad(n,:) = Grad(n,:) - Grad_Dot_N*Cell%Face_Normal(n,f)
             
            ! Evaluate dt_gradP_over_Rho to be used later to calculate the cell-centred gradient.
             dt_gradP_over_Rho(n,f,:) = dt_over_Rho*Grad(n,:)
             
          end do
          
          ! Correct the fluxing velocity with the normal component of the gradient.
          Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) - dt_over_Rho*Grad_Dot_N
       endif

       ! Now Apply BCs to the Fluxing_Velocity;

       ! Free-slip, Fluxing_Velocity = zero.
       where (Boundary_Flag(f,:) == 0 .or. Boundary_Flag(f,:) == 2) &
          Fluxing_Velocity(f,:) = zero

       ! Dirichlet velocity BCs:  Fluxing_Velocity = BC_Vel
       do n = 1,ndim
          where (Boundary_Flag(f,:) == 2) &
             Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + &
                                     BC_Vel(n,f,:)*Cell%Face_Normal(n,f)
       end do

       ! Zero out Fluxing_Velocity's at solid faces.
       where (Solid_Face(f,:)) Fluxing_Velocity(f,:) = zero

    end do FACE_LOOP

    DEALLOCATE (dt_over_Rho)
    DEALLOCATE (Grad_Dot_N)
    DEALLOCATE (Grad)
    DEALLOCATE (fgradx)

    return

  END SUBROUTINE MAC_CORRECTOR



  SUBROUTINE MAC_RHS (Fluxing_Velocity, RHS)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the RHS for the MAC projection, which is the
    !   cell-centered divergence of the face fluxing velocity.
    !=======================================================================
    use constants_module,       only: zero
    use fluid_data_module,      only: Solid_Face
    use kind_module,            only: int_kind, real_kind
    use mesh_module,            only: Cell
    use parameter_module,       only: ncells, nfc
    use projection_data_module,      only: DVol_by_Dt_over_Vol

    implicit none

    ! Argument List
    real(real_kind), dimension(nfc,ncells), intent(IN)  :: Fluxing_Velocity
    real(real_kind), dimension(ncells),     intent(OUT) :: RHS

    ! Local Variables
    integer(int_kind) :: f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Compute the cell-centered divergence of the fluxing velocity as a
    ! discrete sum over faces (Green-Gauss).  Recall that the fluxing
    ! velocity has already been projected along the face normal. This will
    ! be the RHS for MAC projection.

    RHS = zero
    do f = 1, nfc
       where(.not.solid_face(f,:))
           RHS = RHS + Fluxing_Velocity(f,:)*Cell%Face_Area(f)
       endwhere
    end do

    !  Add in the volume change rate term due to the Enthalpy Solution
    RHS = RHS -  DVol_by_Dt_over_Vol*Cell%Volume

!!$    Don't add this in because we are solving for delta-P
!!$    RHS = RHS - Vol_over_RhoCsqDt*Zone%P


    return

  END SUBROUTINE MAC_RHS

  SUBROUTINE NONORTHO_PROJECTION (RHS, Rho_Face, Solution, timer_flag)
    !======================================================================
    ! Purpose(s):
    !
    !   Solve the equation del*Q = S with a general nonorthogonal operator
    !   for the del*Q term.  Here Q is a "pressure flux", Q = 1/rho Grad(P),
    !   where rho is the density and P is the pressure. S is a source term:
    !   S = D/dt for the regular pressure projection and S = D for the MAC
    !   projection, where D is a velocity divergence.
    !======================================================================
    use ArrayAllocate_Module,   only: ARRAYCREATE, ARRAYDESTROY
    use bc_module,              only: BC_Prs
    use constants_module,       only: zero, one
    use cutoffs_module,         only: alittle
    !  FluidRho is needed for scaling of A_ortho and RHS
    use fluid_data_module,      only: FluidRho, isPureImmobile, MinFaceFraction, MinFluidRho
    use kind_module,            only: int_kind, real_kind
    use linear_solution,        only: LINEAR_SOLVER, Ubik_user
    use lnorm_module,           only: L1NORM
    use mesh_module,            only: Cell
    use parameter_module,       only: ncells, ndim, nfc
    use preconditioners,        only: P, PRECONDITION
    use projection_data_module, only: A_Ortho, Boundary_Flag, &
                                      Coeff, UBIK_PRESSURE
    use y_eq_Ax_prs,            only: Y_EQ_AX_PRESSURE
    use time_step_module,       only: dt
    use timing_tree
    use ff_discrete_ops_data,   only: use_ff_support_operators
    use UbikSolve_module

    use fischer_module

    implicit none

    ! Argument List
    character(len=*),                       intent(IN)    :: timer_flag
    real(real_kind), dimension(nfc,ncells), intent(IN)    :: Rho_Face
    real(real_kind), dimension(ncells),     intent(INOUT) :: RHS
    real(real_kind), dimension(ncells),     intent(INOUT) :: Solution


    ! Local Variables
    integer(int_kind) :: f, n, nc
    real(real_kind)   :: X_Dot_N
    real(real_kind)   :: RHSnorm

    ! Convergence criteria used with new PPE
    real(real_kind)   :: original_criterion, convergence_criterion

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! If preconditioning is on, get the orthogonal coefficients.
    if (Ubik_user(UBIK_PRESSURE)%precond /= 0) then

       ! Allocate the orthogonal coeffient array.
       call ARRAYCREATE (A_Ortho, 0, nfc, 1, ncells, 'A_Ortho(0:nfc,ncells)')

       ! Point the preconditioning matrix to A_Ortho
       P => A_Ortho

       ! Set up orthogonal mesh 7-point stencil
       call ORTHO_STENCIL (RHO_FACE = Rho_Face, MATRIX = A_Ortho)

       ! Scale A_Ortho to be consistent with the RHS and mat-vec.
       ! Note 0 index for lower bound of A_Ortho
       do n=1,ncells
          if(FluidRho(n) == zero .or. isPureImmobile(n)) then
              A_Ortho(0,n) = one / (MinFaceFraction * MinFluidRho * Cell(n)%Volume**0.66666)
          else
              do f=0,nfc
                  A_Ortho(f,n) = A_Ortho(f,n)/Cell(n)%Volume
              enddo
          endif
       enddo

    end if

    if ( .not. use_ff_support_operators ) then

       ! Allocate coefficient array.
       call ARRAYCREATE (Coeff, 1, nfc, 1, ncells, 'Coeff(nfc,ncells)')
       Coeff = zero
       
       ! Set Coeff = Area_Face/Rho_Face. If the face has Dirichlet BCs,
       ! then Coeff = -Area_Face*(X*n)/(Rho_Face*L**2), where X is a vector
       ! from the cell to face centroid, n is the unit normal, and L is
       ! the cell halfwidth.
       
       do nc = 1,ncells

          if(isPureImmobile(nc) .or. FluidRho(nc) == zero) cycle
          do f = 1,nfc
             
             ! Coeff = A_f/Rho_f.
             if (Rho_Face(f,nc) == zero) then
                Coeff(f,nc) = zero
                CYCLE
             else
                Coeff(f,nc) = Cell(nc)%Face_Area(f) / Rho_Face(f,nc)
             endif
             
             ! Dirichlet corrections.
             if (Boundary_Flag(f,nc) == 1) then
                X_Dot_N = zero
                do n = 1,ndim
                   X_Dot_N = X_Dot_N + (Cell(nc)%face_centroid(n,f) &
                        - Cell(nc)%Centroid(n))*Cell(nc)%Face_Normal(n,f)
                end do
                X_Dot_N = X_Dot_N/Cell(nc)%Halfwidth(f)**2
                Coeff(f,nc) = -Coeff(f,nc)*X_Dot_N
                RHS(nc) = RHS(nc) + BC_Prs(f,nc)*dt*Coeff(f,nc)
             end if
             
          end do
          
       end do
    end if

    ! Likely not the best place to put this, but the RHS must be divided by
    ! 'dt' so that we're solving for the pressure increment, and not pressure*dt.
    ! Scaling of the RHS here by dt and the cell volume is in conjunction with the
    ! scaled convergence criteria below.
    where(FluidRho(:) == zero .or. isPureImmobile)
         RHS(:) = RHS(:) / (MinFaceFraction * MinFluidRho * Cell(:)%Volume**0.66666)
    elsewhere
         RHS(:) = RHS(:)/(dt*Cell(:)%Volume)
    endwhere

    ! Calculate scaling factor for L2 norm of the residual
    original_criterion = Ubik_eps(Ubik_user(UBIK_PRESSURE)%Control)
    convergence_criterion = original_criterion/(dt*dt)
    call Ubik_set_eps(Ubik_user(UBIK_PRESSURE)%Control, convergence_criterion)

    ! Solve for the pressure correction
    Solution = zero

    ! prevent use of an initial guess of zero if the right-hand side is also zero
    RHSnorm = L1NORM(RHS)
    if (RHSNorm < alittle) then
          Solution = 1.0
    else
          call Fischer_Initial_Guess( Solution, RHS )
    end if

    ! Solve the linear system.
    call start_timer(timer_flag)
    
    call LINEAR_SOLVER(Solution, RHS, Ubik_user(UBIK_PRESSURE), &
                       Y_EQ_AX_PRESSURE, PRECONDITION)

    call stop_timer(timer_flag)

    call Fischer_Update_Space( Solution, Y_EQ_AX_PRESSURE, RHS )

    ! Store off the original convergence criteria for next time-step
    call Ubik_set_eps(Ubik_user(UBIK_PRESSURE)%Control, original_criterion)

    ! Deallocate working arrays.
    if ( .not. use_ff_support_operators ) then
       call ARRAYDESTROY (Coeff, 'Coeff(nfc,ncells)')
    end if

    if (Ubik_user(UBIK_PRESSURE)%precond /= 0) then
       NULLIFY (P)
       call ARRAYDESTROY (A_Ortho, 'A_Ortho(0:nfc,ncells)')
    end if

    return

  END SUBROUTINE NONORTHO_PROJECTION

  SUBROUTINE ORTHO_STENCIL (Rho_Face, RHS, Matrix)
    !======================================================================
    ! Purpose(s):
    !
    !   Discretize the equation del*Q = S with an orthogonal 7-point
    !   operator for the del*Q term.  Here Q is a "pressure flux",
    !   Q = 1/rho Grad(P), where rho is the density and P is the pressure.
    !   S is a source term: S = D/dt for the regular pressure projection
    !   and S = D for the MAC projection, where D is a velocity divergence.
    !
    !   Integrating over a control volume:
    !    ---
    !    \               Af       dX
    !    /   (Pnb - Pp)*----- * ------     = S*V
    !    ---            L*rhof     L
    !   faces
    !
    ! where V is the cell volume, Af the face area, Pnb is the neighbor
    ! cell pressure, Pp is the cell pressure, L the centroid-to-centroid
    ! distance, and rhof the face density.
    !
    ! This is the same operator as the "pseudo-ortho" approximation used in
    ! the full matvec.
    !
    !    This latest implementation due to John Turner and Jim Sicilian
    !    CCS-2   October, 2003
    !
    !======================================================================
    use bc_module,              only: BC_Prs, BC, DIRICHLET, Prs
    use constants_module,       only: zero, one
    use cutoffs_module,         only: alittle
    use fluid_data_module,      only: FluidRho, Solid_Face, isPureImmobile
    use gs_module,              only: EE_GATHER
    use kind_module,            only: int_kind, real_kind
    use mesh_module,            only: Cell, Mesh
    use parameter_module,       only: ncells, nfc, ndim
    use projection_data_module, only: dirichlet_pressure, Vol_over_RhoCsqDt
    use time_step_module,       only: dt

    implicit none

    ! Argument List
    real(real_kind), dimension(nfc,ncells),             intent(IN)    :: Rho_Face
    real(real_kind), dimension(0:nfc,ncells),           intent(INOUT) :: Matrix
    real(real_kind), dimension(ncells),       optional, intent(INOUT) :: RHS

    ! Local Variables
    integer(int_kind)                        :: j, f, n
    real(real_kind)                          :: distance, dX_tmp
    real(real_kind),   dimension(nfc,ncells) :: Face_Coeff
    real(real_kind),   dimension(ndim,nfc,ncells) :: dx_scaled, X_e
    real(real_kind),   dimension(ndim)            :: dx

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Face_Coeff = zero

      ! Quantities for pseudo-ortho operator
      do n = 1,ndim
         call EE_GATHER (DEST=X_e(n,:,:), SRC=Cell(:)%Centroid(n))
         do f = 1,nfc
            where (Mesh%Ngbr_cell(f) == 0 ) X_e(n,f,:)=Cell(:)%Face_Centroid(n,f)
         enddo
      end do

      CELL_LOOP: do j=1,ncells
      FACE_LOOP: do f = 1,nfc
         ! Compute Deltas (Cell to Face Neighbors)
         ! Physical Coordinate Deltas
         dX       = zero
         Distance = zero
         do n = 1,ndim
            dX(n)  = X_e(n,f,j) - Cell(j)%Centroid(n)
            Distance = Distance + dX(n)**2
         end do
         do n = 1,ndim
           dX_tmp = dX(n)/(Distance + alittle)
           dX_scaled(n,f,j) = dX_tmp
         end do
        ! Define face coefficients
        if (Rho_Face(f,j) == zero) then
           Face_Coeff(f,j) = zero
        else
           do n= 1,ndim
              Face_Coeff(f,j) = Face_Coeff(f,j) + Cell(j)%Face_Normal(n,f)*dX_scaled(n,f,j)
           enddo
           Face_Coeff(f,j) = Face_Coeff(f,j) * Cell(j)%Face_Area(f)/Rho_Face(f,j)
        end if
      end do FACE_LOOP
      end do CELL_LOOP


    ! Now load the coefficients at each face; check for BCs
    Matrix = zero
    Matrix(0,:)= -Vol_over_RhoCsqDt / dt
    FACE_INTEGRAL: do f = 1,nfc

       ! Accumulate coefficients on internal faces that are NOT
       ! across from an internal flow obstacle.
       where (Mesh%Ngbr_cell(f) /= 0 .AND. .NOT. Solid_Face(f,:))
          Matrix(0,:) = Matrix(0,:) - Face_Coeff(f,:)
          Matrix(f,:) = Face_Coeff(f,:)
       end where

       ! Check for Dirichlet boundaries; Correct the coefficient.
       DIRICHLET_COEF: if (dirichlet_pressure) then
          where (DIRICHLET (BC%Flag, Prs%Face_bit(f)))
             Matrix(0,:) = Matrix(0,:) - Face_Coeff(f,:)
          end where
       end if DIRICHLET_COEF

       ! Check for Dirichlet boundaries; Correct RHS
       DIRICHLET_RHS: if (PRESENT(RHS) .and. dirichlet_pressure) then
          where (DIRICHLET (BC%Flag, Prs%Face_bit(f)))
             RHS = RHS - Face_Coeff(f,:) * BC_Prs(f,:)
          end where
       end if DIRICHLET_RHS

    end do FACE_INTEGRAL

    ! Finally, take care of void and immobile cells
    do n = 1, ncells
       if(FluidRho(n) == zero .or. isPureImmobile(n)) then
          Matrix(0,n) = one
          do f=1,nfc
             Matrix(f,n) = zero
          enddo
       endif
    enddo

    return

  END SUBROUTINE ORTHO_STENCIL

  SUBROUTINE PROJECTION_CORRECTOR ()
    !=======================================================================
    ! Purpose(s):
    !
    !   The last step in the PROJECTION call, this routine calculates a 
    !   time 'n+1' cell-centered pressure gradient, and applies it to the
    !   cell-centered velocities.
    !=======================================================================
    use bc_module,              only: BC_Prs
    use bc_data_module,         only: BC_Zero
    use bc_operations,          only: Pressure_BC
    use constants_module,       only: zero, one
    use do_interface,           only: DO_Specifier, do_init_ss, &
                                      DO_SOLVE_ORTHO, DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: isPureImmobile, fluidRho, &
                                      Centered_GradP_Dynamic,   &
                                      Momentum_by_Volume,       &
                                      fluidVof
    use kind_module,            only: int_kind, real_kind
    use parameter_module,       only: ndim, ncells
    use time_step_module,       only: dt
    use timing_tree
    use zone_module,            only: Zone
    use fluid_utilities_module, only: CC_GRADP_DYNAMIC

    implicit none

    ! Argument List

    ! Local Variables

    integer(int_kind)                                :: n
    integer(int_kind)                                :: PCSolveTech
    real(KIND = real_kind),  dimension(ncells)       :: Kappa

    type(DO_Specifier),pointer,save                  :: PCorrector_SS =>NULL()

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the projection_corrector timer.
    call start_timer("timer_projection_corrector")

    if(.not.ASSOCIATED(PCorrector_SS))then
      PCSolveTech=DO_SOLVE_LU_LSLR; if(use_ortho_face_gradient)PCSolveTech=DO_SOLVE_ORTHO
      call do_init_ss(PCorrector_SS,SOLVETECH=PCSolveTech,BC_SPEC=Pressure_BC)
    endif

    ! Evaluate and store the volume change rate due to the Enthalpy solve
    ! NNC, 2 May 2012.  It looks like TI_GET_DV_BY_VO returns the relative
    ! density change (not volume change!) that occured from the beginning
    ! of the enthalpy solve to the end.  Any density change due to advection
    ! is entirely different.  The only density change that could occur over
    ! the heat transfer step is due to temperature-dependent densities and/or
    ! phase transformations with densities differing between phases.  Since
    ! we allow neither (and never really have done so correctly), KAPPA = 1.
    !if(heat_conduction) then
    !    Kappa = one - TI_Get_dV_by_Vo()
    !else
        Kappa = one
    !endif

    ! Remove the old time dynamic pressure gradient acceleration from the 
    ! cell centered velocity
    do n = 1,ndim
       Zone%Vc(n) = Zone%Vc(n) + dt * Centered_GradP_Dynamic(n,:)
       Momentum_by_Volume(n,:) = Momentum_by_Volume(n,:) + &
                      Centered_GradP_Dynamic(n,:)*dt*fluidRho*fluidVof/Kappa
       ! This should probably be here to ensure the zero's are consistent w. the 
       ! loop following cc_gradp_dynamic() below.
       where (FluidRho == zero .or. isPureImmobile) 
          Zone%Vc(n) = zero
          Momentum_by_Volume(n,:) = zero
       endwhere
    end do

    call CC_GRADP_DYNAMIC()

    ! Apply the cell-centered gradient to Zone%Vc, and zero out values in 
    ! void and solid cells.    
    do n = 1,ndim
       Zone%Vc(n) = Zone%Vc(n) - Centered_GradP_Dynamic(n,:) * dt
       Momentum_by_Volume(n,:) = Momentum_by_Volume(n,:) - &
                     Centered_GradP_Dynamic(n,:) * dt * fluidRho*fluidVof/Kappa
       where (FluidRho == zero .or. isPureImmobile) 
                Zone%Vc(n) = zero
                Momentum_by_Volume(n,:) = zero
       endwhere
    end do

    ! Reset the pointer for pressures to BC_Zero for the pressure change solution
    BC_Prs => BC_Zero

    ! stop the projection_corrector timer.
    call stop_timer("timer_projection_corrector")

    return

  END SUBROUTINE PROJECTION_CORRECTOR


  SUBROUTINE VF_DIVERGENCE_NORMS (Fluxing_Velocity)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the L1, L2, and Linf norms of the divergence of the
    !   face-centered fluxing velocity after a MAC projection and
    !   correction has been applied. Since the MAC projection is an
    !   exact projection, these norms should be of order epsilon,
    !   where epsilon is the convergence tolerance of the linear MAC
    !   projection system.
    !=======================================================================
    use constants_module,  only: zero
    use fluid_data_module, only: fluidRho
    use fluid_type_module, only: Div_c
    use kind_module,       only: int_kind
    use lnorm_module,      only: L1NORM, L2NORM, LINORM
    use mesh_module,       only: Cell
    use parameter_module,  only: ncells, nfc
    use pgslib_module,     only: PGSLIB_GLOBAL_MAXLOC
    use time_step_module,  only: dt

    implicit none

    ! Argument List
    real(real_kind), dimension(nfc,ncells), intent(IN) :: Fluxing_Velocity

    ! Local Variables
    integer(int_kind)                    :: f
    integer(int_kind), dimension(1)      :: MaxLoc_L
    real(real_kind),   dimension(ncells) :: Velocity_Divergence

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Compute the velocity divergence.
    Velocity_Divergence = zero
    do f = 1, nfc
       Velocity_Divergence = Velocity_Divergence + &
                             Fluxing_Velocity(f,:)*Cell%Face_Area(f)
    end do

    Velocity_Divergence = (Velocity_Divergence/Cell%Volume)*dt

    ! Zero out divergences in void cells
    where (FluidRho == zero) Velocity_Divergence = zero

    !save as a residual

    ! Now get the norms; store them in Div.
    Div_c%V_f%L1   = L1NORM(Velocity_Divergence)
    Div_c%V_f%L2   = L2NORM(Velocity_Divergence)
    Div_c%V_f%Linf = LINORM(Velocity_Divergence)

    ! Get the location of the Linf norm.
    MaxLoc_L = PGSLib_Global_MAXLOC(ABS(Velocity_Divergence))
    Div_c%V_f%Linf_Location = MaxLoc_L(1)

    return

  END SUBROUTINE VF_DIVERGENCE_NORMS

END MODULE PROJECTION_MODULE
