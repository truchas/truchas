!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  public :: PROJECTION, ORTHO_STENCIL

  ! Power for averaging of sound speed by material
  integer, parameter :: navg = 1

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
    use cutoffs_module,         only: alittle
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: fluidRho, Solid_Face, IsPureImmobile, &
        Rho_Face, Rho_Face_n, Fluxing_Velocity
    use fluid_utilities_module, only: FLUIDDENSITYFACE
    use legacy_mesh_api,        only: ncells, ndim, nfc, Mesh, Cell, EE_GATHER
    use projection_data_module, only: dt_gradP_over_Rho, Vol_over_RhoCsqDt, &
        ghc, ghn, dtRhoG_over_Rho, &
        Fcsf_new,dtCsf_over_Rho
    use surface_tension_module, only: CSF_FACE, csf_normal
    use time_step_module,       only: dt
    use fischer_module

    ! Local Variables
    integer :: status
    integer :: n, f
    real(r8), dimension(:,:), allocatable, save :: Scalar_Cell_Face
    real(r8), dimension(:,:), allocatable, save :: Scalar_Cell_Ngbr
    logical, save :: first_time = .true.
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the projection timer.
    call start_timer("timer_projection")

    ALLOCATE (Scalar_Cell_Face(nfc,ncells),       &
        Scalar_Cell_Ngbr(nfc,ncells),       &
        dt_gradP_over_Rho(ndim,nfc,ncells), &
        Vol_over_RhoCsqDt(ncells),          &
        Rho_Face(nfc,ncells),               STAT = status)
    if (status /= 0) call TLS_panic ('PROJECTION: allocation failed')

    if (.not. use_ortho_face_gradient) then
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
    write(*, '("Pre-VELOCITY_TO_FACES FV[",i4,"]: ", 6es20.12)') 771, Fluxing_Velocity(:,771)
    call VELOCITY_TO_FACES()
    write(*, '("Post-VELOCITY_TO_FACES FV[",i4,"]: ", 6es20.12)') 771, Fluxing_Velocity(:,771)

    if (csf_normal) then
      do f=1,nfc
        do n=1,ndim
          ! particular case when void cells. density = zero
          where (Rho_Face(f,:) <= alittle)
            dtCsf_over_Rho(n,f,:)=0
          elsewhere
            dtCsf_over_Rho(n,f,:)=dt*Fcsf_new(n,f,:)/Rho_Face(f,:)
          endwhere
        end do
      end do
    endif

    if (.not. use_ortho_face_gradient) then
      ! Add gravitational acceleration to face velocities and compute Fluxing_Velocity
      call ADD_FACE_BODY_FORCE(dt, dtRhoG_over_Rho)
      do f = 1, nfc
        do n = 1, ndim
          Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + &
              dtRhoG_over_Rho(n,f,:)*Cell(:)%Face_Normal(n,f)
        end do
      end do
    end if
#ifndef NDEBUG
    write(*,*) "<< Pre Correction Fluxing Velocities"
    do n = 1, ncells
      write(*,'("Cell [",i3,"]:",6es15.5)') n, Fluxing_Velocity(:,n)
    end do
#endif
    ! Correct the Fluxing_Velocity's with a MAC projection to make them solenoidal.
    call MAC_PROJECTION(Fluxing_Velocity)

    ! Void cells don't have physical (nor solenoidal) Fluxing_Velocity's
    ! Set Fluxing_Velocity on faces that adjoin 'real' fluid cells to the
    ! corresponding velocity on the other side of the face (corrected for direction)
    call EE_GATHER(Scalar_Cell_Face, -Fluxing_Velocity)
    do f = 1,nfc
      where (FluidRho == 0 .AND. &
          Mesh%Ngbr_cell(f) /= 0 .AND. .not.IsPureImmobile) &
          Fluxing_Velocity(f,:) = Scalar_Cell_Face(f,:)
    end do

    ! Zero out Fluxing_Velocity's at solid faces.
    where (Solid_Face) Fluxing_Velocity = 0

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

    if (.not. use_ortho_face_gradient) then
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
    use bc_module,              only: bndry_vel !BC_Vel
    use fluid_data_module,      only: fluidRho, Solid_Face, &
        void_pressure, IsPureImmobile,  &
        Rho_Face, IsImmobile
    use matl_module,            only: Matl
    use linear_solution,        only: Ubik_user
    use legacy_mesh_api,        only: ncells, ndim, nfc, ncells_tot, Cell
    use parameter_module,       only: nmat, mat_slot
    use pgslib_module,          only: PGSLIB_GLOBAL_SUM, PGSLIB_GLOBAL_ANY
    use projection_data_module, only: mac_projection_iterations,    &
        mac_projection_precond_iter,  &
        Boundary_Flag, UBIK_PRESSURE,   &
        Face_Density,                 &
        dirichlet_pressure,           &
        Vol_over_RhoCsqDt
    use property_data_module,   only: Sound_Speed
    use time_step_module,       only: t, dt
    use zone_module,            only: Zone
    use UbikSolve_module

    ! Argument List
    real(r8), dimension(nfc,ncells), intent(INOUT) :: Fluxing_Velocity

    ! Local Variables
    !    integer :: status
    logical :: Void_Cell_Found
    integer :: f, i, m, n, s, status
    real(r8), dimension(:), allocatable :: RHS, Solution, boundary_fv

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the MAC projection timer.
    call start_timer("timer_projection_mac")
    ALLOCATE (RHS(ncells),      &
        Solution(ncells), &
        boundary_fv(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_PROJECTION: allocation failed')

    ! Zero out Fluxing_Velocity's at solid faces.
    where (Solid_Face) Fluxing_Velocity = 0

    ! Compute cell compressibility expression
    Vol_over_RhoCsqDt = 0
    ! Temporary use of RHS
    RHS = 0
    ! Temporarily store averaged sound speed in Vol_over_RhoCsqDt
    do m = 1,nmat
      if(.not.isImmobile(m)) then
        do s = 1, mat_slot
          where (Matl(s)%Cell%Id == m) &
              RHS = RHS +  Matl(s)%Cell%Vof**navg
        end do
      endif
      ! calculate the contribution of material m to the reciprocal sound speed square
      if(Sound_Speed(m) > 0) then
        do s = 1, mat_slot
          where (Matl(s)%Cell%Id == m) &
              Vol_over_RhoCsqDt = Vol_over_RhoCsqDt + Matl(s)%Cell%Vof**navg/sound_speed(m)**2
        end do
      endif
    enddo
    ! Normalize by the total fluid fraction, Multiply by the Cell Volume and divide by dt
    ! (This puts the array in the most efficient form for Y_EQ_AX_PRS)
    where( RHS > 0 )  Vol_over_RhoCsqDt = Cell%Volume*Vol_over_RhoCsqDt / (dt*RHS)
    RHS = 0
    ! Now divide by the fluid density
    where(FluidRho > 0 ) Vol_over_RhoCsqDt = Vol_over_RhoCsqDt / FluidRho

    Face_Density = Rho_Face

    ! Get the RHS for this system
    call MAC_RHS (Fluxing_Velocity, RHS)

    do f = 1,nfc

      where (Boundary_Flag(f,:) == 0) &
          RHS = RHS - Fluxing_velocity(f,:)*Cell%Face_Area(f)

      !! NNC, Jan 2014.  Time-dependent dirichlet velocity.
      !ORIG: boundary_fv = 0
      !ORIG: do n = 1,ndim
      !ORIG:    where (Boundary_Flag(f,:) == 2) &
      !ORIG:       boundary_fv = boundary_fv + BC_Vel(n,f,:)*Cell%Face_Normal(n,f)
      !ORIG: end do
      do i = 1, ncells
        if (boundary_flag(f,i) == 2) then
          boundary_fv(i) = dot_product(bndry_vel%get(f,i,t), Cell(i)%Face_Normal(:,f))
        else
          boundary_fv(i) = 0.0_r8
        end if
      end do
      where (Boundary_Flag(f,:) == 2) &
          RHS = RHS + (boundary_fv - Fluxing_velocity(f,:))*Cell%Face_Area(f)
    end do

    ! set RHS to zero in void cells because we are solving for the change in pressure
    ! (RHS divided by dt in NONORTHO_PROJECTION)
    ! Scaling of the RHS  requires multiplying by dt here for
    ! global scaling of RHS by dt/Volume.
    where (FluidRho(:) == 0 .or. isPureImmobile(:)) RHS(:) = (void_pressure-Zone(:)%P)
#ifndef NDEBUG
    write(*,*) "<<< Baseline RHS"
    do f = 1, ncells
      write(*,'("Rhs[",i3, "]: ", es15.5)') f, RHS(f)
    end do
#endif
    call NONORTHO_PROJECTION (RHS, Face_Density, Solution, "timer_projection_solver")

    ! Store the number of iterations.
    mac_projection_iterations   = Ubik_iter(Ubik_user(UBIK_PRESSURE)%control)
    mac_projection_precond_iter = Ubik_user(UBIK_PRESSURE)%precond_iter

    ! Apply the solenoidal correction to the face-centered velocities.
    call MAC_CORRECTOR (Solution, Face_Density, Fluxing_Velocity)
#ifndef NDEBUG
    write(*,*) "<<< DP"
    do f = 1, ncells
      write(*,'("DP[",i3, "]: ", es15.5)') f, Solution(f)
    end do
#endif
    ! Computed divergence norms based on the correct fluxing velocity.
    call VF_DIVERGENCE_NORMS (Fluxing_Velocity)

    ! If no dirichlet pressures specified, normalize the solution by
    ! subtracting off the global mean.
    Void_Cell_Found = .false.
    do i= 1,ncells
      if(FluidRho(i) == 0 .and. .not.IsPureImmobile(i)) then
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
    use do_interface,           only: DO_Specifier, do_init_ss, &
                                      do_gradient_face, DO_SOLVE_ORTHO, &
                                      DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: Rho_Face, Solid_Face, Fluxing_Velocity
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell
    use projection_data_module, only: dirichlet_pressure, Boundary_Flag, ghc, ghn, &
                                      Fcsf_new
!-mf jan04
    use time_step_module,       only: dt
    use zone_module,            only: Zone
    use discrete_op_module,     only: DYNAMIC_PRESSURE_FACE_GRADIENT
    use surface_tension_module, only: csf_normal

    ! local variables
    integer :: n, f, status
    logical,  dimension(:),   allocatable :: Mask
    real(r8), dimension(:),   allocatable :: Scalar_Cell_Center
    real(r8), dimension(:),   allocatable :: dt_over_Rho, Grad_Dot_N
    real(r8), dimension(:,:), allocatable :: Grad
    real(r8), dimension(:,:,:), allocatable :: Gradient
    integer :: PSolveTech
    type(DO_Specifier),pointer,save :: Projection_SS2 =>NULL()

   ! First, allocate temporary array space
   ALLOCATE (Mask(ncells),                 &
             Scalar_Cell_Center(ncells),   &
             Gradient(ndim,nfc,ncells),    &
             Grad(ndim,ncells),            &
             dt_over_Rho(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VELOCITY_TO_FACES: allocation failed')

    ALLOCATE (Grad_Dot_N(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('VELOCITY_TO_FACES: Grad_Dot_N(ncells) allocation failed')

    ! Begin by extrapolating the time '*' cell-centered velocity to the faces
    write(*, '("Pre-INTERPOLATE FV[",i4,"]: ", 6es20.12)') 771, Fluxing_Velocity(:,771)
    call INTERPOLATE_VELOCITY_TO_FACES
    write(*, '("Post-INTERPOLATE FV[",i4,"]: ", 6es20.12)') 771, Fluxing_Velocity(:,771)
    ! Use the real Pressures at Dirichlet boundaries in this pressure gradient calculation
    BC_Prs => BC_Pressure

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

    FACE_LOOP : do f = 1,nfc

       ! Store the inverse face density. (needed for the surface tension model)
       dt_over_Rho = 0
       where (Rho_Face(f,:) /= 0) &
            dt_over_Rho = dt/Rho_Face(f,:)

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

    enddo FACE_LOOP

    ! Reset the pointer for pressures to BC_Zero for the pressure change solution
    BC_Prs => BC_Zero

    ! Clean up dynamic memory allocation
    DEALLOCATE(Scalar_Cell_Center)
    DEALLOCATE(Gradient)
    DEALLOCATE(Grad)
    DEALLOCATE(Grad_Dot_N)
    DEALLOCATE(dt_over_Rho)
    DEALLOCATE(Mask)

  END SUBROUTINE VELOCITY_TO_FACES

  SUBROUTINE INTERPOLATE_VELOCITY_TO_FACES
    !=======================================================================
    ! Purpose(s):
    !
    !   Interpolate cell centered velocity vectors to normal velocity
    !   components on faces
    !       Jim Sicilian,   May 2006
    !=======================================================================
    use bc_module,              only: bndry_vel !BC_Vel
    use fluid_data_module,      only: Fluxing_Velocity, Face_Interpolation_Factor,  &
                                      fluidRho, IsPureImmobile, Solid_Face,         &
                                      Centered_GradP_Dynamic
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell, Mesh, EE_GATHER
    use projection_data_module, only: Boundary_Flag
    use time_step_module,       only: t, dt
    use zone_module,            only: Zone

    ! Local Variables
    real(r8), dimension(:,:), allocatable, target :: Fluxing_Velocity_Ngbr
    real(r8), dimension(:,:), pointer             :: fluidRho_Ngbr, Velocity_Component_Ngbr
    real(r8), dimension(:),   allocatable         :: Velocity_Component, Int_factor
    integer :: f, n, c, status

    ! First, allocate temporary array space
    ALLOCATE (Fluxing_Velocity_Ngbr(nfc,ncells),      &
              Int_factor(ncells),                     &
              Velocity_Component(ncells),    STAT = status)
    if (status /= 0) call TLS_panic ('INTERPOLATE_VELOCITY_TO_FACES: allocation failed')

    Velocity_Component_Ngbr => Fluxing_Velocity_Ngbr
    Fluxing_Velocity = 0
    write(*, '("  Face_Interpolation_Factor[",i4,"]: ", 6es20.12)') 771, Face_Interpolation_Factor(:,771)
    do n = 1, ndim
       Velocity_Component(:) = Zone(:)%Vc(n)
       where(fluidRho(:) > 0)
          Velocity_Component(:) = Velocity_Component(:) + dt*Centered_GradP_Dynamic(n,:)
       endwhere
       call EE_GATHER(Velocity_Component_Ngbr,Velocity_Component)
       write(*, '("    VC[",i4,"]: ", es20.12)') 771, Velocity_Component(771)
!!$       write(*, '("    VC[",i4,"] y-/+: ", 2es20.12)') 771, Velocity_Component(21), Velocity_Component(23)
!!$       write(*, '("    VC[",i4,"] y-/+: ", 2es20.12)') 771, Velocity_Component(2), Velocity_Component(42)
       write(*, '("    VCN[",i4,"]: ", 6es20.12)') 771, Velocity_Component_Ngbr(:,771)
       do f = 1, nfc
          where(Boundary_Flag(f,:)<0)
             ! Interpolate at all interior faces
             Int_factor(:) = Face_Interpolation_Factor(f,:)
             Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + Cell(:)%Face_Normal(n,f) *   &
                                    ( (1-Int_factor(:))*Velocity_Component(:)  +       &
                                         Int_factor(:)   *Velocity_Component_Ngbr(f,:)      )
          elsewhere(Boundary_Flag(f,:)==0)
             ! Solid or Symmetry Boundaries
             Fluxing_Velocity(f,:) = 0
          elsewhere(Boundary_Flag(f,:)==1)
             ! Dirichlet Pressure Boundaries
             Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + Cell(:)%Face_Normal(n,f)*Velocity_Component(:)
          elsewhere(Boundary_Flag(f,:)==2)
             ! Dirichlet Velocity Boundaries
             !! NNC, Jan 2014.  Time-dependent dirichlet velocity
             !ORIG: Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + Cell(:)%Face_Normal(n,f)*BC_Vel(n,f,:)
             Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + Cell(:)%Face_Normal(n,f)*bndry_vel%get_all(n,f,t)
          endwhere
       enddo
    enddo
    NULLIFY(Velocity_Component_Ngbr)
    write(*, '("  Pre-SPECIAL CASES FV[",i4,"]: ", 6es20.12)') 771, Fluxing_Velocity(:,771)
    ! Correct Special Cases
    fluidRho_Ngbr => Fluxing_Velocity_Ngbr
    call EE_Gather(fluidRho_Ngbr, fluidRho)
    do f = 1,nfc
       do c = 1, ncells
          if(Solid_Face(f,c)) then
             !  Zero Velocity on Solid Faces
             Fluxing_Velocity(f,c) = 0
          elseif(fluidRho_Ngbr(f,c)==0 .AND. Mesh(c)%Ngbr_cell(f) /= 0 ) then
             ! Use cell center velocity if adjacent cell is void
             Fluxing_Velocity(f,c) = 0
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
       where (FluidRho == 0 .AND. Mesh%Ngbr_cell(f) /= 0 .AND. .not.IsPureImmobile) &
          Fluxing_Velocity(f,:) = Fluxing_Velocity_Ngbr(f,:)
    end do

    DEALLOCATE(Fluxing_Velocity_Ngbr)
    DEALLOCATE(Velocity_Component)
    DEALLOCATE(Int_factor)

  END SUBROUTINE INTERPOLATE_VELOCITY_TO_FACES

  SUBROUTINE MAC_CORRECTOR (Solution, Rho_Face, Fluxing_Velocity)
    !=======================================================================
    ! Purpose(s):
    !
    !   Apply a solenoidal correction, which is the gradient of the
    !   MAC solution field, to the face fluxing velocities.
    !=======================================================================
    use bc_module,              only: BC_Prs, bndry_vel !BC_Vel
    use do_interface,           only: DO_Specifier, do_init_ss, &
                                      do_gradient_face, DO_SOLVE_ORTHO, &
                                      DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: Solid_Face
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell
    use projection_data_module, only: dirichlet_pressure, Boundary_Flag, &
                                      dt_gradP_over_Rho
    use time_step_module,       only: t, dt

    ! Argument List
    real(r8), dimension(ncells),     intent(IN)    :: Solution
    real(r8), dimension(nfc,ncells), intent(IN)    :: Rho_Face
    real(r8), dimension(nfc,ncells), intent(INOUT) :: Fluxing_Velocity

    ! Local Variables
    type(DO_Specifier), pointer, save :: MAC_Cor_SS =>NULL()
    integer :: status
    integer :: f, n, MCSolveTech, i
    real(r8), dimension(:),     allocatable, save :: dt_over_Rho, Grad_Dot_N
    real(r8), dimension(:,:),   allocatable, save :: Grad
    real(r8), dimension(:,:,:), allocatable, save :: fgradx

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (dt_over_Rho(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_CORRECTOR: dt_over_Rho(ncells) allocation failed')
    ALLOCATE (Grad_Dot_N(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_CORRECTOR: Grad_Dot_N(ncells) allocation failed')
    ALLOCATE (Grad(ndim,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_CORRECTOR: Grad(ndim,ncells) allocation failed')
    ALLOCATE (fgradx(ndim,nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('MAC_CORRECTOR: fgradx(ndim,nfc,ncells) allocation failed')

    fgradx = 0

    if(.not.ASSOCIATED(MAC_Cor_SS))then
       MCSolveTech=DO_SOLVE_LU_LSLR; if(use_ortho_face_gradient)MCSolveTech=DO_SOLVE_ORTHO
       call do_init_ss(MAC_Cor_SS,SOLVETECH=MCSolveTech)
    endif

    call DO_GRADIENT_FACE(PHI=Solution, SOLVESPEC=MAC_Cor_SS, GRAD=fgradx)

    ! Loop over faces; compute the face potential gradient (dP_dx,dP_dy,dP_dz),
    ! then correct the face velocities with this gradient.

    FACE_LOOP: do f = 1,nfc

       do n = 1,ndim
          Grad(n,:) = fgradX(n,f,:)
       enddo

       ! Store the inverse face density.
       dt_over_Rho = 0
       where (Rho_Face(f,:) /= 0)  &
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
          where (Solid_Face(f,:)) Grad(n,:) = 0
       end do

       ! Get the normal component of the gradient.
       Grad_Dot_N = 0
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

       ! Now Apply BCs to the Fluxing_Velocity;

       ! Free-slip, Fluxing_Velocity = 0.
       where (Boundary_Flag(f,:) == 0 .or. Boundary_Flag(f,:) == 2) &
          Fluxing_Velocity(f,:) = 0

       ! Dirichlet velocity BCs:  Fluxing_Velocity = BC_Vel
       !! NNC, Jan 2014.  Time-dependent dirichlet velocity.
       !ORIG: do n = 1,ndim
       !ORIG:    where (Boundary_Flag(f,:) == 2) &
       !ORIG:       Fluxing_Velocity(f,:) = Fluxing_Velocity(f,:) + &
       !ORIG:                               BC_Vel(n,f,:)*Cell%Face_Normal(n,f)
       !ORIG: end do
       do i = 1, ncells
         if (Boundary_Flag(f,i) == 2) Fluxing_Velocity(f,i) = Fluxing_Velocity(f,i) + &
                                      dot_product(bndry_vel%get(f,i,t), Cell(i)%Face_Normal(:,f))
       end do

       ! Zero out Fluxing_Velocity's at solid faces.
       where (Solid_Face(f,:)) Fluxing_Velocity(f,:) = 0

    end do FACE_LOOP

    DEALLOCATE (dt_over_Rho)
    DEALLOCATE (Grad_Dot_N)
    DEALLOCATE (Grad)
    DEALLOCATE (fgradx)

  END SUBROUTINE MAC_CORRECTOR


  SUBROUTINE MAC_RHS (Fluxing_Velocity, RHS)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the RHS for the MAC projection, which is the
    !   cell-centered divergence of the face fluxing velocity.
    !=======================================================================
    use fluid_data_module,      only: Solid_Face
    use legacy_mesh_api,        only: ncells, nfc, Cell
    use projection_data_module, only: DVol_by_Dt_over_Vol

    ! Argument List
    real(r8), dimension(nfc,ncells), intent(IN)  :: Fluxing_Velocity
    real(r8), dimension(ncells),     intent(OUT) :: RHS

    ! Local Variables
    integer :: f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Compute the cell-centered divergence of the fluxing velocity as a
    ! discrete sum over faces (Green-Gauss).  Recall that the fluxing
    ! velocity has already been projected along the face normal. This will
    ! be the RHS for MAC projection.

    RHS = 0
    do f = 1, nfc
       where(.not.solid_face(f,:))
           RHS = RHS + Fluxing_Velocity(f,:)*Cell%Face_Area(f)
       endwhere
    end do

    !  Add in the volume change rate term due to the Enthalpy Solution
    RHS = RHS -  DVol_by_Dt_over_Vol*Cell%Volume
#ifndef NDEBUG
    write(*,*) "<<< Pre Boundary Conditions RHS"
    do f = 1, ncells
      write(*,'("Rhs[",i3, "]: ", es15.5)') f, RHS(f)
    end do
#endif

    write(*,'("Fluxing_Velocity[",i4,"]: ", 6es16.8)') 771, Fluxing_Velocity(:,771)
    write(*,'("MAC_RHS[",i4,"]: ", es16.8)') 771, RHS(771)
!!$    Don't add this in because we are solving for delta-P
!!$    RHS = RHS - Vol_over_RhoCsqDt*Zone%P

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
    use bc_module,              only: BC_Prs
    use cutoffs_module,         only: alittle
    !  FluidRho is needed for scaling of A_ortho and RHS
    use fluid_data_module,      only: FluidRho, isPureImmobile, MinFaceFraction, MinFluidRho
    use linear_solution,        only: LINEAR_SOLVER, Ubik_user
    use lnorm_module,           only: L1NORM
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell
    use preconditioners,        only: P, PRECONDITION
    use projection_data_module, only: Boundary_Flag, Coeff, UBIK_PRESSURE
    use y_eq_Ax_prs,            only: Y_EQ_AX_PRESSURE
    use time_step_module,       only: dt
    use UbikSolve_module
    use fischer_module

    ! Argument List
    character(len=*), intent(IN) :: timer_flag
    real(r8), dimension(nfc,ncells), intent(IN)    :: Rho_Face
    real(r8), dimension(ncells),     intent(INOUT) :: RHS
    real(r8), dimension(ncells),     intent(INOUT) :: Solution

    ! Local Variables
    integer  :: f, n, nc
    real(r8) :: X_Dot_N
    real(r8) :: RHSnorm
    real(r8), pointer :: A_ortho(:,:)

    ! Convergence criteria used with new PPE
    real(r8)   :: original_criterion, convergence_criterion

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! If preconditioning is on, get the orthogonal coefficients.
    if (Ubik_user(UBIK_PRESSURE)%precond /= 0) then

      ! Allocate the orthogonal coeffient array.
      allocate(A_Ortho(0:nfc,ncells))

      ! Point the preconditioning matrix to A_Ortho
      P => A_Ortho

      ! Set up orthogonal mesh 7-point stencil
      call ORTHO_STENCIL (RHO_FACE = Rho_Face, MATRIX = A_Ortho)

      ! Scale A_Ortho to be consistent with the RHS and mat-vec.
      ! Note 0 index for lower bound of A_Ortho
      do n=1,ncells
        if(FluidRho(n) == 0 .or. isPureImmobile(n)) then
          A_Ortho(0,n) = 1 / (MinFaceFraction * MinFluidRho * Cell(n)%Volume**0.66666)
        else
          do f=0,nfc
            A_Ortho(f,n) = A_Ortho(f,n)/Cell(n)%Volume
          enddo
        endif
      enddo

    end if

    ! Allocate coefficient array.
    allocate(Coeff(nfc,ncells))
    Coeff = 0

    ! Set Coeff = Area_Face/Rho_Face. If the face has Dirichlet BCs,
    ! then Coeff = -Area_Face*(X*n)/(Rho_Face*L**2), where X is a vector
    ! from the cell to face centroid, n is the unit normal, and L is
    ! the cell halfwidth.

    do nc = 1,ncells

      if(isPureImmobile(nc) .or. FluidRho(nc) == 0) cycle
      do f = 1,nfc

        ! Coeff = A_f/Rho_f.
        if (Rho_Face(f,nc) == 0) then
          Coeff(f,nc) = 0
          CYCLE
        else
          Coeff(f,nc) = Cell(nc)%Face_Area(f) / Rho_Face(f,nc)
        endif

        ! Dirichlet corrections.
        if (Boundary_Flag(f,nc) == 1) then
          X_Dot_N = 0
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

    ! Likely not the best place to put this, but the RHS must be divided by
    ! 'dt' so that we're solving for the pressure increment, and not pressure*dt.
    ! Scaling of the RHS here by dt and the cell volume is in conjunction with the
    ! scaled convergence criteria below.
    where(FluidRho(:) == 0 .or. isPureImmobile)
      RHS(:) = RHS(:) / (MinFaceFraction * MinFluidRho * Cell(:)%Volume**0.66666)
    elsewhere
      RHS(:) = RHS(:)/(dt*Cell(:)%Volume)
    endwhere
#ifndef NDEBUG
    write(*,*) "<<< RHS fed to Ubik"
    do f = 1, ncells
      write(*,'("Rhs[",i3, "]: ", es15.5)') f, RHS(f)
    end do
#endif

    write(*,'("PRESSURE RHS[",i4, "]: ", es15.5)') 771, RHS(771)
    ! Calculate scaling factor for L2 norm of the residual
    original_criterion = Ubik_eps(Ubik_user(UBIK_PRESSURE)%Control)
    convergence_criterion = original_criterion/(dt*dt)
    call Ubik_set_eps(Ubik_user(UBIK_PRESSURE)%Control, convergence_criterion)

    ! Solve for the pressure correction
    Solution = 0

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
    deallocate(Coeff)

    if (Ubik_user(UBIK_PRESSURE)%precond /= 0) then
      P => null()
      deallocate(A_Ortho)
    end if

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
    use cutoffs_module,         only: alittle
    use fluid_data_module,      only: FluidRho, Solid_Face, isPureImmobile
    use legacy_mesh_api,        only: ncells, nfc, ndim, Cell, Mesh, EE_GATHER
    use projection_data_module, only: dirichlet_pressure, Vol_over_RhoCsqDt
    use time_step_module,       only: dt

    ! Argument List
    real(r8), dimension(nfc,ncells),   intent(IN)    :: Rho_Face
    real(r8), dimension(0:nfc,ncells), intent(INOUT) :: Matrix
    real(r8), dimension(ncells), optional, intent(INOUT) :: RHS

    ! Local Variables
    integer :: j, f, n
    real(r8) :: distance, dX_tmp
    real(r8), dimension(nfc,ncells) :: Face_Coeff
    real(r8), dimension(ndim,nfc,ncells) :: dx_scaled, X_e
    real(r8), dimension(ndim) :: dx

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Face_Coeff = 0

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
         dX       = 0
         Distance = 0
         do n = 1,ndim
            dX(n)  = X_e(n,f,j) - Cell(j)%Centroid(n)
            Distance = Distance + dX(n)**2
         end do
         do n = 1,ndim
           dX_tmp = dX(n)/(Distance + alittle)
           dX_scaled(n,f,j) = dX_tmp
         end do
        ! Define face coefficients
        if (Rho_Face(f,j) == 0) then
           Face_Coeff(f,j) = 0
        else
           do n= 1,ndim
              Face_Coeff(f,j) = Face_Coeff(f,j) + Cell(j)%Face_Normal(n,f)*dX_scaled(n,f,j)
           enddo
           Face_Coeff(f,j) = Face_Coeff(f,j) * Cell(j)%Face_Area(f)/Rho_Face(f,j)
        end if
      end do FACE_LOOP
      end do CELL_LOOP


    ! Now load the coefficients at each face; check for BCs
    Matrix = 0
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
       if(FluidRho(n) == 0 .or. isPureImmobile(n)) then
          Matrix(0,n) = 1
          do f=1,nfc
             Matrix(f,n) = 0
          end do
       endif
    end do

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
    use do_interface,           only: DO_Specifier, do_init_ss, &
                                      DO_SOLVE_ORTHO, DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: isPureImmobile, fluidRho, &
                                      Centered_GradP_Dynamic,   &
                                      Momentum_by_Volume,       &
                                      fluidVof
    use legacy_mesh_api,        only: ndim, ncells
    use time_step_module,       only: dt
    use zone_module,            only: Zone
    use fluid_utilities_module, only: CC_GRADP_DYNAMIC

    ! Local Variables
    integer :: n
    integer :: PCSolveTech
    real(r8),  dimension(ncells) :: Kappa

    type(DO_Specifier), pointer, save :: PCorrector_SS =>NULL()

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
    !    Kappa = 1 - TI_Get_dV_by_Vo()
    !else
        Kappa = 1
    !endif

    ! Remove the old time dynamic pressure gradient acceleration from the
    ! cell centered velocity
    do n = 1,ndim
       Zone%Vc(n) = Zone%Vc(n) + dt * Centered_GradP_Dynamic(n,:)
       Momentum_by_Volume(n,:) = Momentum_by_Volume(n,:) + &
                      Centered_GradP_Dynamic(n,:)*dt*fluidRho*fluidVof/Kappa
       ! This should probably be here to ensure the zero's are consistent w. the
       ! loop following cc_gradp_dynamic() below.
       where (FluidRho == 0 .or. isPureImmobile)
          Zone%Vc(n) = 0
          Momentum_by_Volume(n,:) = 0
       endwhere
    end do

    call CC_GRADP_DYNAMIC()

    ! Apply the cell-centered gradient to Zone%Vc, and zero out values in
    ! void and solid cells.
    do n = 1,ndim
       Zone%Vc(n) = Zone%Vc(n) - Centered_GradP_Dynamic(n,:) * dt
       Momentum_by_Volume(n,:) = Momentum_by_Volume(n,:) - &
                     Centered_GradP_Dynamic(n,:) * dt * fluidRho*fluidVof/Kappa
       where (FluidRho == 0 .or. isPureImmobile)
                Zone%Vc(n) = 0
                Momentum_by_Volume(n,:) = 0
       endwhere
    end do

    ! Reset the pointer for pressures to BC_Zero for the pressure change solution
    BC_Prs => BC_Zero

    ! stop the projection_corrector timer.
    call stop_timer("timer_projection_corrector")

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
    use fluid_data_module, only: fluidRho
    use fluid_type_module, only: Div_c
    use lnorm_module,      only: L1NORM, L2NORM, LINORM
    use legacy_mesh_api,   only: ncells, nfc, Cell
    use pgslib_module,     only: PGSLIB_GLOBAL_MAXLOC
    use time_step_module,  only: dt

    ! Argument List
    real(r8), dimension(nfc,ncells), intent(IN) :: Fluxing_Velocity

    ! Local Variables
    integer :: f
    integer, dimension(1) :: MaxLoc_L
    real(r8), dimension(ncells) :: Velocity_Divergence

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Compute the velocity divergence.
    Velocity_Divergence = 0
    do f = 1, nfc
       Velocity_Divergence = Velocity_Divergence + &
                             Fluxing_Velocity(f,:)*Cell%Face_Area(f)
    end do

    Velocity_Divergence = (Velocity_Divergence/Cell%Volume)*dt

    ! Zero out divergences in void cells
    where (FluidRho == 0) Velocity_Divergence = 0

    !save as a residual

    ! Now get the norms; store them in Div.
    Div_c%V_f%L1   = L1NORM(Velocity_Divergence)
    Div_c%V_f%L2   = L2NORM(Velocity_Divergence)
    Div_c%V_f%Linf = LINORM(Velocity_Divergence)

    ! Get the location of the Linf norm.
    MaxLoc_L = PGSLib_Global_MAXLOC(ABS(Velocity_Divergence))
    Div_c%V_f%Linf_Location = MaxLoc_L(1)

  END SUBROUTINE VF_DIVERGENCE_NORMS

END MODULE PROJECTION_MODULE
