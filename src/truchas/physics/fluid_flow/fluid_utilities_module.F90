!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE FLUID_UTILITIES_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables needed to solve the Navier-Stokes equations.
  !
  ! Contains: FLUID_INIT
  !           RESTART_FLOW
  !           FLUID_DEALLOCATE
  !           FLUIDDENSITYFACE
  !
  ! Author(s): Jim Sicilian, LANL CCS-2 (sicilian@lanl.gov)   Jan 2003
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  public :: FLUID_INIT,         &
            FLUID_DEALLOCATE,   &
            FLUIDDENSITYFACE,   &
            CC_GRADP_DYNAMIC,   &
            calcVelLimits

CONTAINS

  SUBROUTINE FLUID_INIT (t)
    !=======================================================================
    ! Purpose:
    !
    ! Allocate and Initialize the volume fraction of each cell that holds a
    ! material that flows.  Also sets up Cell_isnt_Void and Ngbr_isnt_Void
    ! logical arrays.
    !
    ! Allocate and Initialize the storage for the Explicit Pressure Gradient
    ! and Buoyant force.
    !
    !=======================================================================

    use advection_data,         only: Momentum_Delta
    use bc_module,              only: BC, DIRICHLET, DIRICHLET_VEL,       &
                                      FREE_SLIP, INTERNAL_BC, Prs, Vel
    use cutoffs_module,         only: alittle
    use fluid_data_module,      only: fluidVof, fluidVof_n,               &
                                      realfluidVof, cutRho,               &
                                      Cell_isnt_void, Ngbr_isnt_void,     &
                                      Centered_gradP_Dynamic, Rho_Face_n, &
                                      fluid_flow, isPureImmobile,         &
                                      fluidRho, fluidRho_n,               &
                                      avgRho, avgRho_n,                   &
                                      fluidDeltaRho,                      &
                                      Solid_Face, fluid_to_move,          &
                                      Momentum_by_Volume,                 &
                                      Fluxing_Velocity,                   &
                                      Face_Interpolation_Factor,          &
                                      Drag_Coefficient,                   &
                                      Mom_Delta,                          &
                                      courant
    use fluid_type_module,      only: NORMS, DIV_NORMS, Div_c, Div_f
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell, EE_GATHER
    use projection_data_module, only: Boundary_Flag, DVol_by_Dt_over_Vol
    use property_module,        only: FLUID_PROPERTIES
    use restart_variables,      only: restart, have_fluid_flow_data
    use surface_tension_module, only: face_set_ids, csf_boundary, csf_z
    use legacy_mesh_api,        only: mesh_face_set

    real(r8), intent(in) :: t

    ! local variables
    integer :: n, status, f, j, m, nssets, fminz, fmaxz
    logical :: abort
    real(r8), dimension(:,:), allocatable :: Intercell_Distance_Sq, Ngbr_Centroid
    real(r8), dimension(:),   allocatable :: Cell_Centroid_C
    integer, dimension(ncells) :: csf_flag

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE(fluidVof(ncells),                      &
             fluidVof_n(ncells),                    &
             fluidRho(ncells),                      &
             fluidRho_n(ncells),                    &
             fluidDeltaRho(ncells),                 &
             realFluidVof(ncells),                  &
             cutRho(ncells),                        &
             avgRho(ncells),                        &
             avgRho_n(ncells),                      &
             isPureImmobile(ncells),                &
             Momentum_Delta(ndim,ncells),           &
             Mom_Delta(ndim,ncells),                &
             Cell_isnt_void(ncells),                &
             Ngbr_isnt_void(nfc,ncells),            &
             Solid_Face(nfc,ncells),                &
             Centered_GradP_Dynamic(ndim,ncells),   &
             Momentum_by_Volume(ndim,ncells),       &
             Rho_Face_n(nfc,ncells),                &
             Boundary_Flag(nfc,ncells),             &
             DVol_by_Dt_over_Vol(ncells),           &
             Cell_Centroid_C(ncells),               &
             Intercell_Distance_Sq(nfc,ncells),     &
             Ngbr_Centroid(nfc,ncells),             &
             Face_Interpolation_Factor(nfc,ncells), &
             Drag_Coefficient(ndim,ncells),         &
             courant(ncells),                       &
             csf_z(ncells),                         &
             STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_INIT: memory allocation failed')

    ! Precompute and store cell sizes in z direction for those volume cells that
    ! are adjacent to side set on which "csf_boundary"-type tagnential surface
    ! tension is to be applied
    if (csf_boundary) then
      ! Find volume cells adjacent to side set on which "csf_boundary"-type
      ! surface tension is to be applied (tag those with -1.0 as opposed to the
      ! initial 0.0 values)
      csf_flag = 0
      nssets = SIZE(Mesh_Face_Set,1);
      do j = 1, ncells
          do m = 1, nfc
            do n = 1,nssets
              if (ANY(face_set_ids == Mesh_Face_Set(n,m,j))) then
                csf_flag(j) = 1
              end if
            end do
          end do
      end do

      ! Compute and store volume mesh cell size in z direction for cells adjacent
      ! to side set on which "csf_boundary"-type surface tension is to be applied
      csf_z = -1.0
      do j = 1, ncells
        if (csf_flag(j) == 1) then
          ! find cell face ids with minimum and maximum z coordinates of the cell
          ! centroids
          fminz = minloc(Cell(j)%Face_Centroid(3,:), dim=1)
          fmaxz = maxloc(Cell(j)%Face_Centroid(3,:), dim=1)
          ! compute and store size of cell in z direction
          csf_z(j) = abs( Cell(j)%Face_Centroid(3,fmaxz) - &
                          Cell(j)%Face_Centroid(3,fminz) )
       end if
      end do
    end if

    call FLUID_PROPERTIES (abort, t)
    if(abort) then
        fluid_to_move = .false.
    else
        fluid_to_move = fluid_flow
    endif

    ! Initialize the average density for body forces
    avgRho_n = fluidVof*(fluidRho+fluidDeltaRho)

    ! Initialise courant number array
    courant = 0.0_r8

    ! Give initial values to Div_c, and Div_f (they will be printed before computed :-/
    Div_c = DIV_NORMS(NORMS(0.0, 0.0, 0.0, 0))
    Div_f = DIV_NORMS(NORMS(0.0, 0.0, 0.0, 0))

    ! Compute the face boundary flag
    Boundary_Flag = -1
    do f = 1,nfc
       ! Boundary_Flag = 0 at solid boundaries,
       where (FREE_SLIP (BC%Flag, Vel%Face_bit(f)) .or.   &
              INTERNAL_BC (BC%Internal, Vel%Face_bit(f))) &
          Boundary_Flag(f,:) = 0

       ! Boundary_Flag = 1 where Dirichlet pressure BCs apply
       where (DIRICHLET (BC%Flag, Prs%Face_bit(f)))       &
          Boundary_Flag(f,:) = 1

       ! Boundary_Flag = 2 where Dirichlet velocity BCs apply
       where (DIRICHLET_VEL (BC%Flag, Vel%Face_bit(f)))   &
          Boundary_Flag(f,:) = 2
    end do

    ! Define the interpolation factors for transferring cell centered data to faces
    !     (used in INTERPOLATE_TO_FACES)
    Face_Interpolation_Factor = 0.0_r8
    Intercell_Distance_Sq     = 0.0_r8
    do n = 1, ndim
        Cell_Centroid_C(:) = Cell(:)%Centroid(n)
        call EE_GATHER(Ngbr_Centroid, Cell_Centroid_C)
        do f = 1, nfc
            where(Boundary_Flag(f,:) < 0)
              Intercell_Distance_Sq(f,:) = Intercell_Distance_Sq(f,:) + (Ngbr_Centroid(f,:)-Cell_Centroid_C(:))**2
              Face_Interpolation_Factor(f,:) = Face_Interpolation_Factor(f,:) +              &
                                            (Cell(:)%Face_Centroid(n,f)-Cell_Centroid_C(:))*  &
                                            (Ngbr_Centroid(f,:)-Cell_Centroid_C(:))
            endwhere
        enddo
    enddo
    do f = 1, nfc
        where(Intercell_Distance_Sq(f,:)>alittle)   &
            Face_Interpolation_Factor(f,:) = Face_Interpolation_Factor(f,:) / Intercell_Distance_Sq(f,:)
    enddo

    if(fluid_flow .and. restart .and. have_fluid_flow_data) then
       call RESTART_FLOW
    else
       Fluxing_Velocity = 0.0_r8
       Centered_GradP_Dynamic = 0.0_r8
       Rho_Face_n = 0.0_r8
    endif

    fluidRho_n = fluidRho
    fluidVof_n = fluidVof

    !  Initialize Momentum Advection
    Momentum_Delta = 0.0_r8

    !  Initialize volume Change Rate Array
    DVol_by_Dt_over_Vol = 0.0_r8

    !  Initialize old time momentum storage
    Momentum_by_Volume = 0.0_r8

    DEALLOCATE (fluidDeltaRho)
    DEALLOCATE (isPureImmobile)
    DEALLOCATE (Solid_Face)
    DEALLOCATE (Intercell_Distance_Sq)
    DEALLOCATE (Cell_Centroid_C)
    DEALLOCATE (Ngbr_Centroid)

  END SUBROUTINE FLUID_INIT

  SUBROUTINE FLUID_DEALLOCATE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Deallocate the fluidvof, Cell_isnt_Void, Ngbr_isnt_Void and
    !   Centered_GradP_Dynamic  arrays
    !
    !=======================================================================
    use advection_data,         only: Momentum_Delta
    use fluid_data_module,      only: fluidvof, fluidVof_n, Cell_isnt_void, Ngbr_isnt_Void, &
                                      Centered_gradP_Dynamic, Rho_Face_n, fluidRho_n, &
                                      fluidrho, Momentum_by_Volume, Face_Interpolation_Factor, &
                                      Drag_Coefficient, Mom_Delta, courant
    use projection_data_module, only: Boundary_Flag, DVol_by_Dt_over_Vol
    use surface_tension_module, only: csf_z

    ! Local Variables
    integer :: memstat

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! deallocate the Face_Interpolation_Factor array.
    if (ASSOCIATED(Face_Interpolation_Factor)) DEALLOCATE (Face_Interpolation_Factor)

    ! deallocate the fluidvof array.
    if (ASSOCIATED(fluidvof)) DEALLOCATE (fluidvof)

    ! deallocate the fluidvof_n array.
    if (ASSOCIATED(fluidvof_n)) DEALLOCATE (fluidvof_n)

    ! deallocate the fluidRho array
    if (ASSOCIATED(fluidRho)) DEALLOCATE (fluidRho)

    ! deallocate the fluidrho_n array
    if (ASSOCIATED(fluidRho_n)) DEALLOCATE (fluidRho_n)

    ! deallocate the Momentum_Delta array
    if (allocated(Momentum_Delta)) deallocate (Momentum_Delta)

    ! deallocate the Mom_Delta array
    if (allocated(Mom_Delta)) deallocate(Mom_Delta)

    ! deallocate the void cell indication arrays
    if (ASSOCIATED(Cell_isnt_Void)) DEALLOCATE (Cell_isnt_Void)

    ! deallocate the Ngbr_isnt_Void array
    if (ASSOCIATED(Ngbr_isnt_Void)) DEALLOCATE (Ngbr_isnt_Void)

    ! deallocate the Projection / Predictor communiation arrays
    if (allocated(Centered_gradP_Dynamic)) deallocate (Centered_GradP_Dynamic)

    if (allocated(Rho_Face_n)) deallocate (Rho_Face_n)

    if (associated(boundary_flag)) deallocate (Boundary_Flag)

    if (ALLOCATED(DVol_by_Dt_over_Vol)) DEALLOCATE(DVol_by_Dt_over_Vol)
    if (ASSOCIATED(Momentum_by_Volume)) DEALLOCATE(Momentum_by_Volume)
    if (ALLOCATED(Drag_Coefficient)) DEALLOCATE(Drag_Coefficient)
    if (ASSOCIATED(courant)) DEALLOCATE (courant)
    if (ALLOCATED(csf_z)) DEALLOCATE (csf_z)

  END SUBROUTINE FLUID_DEALLOCATE

  SUBROUTINE RESTART_FLOW()
    !=======================================================================
    ! Purpose:
    !
    !  reset saved flow solution quantities at the start of restart calculations
    !
    !    Jim Sicilian, LANL CCS-2 (sicilian@lanl.gov)   Dec 2002
    !
    !=======================================================================
    use bc_module,             only: BC_Prs
    use bc_data_module,        only: BC_Pressure, BC_Zero
    use body_data_module,      only: body_force_face_method
    use discrete_ops_data,     only: use_ortho_face_gradient
    use fluid_data_module,     only: fluidRho, fluidDeltaRho,    &
                                     fluidvof, avgRho_n,           &
                                     Rho_Face, Rho_Face_n
    use legacy_mesh_api,       only: ncells, ndim, nfc

    use projection_data_module, only: dt_gradP_over_Rho, ghc, ghn,dtCsf_over_Rho, &
                                      dtRhoG_over_Rho
    use body_force_module,      only: add_face_body_force, compute_gravityhead

    use time_step_module,      only: dt
    use cutoffs_module,        only: alittle
    use surface_tension_module, only: CSF_FACE, csf_normal

    ! Local Variables
    integer       :: n, f, status

    logical,  dimension(:),   allocatable       :: Mask
    real(r8), dimension(:),   allocatable       :: Scalar_Cell_Center
    real(r8), dimension(:),   allocatable       :: dt_over_Rho, Grad_Dot_N
    real(r8), dimension(:,:), allocatable       :: Grad
    real(r8), dimension(:,:,:), allocatable     :: Gradient
!-mf Jan04
    real(r8), dimension(:,:,:), allocatable     :: Fcsf

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (Rho_Face(nfc,ncells),               &
              dt_gradP_over_Rho(ndim,nfc,ncells), &
              dtRhoG_over_Rho(ndim,nfc,ncells),   &
              STAT = status)
    if (status /= 0) call TLS_panic ('RESTART_FLOW: allocation failed')

    if (use_ortho_face_gradient) then
      ALLOCATE (ghc(nfc,ncells), &
           ghn(nfc,ncells), STAT = status)
      if (status /= 0) call TLS_panic ('RESTART_FLOW: allocation failed')
    end if

    if (csf_normal) then
      ALLOCATE (dtCsf_over_Rho(ndim,nfc,ncells), STAT = status)
      if (status /= 0) call TLS_panic ('RESTART_FLOW: allocation failed')
    endif

    ALLOCATE (Mask(ncells), &
         Scalar_Cell_Center(ncells), &
         Gradient(ndim,nfc,ncells), &
         Grad(ndim,ncells), &
         Grad_Dot_N(ncells), &
         dt_over_Rho(ncells), STAT = status)

    ! Setup the averged density for the body forces
    avgRho_n = fluidVof*(fluidRho+fluidDeltaRho)

    ! Setup the face densities
    call FluidDensityFace

    Rho_Face_n = Rho_Face


    ! Use the real Pressures at Dirichlet boundaries in this pressure gradient calculation
    BC_Prs => BC_Pressure

    ! calculate gravity head ghc and ghn...
    if (use_ortho_face_gradient) then
       if (body_force_face_method) then
          call COMPUTE_GRAVITYHEAD()
       else
          ghc = 0.0
          ghn = 0.0
       end if
    end if

    !-mf Jan04
    if (csf_normal) then
      ALLOCATE (Fcsf(ndim,nfc,ncells), STAT = status)
      if (status /= 0) call TLS_panic ('RESTART_FLOW: Fcsf(ndim,nfc,ncells) allocation failed')
      call CSF_FACE(Fcsf)

      ! dtCsf_over_Rho for use in CC_GRADP_DYNAMIC...

      do f=1,nfc
        do n=1,ndim
          !-particular case when void cells. density = zero
          where (Rho_Face(f,:) <= alittle)
            dtCsf_over_Rho(n,f,:)=0.0_r8
          elsewhere
            dtCsf_over_Rho(n,f,:)=dt*Fcsf(n,f,:)/Rho_Face(f,:)
          endwhere
        enddo
      enddo

    endif

    if (.not. use_ortho_face_gradient) then
      ! dtRhoG_over_Rho for use in CC_GRADP_DYNAMIC...
      call ADD_FACE_BODY_FORCE(dt, dtRhoG_over_Rho)
    end if

    call CC_GRADP_DYNAMIC ()

    ! Reset the pointer for pressures to BC_Zero for the pressure change solution
    BC_Prs => BC_Zero

    DEALLOCATE(Rho_Face, dt_gradP_over_Rho, dtRhoG_over_Rho)

    if (use_ortho_face_gradient) DEALLOCATE(ghc, ghn)

    if (csf_normal) DEALLOCATE(dtCsf_over_Rho)

    DEALLOCATE(Mask, Scalar_Cell_Center, Gradient, Grad, Grad_Dot_N, dt_over_Rho)

    END SUBROUTINE RESTART_FLOW

    SUBROUTINE FLUIDDENSITYFACE
    !=======================================================================
    ! Purpose(s):
    !
    !  Evaluate and store face values of the fluid density for use throughout
    !  the Projection algorithm.
    !
    !    Jim Sicilian (CCS-2)  October, 2002
    !
    !=======================================================================
    use fluid_data_module
    use legacy_mesh_api, only: ncells, nfc, Mesh, Cell, DEGENERATE_FACE, EE_GATHER

    ! Local Variables
    integer :: status
    integer :: f, i
    real(r8) :: weight, weight_ngbr
    real(r8), dimension(:),   allocatable :: CellVolume
    real(r8), dimension(:,:), allocatable :: fluidrho_Ngbr, Volume_Ngbr, fluidvof_Ngbr

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

   ALLOCATE (fluidRho_Ngbr(nfc,ncells), STAT = status)
   if (status /= 0) call TLS_panic ('FluidDensityFace: fluidRho_Ngbr(nfc,ncells) allocation failed')
   ALLOCATE (Volume_Ngbr(nfc,ncells), STAT = status)
   if (status /= 0) call TLS_panic ('FluidDensityFace: Volume_Ngbr(nfc,ncells) allocation failed')
   ALLOCATE (fluidvof_Ngbr(nfc,ncells), STAT = status)
   if (status /= 0) call TLS_panic ('FluidDensityFace: fluidvof_Ngbr(nfc,ncells) allocation failed')
   ALLOCATE (CellVolume(ncells), STAT = status)
   if (status /= 0) call TLS_panic ('FluidDensityFace: CellVolume(ncells) allocation failed')

   ! Evaluate Weighted Face Densities

    ! Move Cell%Volume to a simple array for EE_GATHER
    CellVolume = Cell%Volume

    call EE_GATHER (fluidRho_Ngbr, fluidRho)
    call EE_GATHER (Volume_Ngbr, CellVolume)
    call EE_GATHER (fluidvof_Ngbr, fluidvof)

    do i = 1,ncells
       do f = 1,nfc
          if (Mesh(i)%Ngbr_Cell(f) == DEGENERATE_FACE) then
             rho_face(f,i) = fluidRho(i)
          else if (Solid_Face(f,i)) then
             rho_face(f,i) = fluidRho(i)
          else if (Mesh(i)%Ngbr_Cell(f) == 0) then
             rho_face(f,i) = fluidRho(i)
          else
             weight_ngbr   = Volume_Ngbr(f,i)*fluidvof_Ngbr(f,i)
             weight        = CellVolume(i)*fluidvof(i)
             rho_face(f,i) = (fluidrho(i)*weight + fluidrho_Ngbr(f,i)*weight_ngbr)/    &
                             (weight + weight_ngbr)
          end if
          if(rho_face(f,i) > 0.0_r8) rho_face(f,i) = MAX(MinFaceFraction*MinFluidRho, rho_face(f,i))
       end do
   end do

    DEALLOCATE (fluidRho_ngbr)
    DEALLOCATE (Volume_Ngbr)
    DEALLOCATE (fluidvof_Ngbr)
    DEALLOCATE (CellVolume)

  END SUBROUTINE FLUIDDENSITYFACE

  SUBROUTINE CC_GRADP_DYNAMIC ()
    !=======================================================================
    ! Purpose(s):
    !
    !   The last step in the PROJECTION call, this routine calculates a
    !   time 'n+1' cell-centered pressure gradient, and applies it to the
    !   cell-centered velocities.
    !=======================================================================
    use bc_module,              only: BC, DIRICHLET, Prs, BC_Prs
    use bc_data_module,         only: BC_Pressure, BC_Zero
    use bc_operations,          only: Pressure_BC
    use do_interface,           only: DO_Specifier, do_init_ss, &
                                      do_gradient_face, DO_SOLVE_ORTHO, &
                                      DO_SOLVE_LU_LSLR
    use discrete_ops_data,      only: use_ortho_face_gradient
    use fluid_data_module,      only: Solid_Face, Centered_GradP_Dynamic, Rho_Face
    use legacy_mesh_api,        only: ncells, ndim, nfc, Cell, Mesh
    use projection_data_module, only: Boundary_Flag, dirichlet_pressure, &
                                      dt_gradP_over_Rho, ghc, ghn,dtCsf_over_Rho, &
                                      dtRhoG_over_Rho
    use time_step_module,       only: dt
    use zone_module,            only: Zone
    use discrete_op_module,     only: DYNAMIC_PRESSURE_FACE_GRADIENT
!-mf Jan04
    use surface_tension_module, only: csf_normal

    ! Local Variables                            :: i
    integer :: status
    integer :: f, n
    integer :: PCSolveTech

    type(DO_Specifier), pointer, save :: PCorrector_SS =>NULL()

    logical,  dimension(:), allocatable :: Mask
    real(r8), dimension(:), allocatable :: cell_vol, Weight, Scalar_Cell_Center, Grad_Dot_N
    real(r8), dimension(:,:), allocatable   :: Grad
    real(r8), dimension(:,:,:), allocatable :: Gradient

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the projection_corrector timer.
    call start_timer("timer_projection_corrector")

    ALLOCATE (Mask(ncells), &
              cell_vol(ncells), &
              Grad(ndim,ncells), &
              Weight(ncells), &
              Scalar_Cell_Center(ncells), &
              Gradient(ndim,nfc,ncells), &
              Grad_Dot_N(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('CC_GRADP_DYNAMIC: allocation failed')

    if(.not.ASSOCIATED(PCorrector_SS))then
       PCSolveTech=DO_SOLVE_LU_LSLR; if(use_ortho_face_gradient)PCSolveTech=DO_SOLVE_ORTHO
       call do_init_ss(PCorrector_SS,SOLVETECH=PCSolveTech,BC_SPEC=Pressure_BC)
    endif

    ! Use the real Pressures at Dirichlet boundaries in this pressure gradient calculation
    BC_Prs => BC_Pressure

    ! Now add in the explicit pressure gradient and buoyant term at this face
    Scalar_Cell_Center = Zone%P
    if (.not. use_ortho_face_gradient) then
       ! calculate the gradient of the total pressure.
       call DO_GRADIENT_FACE(PHI=Scalar_Cell_Center, SOLVESPEC=PCorrector_SS, GRAD=Gradient)
    else
       ! calculate the gradient of the dynamic pressure
       ! ie; the total pressure less the gravity  head.
       call DYNAMIC_PRESSURE_FACE_GRADIENT(Gradient, Scalar_Cell_Center, ghc, ghn)
    end if

    dt_gradP_over_Rho = 0.0_r8

    FACE_LOOP : do f = 1,nfc

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

       ! Kill the gradient at solid faces.
       do n = 1,ndim
          where (Solid_Face(f,:)) Grad(n,:) = 0.0_r8
       end do

       ! Get the normal component of the gradient.
       Grad_Dot_N = 0.0_r8
       do n = 1,ndim
          Grad_Dot_N = Grad_Dot_N + Grad(n,:)*Cell%Face_Normal(n,f)
       end do


       ! Correct the face velocities.
       do n = 1,ndim
          ! Kill the normal component of the gradient where appropriate.
          where (Boundary_Flag(f,:) == 0 .or. Boundary_Flag(f,:) == 2) &
               Grad(n,:) = Grad(n,:) - Grad_Dot_N*Cell%Face_Normal(n,f)
       end do

       ! Store the inverse face density.
       do n = 1,ndim
          where (Rho_Face(f,:) /= 0.0_r8) &
               dt_gradP_over_Rho(n,f,:) = dt*Grad(n,:)/Rho_Face(f,:)
       enddo

    enddo FACE_LOOP

    Centered_GradP_Dynamic = 0.0_r8
    Grad = 0.0_r8

    if (.not. use_ortho_face_gradient) then
      ! we are working with the total pressure...
      do n = 1,ndim
         cell_vol = 0.0_r8
         do f = 1,nfc
            Mask = (Mesh%Ngbr_cell(f) /= 0 .OR. DIRICHLET (BC%Flag, Prs%Face_bit(f))) &
                    .AND. .NOT. Solid_Face(f,:)
            where (Mask)
               cell_vol = cell_vol + ABS(Cell%Face_Area(f)*Cell%Face_normal(n,f))
               Grad(n,:) = Grad(n,:) + (dt_gradP_over_Rho(n,f,:) - dtRhoG_over_Rho(n,f,:))&
                                     * ABS(Cell%Face_Area(f)*Cell%Face_normal(n,f))
            end where
            !-mf Jan04
            if (csf_normal) then
              where (Mask)
                 Grad(n,:) = Grad(n,:) - dtCsf_over_Rho(n,f,:) &
                                       * ABS(Cell%Face_Area(f)*Cell%Face_normal(n,f))
              end where
            endif

         end do
         where (cell_vol > 0.0_r8)
            Grad(n,:) = Grad(n,:)/cell_vol
            Centered_GradP_Dynamic(n,:) = Centered_GradP_Dynamic(n,:) + Grad(n,:) / dt
         endwhere
       end do
    else
      ! we are working with the dynamic pressure...
      do n = 1,ndim
         cell_vol = 0.0_r8
         do f = 1,nfc
            Mask = (Mesh%Ngbr_cell(f) /= 0 .OR. DIRICHLET (BC%Flag, Prs%Face_bit(f))) &
                    .AND. .NOT. Solid_Face(f,:)
            where (Mask)
               cell_vol = cell_vol + ABS(Cell%Face_Area(f)*Cell%Face_normal(n,f))
               Grad(n,:) = Grad(n,:) + (dt_gradP_over_Rho(n,f,:))&
                                     * ABS(Cell%Face_Area(f)*Cell%Face_normal(n,f))
            end where
            !-mf Jan04
            if (csf_normal) then
              where (Mask)
                 Grad(n,:) = Grad(n,:) - dtCsf_over_Rho(n,f,:)&
                                       * ABS(Cell%Face_Area(f)*Cell%Face_normal(n,f))
              end where
            endif
         end do
         where (cell_vol > 0.0_r8)
             Grad(n,:) = Grad(n,:)/cell_vol
             Centered_GradP_Dynamic(n,:) = Centered_GradP_Dynamic(n,:) + Grad(n,:) / dt
          endwhere
       end do
    end if

    write(*,*) "<< CC_GRADP_DYNAMIC"
    write(*, '("P[",i3,"]: ", es20.12)') 771, Zone(771)%P
    write(*, '("GradP_rho_cc[",i3,"]: ", 3es20.12)') 771, Centered_GradP_Dynamic(:,771)

    ! Reset the pointer for pressures to BC_Zero for the pressure change solution
    BC_Prs => BC_Zero

    DEALLOCATE(Mask, cell_vol, Grad, Weight, Scalar_Cell_Center, Gradient, Grad_Dot_N)

    ! Stop the projection_corrector timer.
    call stop_timer("timer_projection_corrector")

  END SUBROUTINE CC_GRADP_DYNAMIC

  subroutine calcVelLimits()
!===============================================================================
! Purpose:
! Compute the min/max velocities to report to the screen
!
! Author(s): M. A. Christon, LANL CCS-2 (christon@lanl.gov)
!
!===============================================================================
    use legacy_mesh_api,     only: ndim, ncells
    use fluid_data_module,   only: minVel, maxVel
    use pgslib_module,       only: pgslib_global_minval, pgslib_global_maxval
    use zone_module,         only: Zone

    ! Local Variables
    integer :: i, j

    minVel(:) =  huge(0.0_r8)
    maxVel(:) = -huge(0.0_r8)
    do i = 1, ncells
       do j = 1, ndim
          if( Zone(i)%Vc(j) .le. minVel(j)) then
             minVel(j) = Zone(i)%Vc(j)
          endif
          if( Zone(i)%Vc(j) .ge. maxVel(j)) then
             maxVel(j) = Zone(i)%Vc(j)
          endif
       end do
    end do

    ! Handle the min/max values in parallel
    do j = 1, ndim
      minVel(j) = pgslib_global_minval(minvel(j))
      maxVel(j) = pgslib_global_maxval(maxvel(j))
   end do

  end subroutine calcVelLimits

END MODULE FLUID_UTILITIES_MODULE
